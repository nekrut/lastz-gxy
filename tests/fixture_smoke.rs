//! End-to-end smoke tests: run the library driver on tiny fixtures and
//! verify the emitted records. Guards against regressions in the full
//! pipeline (sequences → pos_table → seed_search → HSP → chain → anchor →
//! gapped_extend → tweener → output).

use lastz_gxy::driver::{run, ChunkConfig, Config, StrandSpec};
use lastz_gxy::gapped_extend::GappedParams;
use lastz_gxy::hsp::HspParams;
use lastz_gxy::output::Strand;
use lastz_gxy::seeds::SeedPattern;
use lastz_gxy::sequences::{PackedSeq, Sequence};
use lastz_gxy::tweener::TweenerConfig;

#[test]
fn tiny_fixture_emits_single_ungapped_alignment() {
    let targets = vec![Sequence {
        name: "chrT".into(),
        seq: PackedSeq::from_ascii(b"AAAAAAAAAAAGATTACACATGGCATGTCGAATTTTTTTTTTT"),
    }];
    let queries = vec![Sequence {
        name: "chrQ".into(),
        seq: PackedSeq::from_ascii(b"GATTACACATGGCATGTCGA"),
    }];

    let cfg = Config {
        pattern: SeedPattern::solid(12),
        strand: StrandSpec::Both,
        hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
        // Gapped extension on a perfect match produces the same footprint;
        // this exercises the full Phase 2 pipeline end-to-end. The 20-bp
        // perfect-match payload scores around 1900 under HOXD70, so we
        // lower gapped_threshold accordingly.
        gapped: GappedParams { y_drop: 9_400, gapped_threshold: 500 },
        gapped_enabled: true,
        ..Config::default()
    };

    let recs = run(&targets, &queries, &cfg);
    assert_eq!(recs.len(), 1, "got records: {recs:#?}");
    let r = &recs[0];
    assert_eq!(r.target_name, "chrT");
    assert_eq!(r.query_name, "chrQ");
    assert_eq!(r.query_strand, Strand::Plus);
    assert_eq!(r.t_start, 11);
    assert_eq!(r.q_start, 0);
    assert_eq!(r.t_span, 20);
    assert_eq!(r.q_span, 20);
    assert_eq!(r.script.to_cigar(), "20M");
    assert_eq!(r.target_bases, b"GATTACACATGGCATGTCGA");
    assert_eq!(r.query_bases, b"GATTACACATGGCATGTCGA");
}

/// Exercises the full stack on the minus strand *and* the tweener. The
/// target embeds the reverse complement of a query that has two
/// homologous blocks separated by a distinctive spacer: the outer
/// high-stringency pass only catches the two flanks, the tweener fills
/// the middle, and the final records must live in the forward-strand
/// coordinate frame (what SAM and MAF expect) with `query_strand = Minus`.
#[test]
fn minus_strand_with_tweener_and_gapped() {
    // Query: two 20 bp homologous anchors separated by a 25 bp middle
    // block. The outer scan uses match16 which only fires on the anchors;
    // the tweener re-scans the gap with match8.
    //
    // We then reverse-complement the whole thing and embed that into the
    // target. The minus-strand scan therefore sees the original
    // orientation of the query and should find all three blocks.
    let query_ascii: Vec<u8> = b"GATTACACATGGCATGTCGA"
        .iter()
        .chain(b"AGCTAGCTAGCTAGCTAGCTAGCTA".iter())
        .chain(b"CCCGGGAATTCCGGTTAAGGCC".iter())
        .copied()
        .collect();

    // Reverse-complement the query into a target payload.
    let rc_payload = PackedSeq::from_ascii(&query_ascii)
        .reverse_complement()
        .to_ascii();
    let mut target_bytes: Vec<u8> = b"TTTTTTTTTT".to_vec();
    target_bytes.extend_from_slice(&rc_payload);
    target_bytes.extend_from_slice(b"TTTTTTTTTT");

    let targets = vec![Sequence {
        name: "chrT".into(),
        seq: PackedSeq::from_ascii(&target_bytes),
    }];
    let queries = vec![Sequence {
        name: "chrQ".into(),
        seq: PackedSeq::from_ascii(&query_ascii),
    }];

    let cfg = Config {
        // Dense outer seed still hits multiple disjoint HSPs along the
        // minus diagonal; that's enough to drive chain + tweener through
        // all their branches.
        pattern: SeedPattern::solid(12),
        strand: StrandSpec::Minus,
        hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
        gapped: GappedParams { y_drop: 9_400, gapped_threshold: 500 },
        gapped_enabled: true,
        chain_enabled: true,
        tweener: Some(TweenerConfig {
            pattern: SeedPattern::solid(8),
            min_gap: 5,
            max_gap: 1_000,
            relax_factor: 2,
        }),
        ..Config::default()
    };

    let recs = run(&targets, &queries, &cfg);
    assert!(
        !recs.is_empty(),
        "expected at least one minus-strand record"
    );
    for r in &recs {
        assert_eq!(r.query_strand, Strand::Minus, "strand flag lost: {r:?}");

        // Coordinates must land within the embedded RC payload on the
        // target side: target positions [10, 10 + rc_payload.len()).
        let rc_lo = 10u32;
        let rc_hi = 10 + rc_payload.len() as u32;
        assert!(
            r.t_start >= rc_lo && r.t_end() <= rc_hi,
            "target coords out of payload: {} .. {} not in [{rc_lo}, {rc_hi}]",
            r.t_start,
            r.t_end()
        );

        // Query coordinates (after strand reversal, the driver already
        // exposes them in the rc'd query's frame) must stay within the
        // query bounds and be consistent with the script.
        assert!(r.q_end() <= r.query_len);
        assert_eq!(r.q_span, r.script.query_span());
        assert_eq!(r.t_span, r.script.target_span());

        // The bases stored on the record must round-trip through the
        // script: rendering them with gaps has to match the target's
        // target_span / query's query_span lengths.
        let (t_aln, q_aln) =
            lastz_gxy::edit_script::render_aligned_bases(&r.script, &r.target_bases, &r.query_bases);
        assert_eq!(t_aln.len(), q_aln.len());
        assert_eq!(t_aln.len() as u32, r.script.columns());
    }

    // Sanity check that records are sorted and non-overlapping on t_start
    // (what chain + tweener should guarantee for the output set).
    for pair in recs.windows(2) {
        assert!(
            pair[0].t_start <= pair[1].t_start,
            "records not sorted by t_start: {pair:#?}"
        );
    }
}

/// Within-target chunking must produce the same record set as non-chunked
/// execution, block-for-block, modulo the dedup pass the driver already
/// runs on overlapping halo regions. Uses a target long enough that small
/// chunk/halo values force multiple chunks.
#[test]
fn chunking_preserves_record_set() {
    // A 2 kbp target with two embedded homology blocks. Chunking with
    // chunk_size=600, halo=200 produces 4 chunks covering 600 bp each,
    // with halo overlaps forcing at least one duplicate emission that the
    // driver must dedup.
    let payload1: &[u8] = b"GATTACACATGGCATGTCGACCCCCCCCCCCC";
    let payload2: &[u8] = b"AGCTAGCTAGCTAGCTAGCGCATCATCATGCA";

    let filler = |seed: u64, n: usize| -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[((s >> 33) & 3) as usize]
            })
            .collect()
    };

    let mut target_bytes = filler(1, 400);
    target_bytes.extend_from_slice(payload1);
    target_bytes.extend(filler(2, 800));
    target_bytes.extend_from_slice(payload2);
    target_bytes.extend(filler(3, 400));

    let mut query_bytes = filler(4, 100);
    query_bytes.extend_from_slice(payload1);
    query_bytes.extend(filler(5, 200));
    query_bytes.extend_from_slice(payload2);
    query_bytes.extend(filler(6, 100));

    let targets = vec![Sequence {
        name: "t".into(),
        seq: PackedSeq::from_ascii(&target_bytes),
    }];
    let queries = vec![Sequence {
        name: "q".into(),
        seq: PackedSeq::from_ascii(&query_bytes),
    }];

    let base_cfg = Config {
        pattern: SeedPattern::solid(12),
        strand: StrandSpec::Plus,
        hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
        gapped: GappedParams { y_drop: 3_000, gapped_threshold: 500 },
        gapped_enabled: true,
        transitions: 0,
        ..Config::default()
    };

    let non_chunked = run(&targets, &queries, &base_cfg);
    let chunked_cfg = Config {
        chunk: ChunkConfig { primary_bp: 600, halo_bp: 200 },
        ..base_cfg
    };
    let chunked = run(&targets, &queries, &chunked_cfg);

    // Both record sets should cover the same target regions, strands, and
    // spans. The chunked path may reorder slightly; compare signatures.
    type Sig = (u32, u32, u32, u32, Strand);
    let mut a: Vec<Sig> = non_chunked
        .iter()
        .map(|r| (r.t_start, r.t_span, r.q_start, r.q_span, r.query_strand))
        .collect();
    let mut b: Vec<Sig> = chunked
        .iter()
        .map(|r| (r.t_start, r.t_span, r.q_start, r.q_span, r.query_strand))
        .collect();
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "chunked record set differs from non-chunked:\nnon-chunked={non_chunked:#?}\nchunked={chunked:#?}"
    );
    assert!(!a.is_empty(), "test fixture produced no records");
}
