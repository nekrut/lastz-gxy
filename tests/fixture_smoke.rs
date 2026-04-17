//! End-to-end smoke test: runs the library driver on a tiny fixture and
//! checks that the expected single HSP is emitted. Guards against
//! regressions in the Phase 1 pipeline (sequences → pos_table →
//! seed_search → driver → output).

use lastz_gxy::driver::{run, Config, StrandSpec};
use lastz_gxy::gapped_extend::GappedParams;
use lastz_gxy::hsp::HspParams;
use lastz_gxy::output::Strand;
use lastz_gxy::seeds::SeedPattern;
use lastz_gxy::sequences::{PackedSeq, Sequence};

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
