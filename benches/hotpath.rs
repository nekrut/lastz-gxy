//! Microbenchmarks for the hot path stages.
//!
//! Establishes a reproducible baseline so future SIMD and GPU work can
//! report speedups against a fixed reference. Each bench uses a
//! deterministic synthetic input so numbers are comparable across machines.

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use lastz_gxy::anchor::choose_anchor;
use lastz_gxy::chain::best_chain;
use lastz_gxy::edit_script::{EditOp, EditScript};
use lastz_gxy::gapped_extend::{extend as gapped_extend, GappedParams};
use lastz_gxy::hsp::{extend_hit, extend_hit_scalar, Hsp, HspParams};
use lastz_gxy::output::{maf::MafWriter, paf::PafWriter, Record, Strand};
use lastz_gxy::pos_table::PosTable;
use lastz_gxy::scoring::ScoringMatrix;
use lastz_gxy::seed_search::{search, SearchParams};
use lastz_gxy::seeds::SeedPattern;
use lastz_gxy::sequences::PackedSeq;

/// Deterministic pseudorandom nucleotide buffer.
fn synthetic_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        out.push(b"ACGT"[((state >> 33) & 3) as usize]);
    }
    out
}

fn bench_pos_table_build(c: &mut Criterion) {
    let seq = PackedSeq::from_ascii(&synthetic_dna(1_000_000, 1));
    let pat = SeedPattern::twelve_of_nineteen();
    c.bench_function("pos_table/1Mbp/12of19", |b| {
        b.iter(|| {
            let _ = PosTable::build(black_box(&seq), black_box(&pat), 1);
        })
    });
}

fn bench_seed_search(c: &mut Criterion) {
    let target = PackedSeq::from_ascii(&synthetic_dna(1_000_000, 1));
    let query = PackedSeq::from_ascii(&synthetic_dna(100_000, 2));
    let pat = SeedPattern::solid(12);
    let table = PosTable::build(&target, &pat, 1);
    let matrix = ScoringMatrix::hoxd70();
    let params = SearchParams {
        step: 1,
        hsp: HspParams::default(),
        transitions: 1,
        entropy_threshold: None,
    };
    c.bench_function("seed_search/1Mbp_vs_100kbp/match12", |b| {
        b.iter(|| {
            let hsps = search(
                black_box(&table),
                black_box(&target),
                black_box(&query),
                black_box(&matrix),
                black_box(&params),
            );
            black_box(hsps.len());
        })
    });
}

fn bench_hsp_extend(c: &mut Criterion) {
    // Sequences with a planted homologous run so the extender actually
    // runs a non-trivial number of columns before x-drop.
    let payload = synthetic_dna(2_000, 42);
    let mut target = synthetic_dna(4_000, 1);
    let mut query = synthetic_dna(4_000, 2);
    target.splice(1_000..1_000, payload.iter().copied());
    query.splice(1_000..1_000, payload.iter().copied());
    let target = PackedSeq::from_ascii(&target);
    let query = PackedSeq::from_ascii(&query);
    let matrix = ScoringMatrix::hoxd70();
    let params = HspParams { x_drop: 910, hsp_threshold: 0 };
    let mut grp = c.benchmark_group("hsp/extend_hit");
    grp.bench_function("dispatch", |b| {
        b.iter(|| {
            let h = extend_hit(
                black_box(&target),
                black_box(&query),
                black_box(&matrix),
                black_box(&params),
                1_000,
                1_000,
                12,
            );
            black_box(h);
        })
    });
    grp.bench_function("scalar", |b| {
        b.iter(|| {
            let h = extend_hit_scalar(
                black_box(&target),
                black_box(&query),
                black_box(&matrix),
                black_box(&params),
                1_000,
                1_000,
                12,
            );
            black_box(h);
        })
    });
    grp.finish();
}

fn bench_gapped_extend(c: &mut Criterion) {
    let target = PackedSeq::from_ascii(&synthetic_dna(2_000, 1));
    let query = {
        // Introduce a deletion so the DP actually exercises traceback.
        let mut bytes = synthetic_dna(2_000, 1);
        bytes.remove(1_000);
        bytes.remove(1_500);
        PackedSeq::from_ascii(&bytes)
    };
    let matrix = ScoringMatrix::hoxd70();
    let hsp = Hsp { t_start: 500, q_start: 500, length: 12, score: 0 };
    let anchor = choose_anchor(hsp, &target, &query, &matrix, 31);
    let params = GappedParams { y_drop: 9_400, gapped_threshold: 0 };
    c.bench_function("gapped_extend/2kbp_with_indel", |b| {
        b.iter(|| {
            let aln = gapped_extend(
                black_box(anchor),
                black_box(&target),
                black_box(&query),
                black_box(&matrix),
                black_box(&params),
            );
            black_box(aln);
        })
    });
}

fn bench_chain(c: &mut Criterion) {
    // 200 HSPs scattered across a target; chain should pick a subset.
    let mut hsps: Vec<Hsp> = Vec::new();
    for i in 0..200u32 {
        hsps.push(Hsp {
            t_start: i * 100 + (i.wrapping_mul(37)) % 25,
            q_start: i * 100 + (i.wrapping_mul(41)) % 25,
            length: 50,
            score: 2_000 + (i as i32 % 50),
        });
    }
    c.bench_function("chain/best_of_200", |b| {
        b.iter(|| {
            let out = best_chain(black_box(&hsps));
            black_box(out.len());
        })
    });
}

fn bench_output(c: &mut Criterion) {
    let mut script = EditScript::new();
    script.push(EditOp::Match, 200);
    script.push(EditOp::InsertQuery, 2);
    script.push(EditOp::Match, 100);
    script.push(EditOp::DeleteQuery, 3);
    script.push(EditOp::Match, 150);
    let rec = Record {
        target_name: "chrT".into(),
        target_len: 1_000_000,
        query_name: "chrQ".into(),
        query_len: 1_000,
        query_strand: Strand::Plus,
        t_start: 500_000,
        q_start: 100,
        t_span: script.target_span(),
        q_span: script.query_span(),
        score: 4_200,
        script,
        target_bases: synthetic_dna(455, 3),
        query_bases: synthetic_dna(452, 4),
    };

    c.bench_function("output/maf/write_record", |b| {
        b.iter(|| {
            let mut out = Vec::with_capacity(2048);
            let mut w = MafWriter::new(&mut out);
            w.write_record(black_box(&rec)).unwrap();
            black_box(out.len());
        })
    });
    c.bench_function("output/paf/write_record", |b| {
        b.iter(|| {
            let mut out = Vec::with_capacity(1024);
            let mut w = PafWriter::new(&mut out);
            w.write_record(black_box(&rec)).unwrap();
            black_box(out.len());
        })
    });
}

criterion_group!(
    benches,
    bench_pos_table_build,
    bench_seed_search,
    bench_hsp_extend,
    bench_gapped_extend,
    bench_chain,
    bench_output
);
criterion_main!(benches);
