//! End-to-end smoke test: runs the library driver on a tiny fixture and
//! checks that the expected single HSP is emitted. Guards against
//! regressions in the Phase 1 pipeline (sequences → pos_table →
//! seed_search → driver → output).

use lastz_gxy::driver::{run, Config, StrandSpec};
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
        ..Config::default()
    };

    let recs = run(&targets, &queries, &cfg);
    assert_eq!(recs.len(), 1, "got records: {recs:#?}");
    let r = &recs[0];
    assert_eq!(r.target_name, "chrT");
    assert_eq!(r.query_name, "chrQ");
    assert_eq!(r.query_strand, Strand::Plus);
    assert_eq!(r.hsp.t_start, 11);
    assert_eq!(r.hsp.q_start, 0);
    assert_eq!(r.hsp.length, 20);
    assert_eq!(r.target_bases, b"GATTACACATGGCATGTCGA");
    assert_eq!(r.query_bases, b"GATTACACATGGCATGTCGA");
}
