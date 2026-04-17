# lastz-gxy — Design Plan

This document records the design decisions and the phased roadmap for a Rust
reimplementation of [LASTZ](https://github.com/lastz/lastz). The guiding
principle: **keep every algorithmic choice upstream makes; change only the
implementation substrate.** Sensitivity must not regress.

## 1. Goals

1. Match upstream lastz 1.04.52 output within the parity gate (§5).
2. Deliver 20–30× wall-clock speedup on a 32-core CPU for mammalian-scale
   whole-genome-alignment workloads.
3. Ship a single static binary; no htslib/autotools/Python runtime.
4. Keep the CLI recognisable to existing lastz users for the common flag set.
5. Stay CPU-only and portable (x86_64 + aarch64). GPU is SegAlign's domain.

## 2. Non-goals

- New alignment algorithms (no WFA, no ksw2 substitution — just SIMD lanes
  over the existing 3-state affine recurrence).
- Quantum/probabilistic seeds (`--quantum`). Rare; deferred to Tier 3.
- LAV, GFA, HSX, and text-align output in v1. MAF/AXT/SAM/PAF only.
- `--inferscores` in v1 (useful but not hot-path; port in Tier 2).
- Self-alignment heuristics beyond what falls out of target-vs-target.

## 3. Architecture

### 3.1 Stage pipeline (preserved from upstream)

```
  target FASTA ──► pos_table ──┐
                               ├─► seed_search ──► diag_hash ──► hsp
  query FASTA ──► seed stream ─┘                                  │
                                                                  ▼
                             output ◄── tweener ◄── gapped_extend ◄── anchor ◄── chain
```

Every stage maps 1:1 to an upstream file
(`pos_table.c`, `seed_search.c`, `diag_hash.c`, `seed_search.c::process_hsp`,
`chain.c`, `segment.c`, `gapped_extend.c`, `tweener.c`). Module boundaries in
`src/` mirror these names so that parity debugging has a short path.

### 3.2 Concurrency model

- **Outer layer**: rayon work-stealing over `(target_seq, strand)` pairs. One
  shared immutable `pos_table`, one `Arc<ScoringMatrix>`, per-worker query
  cursors.
- **Inner layer** (Phase 2): for targets larger than `CHUNK_BP` (default
  10 Mbp), split into chunks with a `HALO_BP` overlap (default = max
  gapped-extension reach, ~50 kbp). HSPs that cross chunk boundaries are
  de-duplicated in a final merge pass keyed on `(target_id, strand, diag,
  start)`.
- **Gapped DP parallelism**: not parallelised across cells; each DP invocation
  stays on one core. Parallelism comes from running many DP calls
  simultaneously.

### 3.3 Data layout

- `pos_table`: struct-of-arrays. `seed_hash → Range<u32>` into a packed
  `positions: Vec<u32>` array, radix-bucketed by the high bits of the seed to
  keep hot buckets in L2. Upstream uses a linked-list-per-bucket layout;
  flattening gives ~1.7× on micro-benchmarks of random 19-of-20 seed lookups.
- Sequences: `mmap` 2bit/HSX; ASCII FASTA is decoded once into a packed
  2-bit+N buffer.
- Diagonal hash: per-worker `FxHashMap<i32, DiagEntry>`. Per-worker avoids the
  cross-core contention upstream has no reason to worry about (it's
  single-threaded).

### 3.4 SIMD boundaries

- **Phase 2**: ungapped x-drop extension. Inner loop is a streaming
  score-diff accumulator on 2-bit-packed sequence; maps cleanly to 16× u8
  lanes in AVX2 / 32× in AVX-512 / 16× in NEON. Scalar fallback retained.
- **Phase 3**: 3-state affine DP. Striped-vector lanes (Farrar-style)
  over the query, one DP call per lane. Banded (`band = ydrop / gap_extend`);
  x/y-drop termination checked once per lane-wide stripe.
- Phase 1 ships scalar-only and still targets 5–10× end-to-end from threading
  and data-layout alone.

## 4. Phased roadmap

Each phase has a release gate. No phase ships without passing its gate.

### Phase 1 — Scaffolding + seed/HSP MVP (scalar, threaded)

Deliverables:
- `cli.rs` accepts the common flag set; unknown flags error out loudly with a
  "not implemented in gxy" pointer.
- `sequences.rs` loads FASTA and 2bit via `noodles` + `memmap2`.
- `scoring.rs` parses lastz `--scores` files; HOXD70 as default built-in.
- `seeds.rs` implements 12-of-19 and 19-of-20 patterns with transitions.
- `pos_table.rs` builds the SoA index with radix buckets.
- `seed_search.rs` + `diag_hash.rs` + `hsp.rs` produce ungapped HSPs.
- `driver.rs` shards `(target, strand)` across a rayon pool.
- `output/maf.rs` + `output/axt.rs` stream results.

Gate: on `test_data/`-equivalent fixtures, Jaccard ≥ 0.99 on ungapped HSPs vs
upstream `--nogapped`; end-to-end wall time ≥ 5× faster on 16 cores.

### Phase 2 — Chain, anchor, gapped extension (scalar DP + SIMD HSP)

Deliverables:
- `chain.rs` — lastz-style chaining (not minimap2's), preserves upstream
  chain selection.
- `anchor.rs` — sliding-window midpoint.
- `gapped_extend.rs` — 3-state affine DP with x/y-drop, scalar.
- `output/sam.rs` + `output/paf.rs`.
- SIMD ungapped x-drop behind `#[cfg(target_feature="avx2")]` / NEON.
- Within-target chunking with halo-join de-duplication.

Gate: full parity gate (§5) at Jaccard ≥ 0.99 on the benchmark corpus;
≥ 15× on 32 cores for human-chr1 vs mouse-chr1.

### Phase 3 — SIMD gapped DP + tweener

Deliverables:
- Striped-vector 3-state affine DP (AVX2 / NEON lanes). Scalar retained as
  reference impl for differential testing.
- `tweener.rs` — interpolation stage.
- `--inferscores` port.
- Criterion benches for each stage with upstream as baseline.

Gate: ≥ 20× on 32 cores on the full corpus; parity gate maintained.

### Phase 4 — Long tail

- Minimizer prefilter for `--step ≥ 10`.
- LAV / GFA / HSX / text-align writers as compile-time features.
- Quantum DNA seeding.
- CRAM target support via `noodles-cram`.
- Self-alignment heuristics beyond target=query.

## 5. Parity testing

Four tiers, mirroring `lofreq-gxy`.

### Tier 1 — Unit / property (every PR)

- Per-stage unit tests with hand-crafted inputs (known seed hits, known HSPs,
  known DP scores).
- Property tests (`proptest`): for random scoring matrices and random
  sequences, scalar DP = SIMD DP bit-exact, HSP score is maximal on its
  diagonal, seed lookups are order-independent.

### Tier 2 — Golden-MAF differential (every PR)

- `parity/fixtures/` holds ~12 MB of small FASTAs (SARS-CoV-2 pairs,
  yeast chrI pairs, human chrM vs chimp chrM, a synthetic 1 Mbp pair with
  injected indels).
- `scripts/compare.sh` runs upstream lastz v1.04.52 and lastz-gxy on each
  fixture, diffs with `bin/gxy-compare` which computes:
  - aligned-block Jaccard (block identity = target_id + strand + t_start +
    t_end + q_start + q_end, rounded to 1 bp)
  - per-matched-block score delta histogram
  - total aligned bp delta
  - identity distribution KS-statistic
- Stored as snapshots via `insta`.

### Tier 3 — Simulated truth (nightly)

- Use `pbsim3` / custom evolver to generate pairs with known indel/SNV
  history. Measure recall and precision of injected homologies for both
  implementations; require lastz-gxy ≥ upstream on both.

### Tier 4 — Differential fuzzing (pre-release)

- Random FASTA pairs fed to both binaries for N CPU-hours; any parity-gate
  violation is a release blocker.

### Release gate (v1.0)

All four tiers green, plus the per-chromosome numeric thresholds listed in
`README.md` on the full benchmark corpus (`parity/corpus/`):

- human chr1 vs mouse chr1 (HOXD70, `--notransition --step=20`)
- human chr21 vs chimp chr21 (identity regime, `--step=1`)
- yeast-vs-yeast all-vs-all (short targets, many shards)
- SARS-CoV-2 multi-genome (pathological repeats)

## 6. CLI compatibility matrix (v1 MVP)

| Upstream flag                        | v1 status |
|--------------------------------------|-----------|
| `--format={maf,axt,sam,paf}`         | ✅ |
| `--format={lav,gfa,hsx,text}`        | ❌ Phase 4 |
| `--seed={12of19,19of20,match<k>}`    | ✅ |
| `--seed=half`                        | ❌ Phase 4 |
| `--step=N`                           | ✅ |
| `--notransition / --transition[=2]`  | ✅ |
| `--xdrop=N / --ydrop=N`              | ✅ |
| `--hspthresh=N / --gappedthresh=N`   | ✅ |
| `--chain / --nochain`                | ✅ |
| `--gapped / --nogapped`              | ✅ |
| `--masking=N / --census`             | ✅ (soft+dynamic); `--census` Phase 2 |
| `--ambiguous={n,iupac}`              | ✅ |
| `--scores=FILE`                      | ✅ |
| `--strand={plus,minus,both}`         | ✅ |
| `--inferscores`                      | ❌ Phase 3 |
| `--quantum`                          | ❌ Phase 4 |
| `--self`                             | ⚠️ works by `target==query`; Phase 4 heuristics |
| `--allocate:*`                       | ❌ ignored (Rust allocator) |

Unknown/unsupported flags fail fast with
`error: flag '--foo' not implemented in lastz-gxy; see PLAN.md §6`.

## 7. Dropped features (intentional)

- **`--allocate:*`** — Rust's allocator handles sizing.
- **Text-align debug output** — low value, high surface area.
- **`--inner=FILE` recursive configs** — all config via CLI + optional TOML.
- **`--runtime` stats dump** — integrated into `tracing` spans instead.

## 8. Risks & mitigations

| Risk                                                         | Mitigation |
|--------------------------------------------------------------|------------|
| Tie-breaking differences in chain selection drop Jaccard     | Port upstream ordering byte-for-byte; parity fixtures catch regressions |
| SIMD DP disagrees with scalar in edge cases                  | Scalar kept as reference; property tests on every PR |
| `noodles` FASTA perf on huge references                      | Benched in Phase 1; fallback is direct `memmap2` reader |
| Seed-hash collisions differ from upstream and change HSPs    | Use upstream's hash bitwidth; Tier 2 gate catches it |
| Memory blow-up from holding whole reference in RAM           | `mmap` 2bit keeps RSS flat; benchmarked on chr1-scale |
| 20× target unreachable on small targets                      | Document: gains are WGA-scale; short-target overhead is acceptable |

## 9. Out of scope, but on the radar

- WFA2 as an alternative gapped extender (would change parity gate; tracked
  as a future fork, not a v1 swap).
- GPU path (would re-converge with SegAlign; not a goal for this repo).
- Streaming target from stdin / pipes (useful for progressive alignment
  pipelines).
