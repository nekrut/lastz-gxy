# lastz-gxy

A Rust reimplementation of [LASTZ](https://github.com/lastz/lastz) focused on
**throughput on commodity multi-core CPUs — with an optional GPU backend for
the embarrassingly-parallel seeding and HSP stages** — while preserving
upstream sensitivity. Algorithmic choices (seed patterns, HSP extension,
chaining, 3-state affine gapped DP, interpolation) are carried over unchanged;
the speed comes from parallelism, cache-friendly data structures, SIMD on the
CPU hot-path, and a portable GPU compute path for seed+HSP — not from relaxing
the alignment model.

> Status: **Phase 3 in progress** — AVX2 SIMD HSP x-drop landed with a
> proptest-enforced bit-exact parity gate against the scalar reference.
> Full pipeline: seed → HSP (SIMD) → chain → anchor → gapped 3-state affine
> DP → tweener → MAF/PAF/SAM output. See [PLAN.md](PLAN.md) for the phased
> roadmap and parity gate.

## Quickstart

```bash
cargo build --release

# Gapped alignment (default), MAF output:
./target/release/lastz-gxy target.fa query.fa \
    --seed 12of19 --hspthresh 3000 --gappedthresh 3000 --format maf

# SAM output (headers + CIGAR + NM/AS tags):
./target/release/lastz-gxy target.fa query.fa --format sam

# Full sensitivity pipeline: chain + tweener interpolation + PAF:
./target/release/lastz-gxy target.fa query.fa \
    --chain --inner --inner-seed match12 --format paf

# Fast HSP-only pass:
./target/release/lastz-gxy target.fa query.fa --seed match12 --nogapped
```

## Parity harness

Run the harness against upstream lastz 1.04.52:

```bash
parity/scripts/build-upstream.sh      # clones + builds pinned lastz tag
cargo build --release                 # builds lastz-gxy + gxy-compare
parity/scripts/compare.sh parity/corpus/pseudocat.fa parity/corpus/pseudopig.fa
```

`gxy-compare` prints block-level Jaccard, per-block score delta, aligned-bp
delta, per-side-only signatures, and whether the PLAN.md release gate passes.
Current state on the `pseudocat.fa × pseudopig.fa` fixture (upstream's own
test data):

| Metric                                      | Value     |
|---------------------------------------------|-----------|
| Upstream blocks                             | 14        |
| lastz-gxy blocks                            | 31        |
| Jaccard                                     | **0.22**  |
| **Score delta on shared blocks**            | **0** (median + max) |
| Aligned bp delta                            | +111 %    |
| Release gate                                | **FAIL**  |

The zero score-delta on shared blocks says the gapped DP matches upstream
bit-for-bit *when we agree on the block boundary*. The +111 % aligned-bp
and the low Jaccard are about *which* blocks we emit — directly attributable
to three known gaps: no soft-masking (PLAN.md §6 `--masking`), no `--census`
low-complexity suppression, and subtle chain/tweener differences. Closing
those is what moves Jaccard from 0.22 toward the 0.99 release gate.

## Benchmarks

Reproducible micro-benchmarks via `cargo bench` (criterion). Current
single-thread baselines on a commodity x86_64 box:

| Stage                                 | Input                      | Time   |
|---------------------------------------|----------------------------|-------:|
| `pos_table` build                     | 1 Mbp, 12of19              | ~TBD   |
| `seed_search` (HSP ungapped)          | 1 Mbp × 100 kbp, match12   | ~9 ms  |
| `hsp::extend_hit`                     | per hit, 10 kbp box        | ~238 ns|
| `gapped_extend` (with 2 indels)       | 2 kbp                      | ~40 ms |
| `chain::best_chain`                   | 200 HSPs                   | ~15 µs |
| `maf::write_record`                   | 450-col block              | ~1.3 µs|
| `paf::write_record`                   | 450-col block              | ~1.0 µs|

These are the reference numbers Phase 3 SIMD and Phase 2.5 GPU work will
report speedups against.

Still to come (Phase 3+): SIMD ungapped x-drop, striped-vector gapped DP,
`wgpu` GPU backend, 2bit/HSX readers, within-target chunking.
`cargo test` exercises 90 unit tests and 1 end-to-end integration fixture.

## Why another lastz?

UCSC/Harris lastz remains the reference for pairwise mammalian-scale
alignment, but in 2026 the upstream code is:

1. **Single-threaded.** Pipelines rely on external splitters (`run_lastz.py`,
   Galaxy's `lastz_wrapper`, Snakemake shards) to saturate cores. Intra-process
   parallelism across query chunks, targets, or DP cells does not exist.
2. **Scalar.** The 3-state affine gapped extension and ungapped x-drop loops
   are portable C with no SIMD. Modern AVX2/AVX-512/NEON implementations
   (KSW2, WFA2, block-aligner) are 4–10× faster on the same recurrence.
3. **Monolithic.** One binary handles FASTA/2bit/HSX parsing, masking,
   seeding, HSP, chaining, DP, and eight output formats in ~50 kLOC of C with
   a hand-rolled build. Static binaries and reproducible builds are difficult.
4. **GPU variants are clanky and disjoint.** [SegAlign](https://github.com/gsneha26/SegAlign)
   and [KegAlign](https://github.com/gsneha26/KegAlign) are separate forks
   that diverge on output semantics, require a CUDA SDK, and lock you out of
   AMD/Apple hardware. The GPU story should be a feature flag on the main
   binary, not a fork.

`lastz-gxy` targets the gap between upstream lastz (1 core, portable) and
SegAlign (GPU-only, niche): a **single Rust binary that matches lastz output
within a parity gate, uses every CPU core by default, and optionally offloads
seeding + HSP to any GPU that speaks Vulkan/Metal/DX12 (via `wgpu`) or CUDA
(opt-in feature)**.

## Non-goals

- **New alignment algorithms.** Seed, HSP, chain, anchor, align, interp all
  match upstream. No WFA substitution in v1; banded SIMD DP is bolted onto the
  existing recurrence.
- **Gapped DP on the GPU.** Variable-length, branch-heavy traceback is where
  SegAlign's complexity comes from. Gapped extension stays on CPU in v1 — GPU
  covers seeding + ungapped HSP only (the 70 % of upstream wall time).
- **CUDA-only GPU path.** The default GPU backend is [`wgpu`](https://github.com/gfx-rs/wgpu)
  (Vulkan / Metal / DX12 / WebGPU) so the same binary runs on NVIDIA, AMD,
  Intel, and Apple Silicon. A `--features cuda` build via
  [`cudarc`](https://github.com/coreylowman/cudarc) is opt-in for peak
  NVIDIA throughput; the `wgpu` path remains authoritative for parity.
- **Quantum DNA (probabilistic seeds).** Rarely used; deferred behind a
  feature flag.
- **LAV/GFA/HSX writers in v1.** MAF, AXT, SAM, and PAF cover ~all modern
  pipelines; the others are straightforward to add later.
- **Drop-in replacement for every lastz flag.** CLI-compatible for common
  flags (`--format`, `--step`, `--seed`, `--notransition`, `--xdrop`,
  `--ydrop`, `--gappedthresh`, `--hspthresh`, `--chain`, `--masking`,
  `--ambiguous`, `--strand`, `--scores`). Long-tail flags tracked in PLAN.md.

## Parity guarantee (target)

Release gate for v1:

| Metric                                         | Threshold        |
|------------------------------------------------|------------------|
| Aligned-block Jaccard vs upstream (per chrom)  | ≥ 0.99           |
| Per-block score delta (median)                 | 0                |
| Per-block score delta (max, outside tie zone)  | ≤ 1              |
| Total aligned bp delta                         | ≤ 0.1 %          |
| Identity distribution KS-statistic             | ≤ 0.01           |
| GPU HSPs vs CPU HSPs (same input, same seed)   | **bit-exact set** |

The 1% Jaccard slack absorbs tie-breaking in seed ordering, diagonal-hash
eviction, and DP traceback — not sensitivity loss. The GPU backend is held to
a stricter bar: for the same `pos_table` and seed stream, the set of HSPs
(tuples of `(t_id, q_id, strand, t_start, t_end, q_start, q_end, score)`)
emitted by `--backend=gpu` must equal the set emitted by `--backend=cpu`.
Any divergence is a GPU-backend bug, not a parity-gate slack. See
[PLAN.md §Parity testing](PLAN.md#parity-testing) for the four-tier harness.

## Performance target

On a 32-core workstation (Zen 4, AVX-512, RTX 4090), for `human chr1 vs mouse
chr1` with `--step=20 --notransition --format=maf`:

| Implementation                              | Wall time | Speedup |
|---------------------------------------------|----------:|--------:|
| upstream lastz 1.04.52                      | baseline  |   1.0×  |
| lastz-gxy (CPU, 32 threads)                 | target    | 20–30×  |
| lastz-gxy (`--backend=gpu`, wgpu)           | target    | 50–80×  |
| lastz-gxy (`--backend=gpu --features cuda`) | target    | 80–120× |

CPU target reflects near-linear scaling across the seeding and HSP stages
(the current single-thread bottleneck per `perf record` on upstream) plus
2–3× on gapped DP from SIMD. GPU target assumes seeding + HSP fully offloaded
(consistent with SegAlign's ~100× on those stages) while chain → DP → tweener
remain on CPU, so Amdahl bounds the end-to-end speedup at ~10× the
CPU-fraction ratio. Hard numbers published once the MVP lands; see
`benches/` for reproducible micro-benchmarks.

## Architecture sketch

```
lastz-gxy/
├── src/
│   ├── cli.rs            - clap CLI, upstream flag compatibility layer
│   ├── sequences.rs      - FASTA / 2bit / HSX loader (noodles + mmap)
│   ├── masking.rs        - soft/hard masking, dynamic masking
│   ├── scoring.rs        - HOXD70 / user matrices, --inferscores
│   ├── seeds.rs          - spaced-seed patterns (12-of-19, 19-of-20, twin)
│   ├── pos_table.rs      - bucketed seed->positions index (SoA, radix)
│   ├── seed_search.rs    - query-side seed iteration, transition handling
│   ├── diag_hash.rs      - lock-free diagonal hash for hit dedup
│   ├── hsp.rs            - ungapped x-drop extension (SIMD inner loop)
│   ├── gpu/
│   │   ├── mod.rs        - backend trait, feature-gated dispatch
│   │   ├── wgpu/         - portable backend (Vulkan/Metal/DX12)
│   │   │   ├── seed.wgsl - seed lookup + diag packing kernel
│   │   │   └── hsp.wgsl  - ungapped x-drop kernel
│   │   └── cuda/         - optional CUDA backend (cudarc, --features cuda)
│   │       ├── seed.cu
│   │       └── hsp.cu
│   ├── chain.rs          - HSP chaining (lastz-style, not minimap2)
│   ├── anchor.rs         - sliding-window midpoint anchor selection
│   ├── gapped_extend.rs  - 3-state affine DP, banded + SIMD lanes
│   ├── tweener.rs        - inter-alignment interpolation stage
│   ├── output/
│   │   ├── maf.rs        - streaming MAF writer
│   │   ├── axt.rs
│   │   ├── sam.rs
│   │   └── paf.rs
│   ├── driver.rs         - pipeline: target-sharded rayon scheduler
│   └── bin/
│       └── gxy-compare.rs - MAF/AXT concordance tool for parity gate
├── benches/              - criterion microbenchmarks (hot path)
├── parity/
│   ├── fixtures/         - small FASTAs + golden MAFs from upstream
│   └── compare/          - release-gate harness
├── scripts/
│   ├── build-upstream.sh - pins lastz v1.04.52 as ground truth
│   └── compare.sh        - runs both, diffs with gxy-compare
├── Cargo.toml
├── PLAN.md
└── README.md
```

## Performance levers

| Lever                                               | Expected gain | Phase |
|-----------------------------------------------------|--------------:|------:|
| Rayon work-stealing across (target × strand) pairs  |        5–10×  |    1  |
| Within-target chunking with halo-joined HSPs        |        2–4×   |    2  |
| SoA `pos_table` with radix-bucketed seed hash       |        1.5–2× |    1  |
| SIMD ungapped x-drop (AVX2/NEON)                    |        2–3×   |    2  |
| Banded SIMD 3-state affine DP (KSW2-style lanes)    |        3–5×   |    3  |
| Lock-free diagonal hash (per-shard, epoch-reclaim)  |        1.3×   |    2  |
| `mmap` 2bit/HSX reference; zero-copy query slice    |        1.2×   |    1  |
| Streaming MAF/AXT writer, bounded channel           |     I/O-bound |    1  |
| GPU seed lookup + diag packing (`wgpu` compute)     |     30–60× ⁶  |  2.5  |
| GPU ungapped x-drop (warp-parallel over hits)       |     40–80× ⁶  |  2.5  |
| CUDA backend (opt-in, same kernels via `cudarc`)    |       1.3–1.6×⁷|    3  |
| Minimizer prefilter for `--step ≥ 10` seeds         |        1.5×   |    4  |

⁶ Over single-thread upstream, on the seed+HSP stages only. End-to-end gain
  is Amdahl-bounded by the CPU-resident chain/DP/tweener tail.
⁷ Over the `wgpu` backend on the same NVIDIA GPU.

Multiplicative gains do not compose linearly; the 20–30× CPU target and
50–80× GPU target are conservative estimates based on profiling upstream on a
16 Mbp vs 16 Mbp pair where seeding + HSP is ~70 % of wall time and gapped DP
~25 %.

## Dependencies

- [`noodles`](https://github.com/zaeleus/noodles) — FASTA, SAM, BAM
- [`rayon`](https://github.com/rayon-rs/rayon) — work-stealing thread pool
- [`clap`](https://github.com/clap-rs/clap) — CLI
- [`memmap2`](https://github.com/RazrFalcon/memmap2) — reference mmap
- [`wgpu`](https://github.com/gfx-rs/wgpu) — portable GPU compute (default GPU backend)
- [`cudarc`](https://github.com/coreylowman/cudarc) — optional CUDA backend (`--features cuda`)
- [`criterion`](https://github.com/bheisler/criterion.rs) — benchmarks
- [`insta`](https://github.com/mitsuhiko/insta) — snapshot tests

No C dependencies on the default build. `cargo build --release` produces a
single binary; the `wgpu` GPU path uses whatever Vulkan/Metal/DX12 driver the
host already has. `cargo build --release --features cuda` additionally links
against a system CUDA runtime.

## License

MIT. Same as the `lofreq-gxy` sibling project.

## Related work

- [lastz/lastz](https://github.com/lastz/lastz) — upstream, C, single-thread
- [SegAlign](https://github.com/gsneha26/SegAlign) — GPU port of lastz seeding/HSP
- [KegAlign](https://github.com/gsneha26/KegAlign) — SegAlign successor, MAF output
- [minimap2](https://github.com/lh3/minimap2) — different niche (long reads, not mammalian WGA)
- [nekrut/lofreq-gxy](https://github.com/nekrut/lofreq-gxy) — sibling project, same rewrite philosophy
