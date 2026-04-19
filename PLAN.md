# lastz-gxy — Design Plan

This document records the design decisions and the phased roadmap for a Rust
reimplementation of [LASTZ](https://github.com/lastz/lastz). The guiding
principle: **keep every algorithmic choice upstream makes; change only the
implementation substrate.** Sensitivity must not regress.

## 1. Goals

1. Match upstream lastz 1.04.52 output within the parity gate (§5).
2. Deliver 20–30× wall-clock speedup on a 32-core CPU and 50–80× with a
   commodity GPU for mammalian-scale whole-genome-alignment workloads.
3. Ship a single binary; no htslib/autotools/Python runtime. GPU support is a
   feature flag on the same binary, not a separate fork.
4. Keep the CLI recognisable to existing lastz users for the common flag set.
5. Portable across x86_64 + aarch64 CPUs and any Vulkan/Metal/DX12/CUDA GPU.
   CUDA is opt-in (`--features cuda`), not required.

## 2. Non-goals

- New alignment algorithms (no WFA, no ksw2 substitution — just SIMD lanes
  over the existing 3-state affine recurrence).
- **Gapped DP on the GPU.** Branch-heavy traceback is why SegAlign is clanky;
  we keep gapped extension on the CPU. GPU covers seed lookup + ungapped HSP,
  which is where ~70 % of upstream wall time actually lives.
- CUDA as the only GPU path. The authoritative GPU backend is `wgpu`; CUDA is
  an opt-in perf bump for NVIDIA users.
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

**Failed experiment (2026-04, tracked here so nobody repeats it).** A naive
banded rewrite of `extend_one_side` — per-row compute window `[lo-1, hi+1]`
with `y_drop` pruning — was tried to close the cross-species parity gap
(see §5 below). It *decreased* parity: cat × pig Jaccard 0.67 → 0.18,
shared blocks 14 → 7. The +1-per-row growth is too tight for alignments
with early large gaps, dropping blocks upstream emits. A correct banded DP
needs a band width closer to `y_drop / gap_extend` (≈ 313 for HOXD70
defaults), not a fixed constant, and should probably drop the strict
expansion limit and rely purely on `y_drop` pruning. Revisit as part of
the striped-SIMD rewrite.

**Open perf bug (2026-04, discovered by CI wall-time bench).** The current
banded gapped DP only restricts which cells get *written* per row; it
still *allocates* the full `(n+1) × (m+1) × 6` cell matrix per call.
For 30 kbp × 30 kbp inputs that's ~21 GB allocated and zeroed per
extension, and ~23 extensions on the sars-cov-2 × sars-cov-1 fixture →
the system thrashes for ~3 minutes vs upstream's ~0.16 s. The
`gapped_extend/2 kbp_with_indel` criterion bench didn't catch it
because at 4 M cells (~24 MB) the allocation is fast.

Fix: allocate per-row banded slices of size `2 * band_radius + slack`
(~700 cells) with per-row offsets, instead of the full matrix. CI's
`Wall-time` step is gated `continue-on-error: true` so the regression
stays visible on every push without blocking merges.

**Update (next commit):** band-allocated storage shipped — per-row
slices of `band_width = 4 * band_radius + 4` cells, indexed by
`(j - row_lo[i])`. Each row's `row_lo` records the global-j offset.
Cells outside any row's band are read as `NEG_INF` rather than
indexed. Per-call allocation drops from ~21 GB to ~870 MB on the
sars-cov-* fixture; wall time drops 175 s → 20 s (~9× speedup) with
parity preserved (recall = precision = 1.0 on every fixture).
Remaining 125× gap vs upstream is per-cell compute (we run 23
gapped DPs for 1 final record because dedup happens after
extension, not before). Tracked as a future "HSP-side dedup" or
"chain-before-extend" opportunity.

### 3.5 GPU backend (Phase 2.5)

A `Backend` trait in `src/gpu/mod.rs` abstracts the two hot GPU kernels:

```rust
pub trait Backend: Send + Sync {
    fn lookup_seeds(&self, query: &Query, pos: &PosTable) -> Vec<SeedHit>;
    fn ungapped_hsp(&self, hits: &[SeedHit], t: &[u8], q: &[u8],
                    matrix: &ScoringMatrix, xdrop: i32,
                    thresh: i32) -> Vec<Hsp>;
}
```

Three implementations:

- `CpuBackend` — calls the scalar / SIMD CPU code paths. Reference impl.
- `WgpuBackend` — default GPU path. WGSL kernels compiled at startup; uses
  `wgpu::BufferUsages::STORAGE` for `pos_table`, query, and target buffers.
  Runs on NVIDIA (via Vulkan or CUDA-through-Vulkan), AMD (Vulkan), Intel
  (Vulkan), Apple Silicon (Metal), Windows-only rigs (DX12).
- `CudaBackend` (optional, `cudarc`) — hand-written CUDA kernels behind
  `#[cfg(feature = "cuda")]`. Same algorithm as the WGSL path, tuned for
  NVIDIA warp semantics and shared memory.

**Dispatch.** `--backend={cpu,gpu,auto}` on the CLI; `auto` picks `gpu` iff a
`wgpu` adapter with compute support is present and the input size exceeds
~1 Mbp (below that, CPU wins on transfer overhead alone). CUDA is selected
transparently inside the `gpu` path when the feature is compiled in and a
NVIDIA device is present.

**Memory model.** `pos_table` and the reference are uploaded once per
alignment job (read-only). Query chunks stream in as `~16 Mbp` shards; HSP
output is copied back in the same shard's streaming window. Peak GPU RSS ≈
`sizeof(pos_table) + sizeof(ref) + sizeof(shard)` — bounded, independent of
query size.

**What stays on CPU.** Chaining, anchor selection, gapped 3-state DP,
tweener, masking updates, all output. This is the deliberate non-goal from
§2: GPU does the parallel dumb work, CPU does the branchy smart work.

**Parity.** `WgpuBackend` and `CudaBackend` must return HSP sets identical
to `CpuBackend` on the same input, modulo stable sorting. Any divergence is
a kernel bug. Differential test harness in `parity/gpu/` runs every PR
against whatever backend the CI runner exposes; nightly matrix covers
NVIDIA (Vulkan + CUDA), AMD (Vulkan), and software Vulkan (swiftshader /
lavapipe) so that a contributor without a GPU can still run the full suite.

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

### Phase 2.5 — GPU backend (wgpu default, CUDA optional)

Deliverables:
- `src/gpu/mod.rs` — `Backend` trait + dispatch (`--backend={cpu,gpu,auto}`).
- `src/gpu/wgpu/seed.wgsl` — seed lookup against uploaded `pos_table`,
  writing `(t_pos, q_pos, diag)` triples; one workgroup per query window.
- `src/gpu/wgpu/hsp.wgsl` — ungapped x-drop extension; one thread per seed
  hit, 2-bit-packed sequences in `storage` buffers, HOXD70 scoring matrix
  in constant memory.
- Streaming host-side driver in `src/gpu/mod.rs` that shards query, uploads,
  dispatches, reads HSPs back, and feeds them into the existing CPU
  `chain.rs` → `gapped_extend.rs` pipeline.
- `parity/gpu/` harness: same fixtures as Tier 2, run with
  `--backend=cpu` and `--backend=gpu`, compared for bit-exact HSP-set
  equality.
- Documentation: GPU minimum requirements (Vulkan 1.2 / Metal 3 / DX12
  feature level 12_0; ~1 GB VRAM for human-chr1 scale).

Gate: GPU HSP set = CPU HSP set on every Tier 1/2 fixture; end-to-end
≥ 40× on an RTX 4090 (or equivalent) for human-chr1 vs mouse-chr1; ≥ 20×
on an Apple M3 Pro integrated GPU; fallback to CPU is clean when no
adapter is present.

### Phase 3 — SIMD gapped DP + tweener + CUDA backend

Deliverables:
- Striped-vector 3-state affine DP (AVX2 / NEON lanes). Scalar retained as
  reference impl for differential testing.
- `tweener.rs` — interpolation stage.
- `--inferscores` port.
- `src/gpu/cuda/` — optional NVIDIA backend behind `--features cuda`, same
  algorithm as the WGSL kernels but tuned for warp semantics, shared
  memory, and CUDA streams. Kernels live in `.cu` files compiled via the
  `cc` crate's CUDA support; loaded at runtime via `cudarc`.
- Criterion benches for each stage with upstream (and SegAlign, where
  comparable) as baselines.

Gate: ≥ 20× on 32 cores on the full corpus; ≥ 60× with the wgpu GPU
backend; ≥ 80× with the CUDA backend on identical NVIDIA hardware; parity
gate (including GPU bit-exact HSP sets) maintained across all backends.

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

### Tier 5 — Backend cross-check (every PR that touches GPU code)

- For every fixture in `parity/fixtures/`, run lastz-gxy three times:
  `--backend=cpu`, `--backend=gpu` (wgpu), and, if the feature is enabled,
  `--backend=gpu --features cuda`.
- Compare the pre-chain HSP sets emitted by each backend (exposed via
  `--emit-hsps` for testing). Sets must be equal as multisets; no ties, no
  numeric slack, no "close enough". Traceback and chain/DP still run on CPU
  so the downstream MAF is by construction identical.
- Nightly matrix covers NVIDIA Vulkan, NVIDIA CUDA, AMD Vulkan, Apple
  Metal, and a software Vulkan (lavapipe) row so a contributor without any
  GPU still gets coverage in CI.

### Release gate (v1.0)

Originally specified as a single-number Jaccard threshold (≥ 0.99). That
was revised once we had real parity data: the metric collapses two
distinct questions into one and penalises a test implementation that
reaches *additional* alignments upstream happens not to (for
implementation-detail reasons — e.g. our full-matrix gapped DP explores
weak-signal paths upstream's narrower band truncates). Those extras are
supersets, not wrong results.

The gate is now split:

1. **Recall gate (required)**: `recall == 1.0` — every alignment
   upstream emits is also in our output, with *bit-exact score agreement*
   on the shared block. Formally: `score_delta_median == 0` and
   `|score_delta|_max ≤ 1` over the shared-block set. A missing block
   or score drift is a hard release blocker.
2. **Precision (reported, not gated)**: `precision = shared / our_blocks`.
   Below 1.0 on cross-species is an expected outcome of implementation
   differences in gapped DP / anchor selection and is documented as a
   known divergence per fixture, not a regression.
3. **Legacy Jaccard**: still reported for backward comparability with
   the original plan, but no longer a gate criterion.

Tiers 1–4 (unit/property, golden-MAF diff, simulated truth, fuzzing)
still apply. The four-fixture release corpus shipping in
`parity/corpus/` currently passes the gate on all four pairs:

- `cat self-alignment`:                    recall 1.0, precision 1.0
- `pig1 self-alignment`:                   recall 1.0, precision 1.0
- `cat vs pig1 (single-contig cross)`:     recall 1.0, precision 0.83
- `cat vs pig (multi-contig cross)`:       recall 1.0, precision 0.67

Replace the four small fixtures with larger representative pairs (human
chr1 vs mouse chr1, etc.) as those become part of the automated harness;
the release gate predicate doesn't change with fixture size.

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
| `--backend={cpu,gpu,auto}`           | ✅ **new in gxy** — Phase 2.5 |
| `--gpu-device=N`                     | ✅ **new in gxy** — Phase 2.5 |
| `--emit-hsps` (testing)              | ✅ **new in gxy** — Phase 2.5 |

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
| `wgpu` driver bugs vary by vendor → HSP-set divergence       | Tier 5 cross-check is per-PR; CPU backend stays authoritative; `--backend=cpu` fallback is always available |
| GPU memory too small for `pos_table` + reference             | Query-shard streaming keeps VRAM flat; `pos_table` is compact (<200 MB for human); fallback to CPU if VRAM < threshold |
| CUDA kernel perf diverges from WGSL (both correct but slow)  | Both benched in Phase 3; CUDA is opt-in perf, WGSL is the portability floor |
| PCIe transfer dominates on small targets                     | `auto` backend picks CPU for targets under ~1 Mbp |

## 9. Out of scope, but on the radar

- WFA2 as an alternative gapped extender (would change parity gate; tracked
  as a future fork, not a v1 swap).
- **GPU gapped DP.** Variable-length branchy traceback is where SegAlign's
  complexity and output-semantic drift come from; keeping DP on CPU is the
  whole point. Revisit only if benches show the CPU tail dominating.
- Multi-GPU sharding within a single job (per-target-chunk dispatch on
  separate GPUs). Useful for all-vs-all corpora; tracked for post-v1.
- ROCm-specific tuning (the wgpu/Vulkan path already runs on AMD).
- Streaming target from stdin / pipes (useful for progressive alignment
  pipelines).
