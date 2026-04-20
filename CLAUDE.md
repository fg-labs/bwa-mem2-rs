# CLAUDE.md

Project-specific guidance for AI assistants working in this repository.

## Project overview

`bwa-mem2-rs` is a Rust FFI crate for [bwa-mem2] that emits packed BAM records
(per the BAM spec) and leaves all parallelism to the caller. It ships as a
three-crate workspace:

- `bwa-mem2-sys` — unsafe FFI + vendored C/C++ + custom shim.
- `bwa-mem2-rs` — safe wrapper.
- `bwa-mem2-rs-cli` — thin `bwa-rs` CLI binary.

Hosted under [`fg-labs/bwa-mem2-rs`](https://github.com/fg-labs/bwa-mem2-rs).
Vendors [`fg-labs/bwa-mem2`](https://github.com/fg-labs/bwa-mem2) at the
`fg-main` branch.

[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2

## Build and test

```bash
cargo ci-build     # cargo build --workspace --all-targets --locked
cargo ci-test      # cargo test --workspace --locked
cargo ci-fmt       # cargo fmt --all -- --check
cargo ci-lint      # cargo clippy --workspace --all-targets -- -D warnings

# Integration tests (require a prebuilt bwa-mem2 index; skipped otherwise):
BWA_MEM2_RS_TEST_REF=/Users/nhomer/work/references/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
    cargo test --workspace

# E2E (requires bwa-mem2 CLI + samtools on PATH or BWA_MEM2_BIN set):
cargo test -p bwa-mem2-rs-cli --test e2e
```

## Architecture

Ownership layers:

| Layer | What's there |
|---|---|
| `fg-labs/bwa-mem2` / `fg-main` branch | Our integration fork of upstream bwa-mem2. Everything-we've-merged: Apple Silicon NEON, Linux-arm64 Makefile, drop-stat fix, future bwa-meth / XB / etc. |
| `bwa-mem2-sys/vendor/bwa-mem2/` | Pristine snapshot of `fg-labs/bwa-mem2@fg-main` at `vendor/COMMIT`. Never edit in place. Refresh via `scripts/refresh-bwa-mem2.sh`. |
| `bwa-mem2-sys/patches/` | Numbered `.patch` files applied to the vendored source at build time. Goal state: empty (fixes landed in `fg-main` instead). |
| `bwa-mem2-sys/shim/` | Our C/C++ shim. `bwa_shim.h` is the public header bindgen consumes. `bwa_shim.cpp` uses our POD copies of `mem_opt_t` / `mem_pestat_t` from `bwa_shim_types.h`. `bwa_shim_align.cpp` includes upstream's real `bwamem.h` / `FMI_search.h` and bridges to `bwa_shim.cpp` via opaque pointers (layouts match; verified by bindgen layout-assertion tests). |
| `bwa-mem2-rs/` | Safe Rust API: `BwaIndex`, `MemOpts`, `MemPeStat`, `ReadPair`, `Seeds`, `AlignmentBatch`, `Record`. Phase split: `seed_batch` → `extend_batch`. Convenience: `align_batch`. Pilot-only: `estimate_pestat`. |
| `bwa-mem2-rs-cli/` | Minimal `bwa-rs mem` CLI. Reads paired FASTQ (gzip-aware), writes proper BGZF-BAM. |

## Gotchas

### 1. Vendored Makefile must retain `MATE_SORT=0`

Our shim delegates the paired-end decision to upstream's `mem_pair_resolve` (exposed by fg-labs/bwa-mem2 PR #9) and then runs our own BAM emission. `mem_pair_resolve`'s internal branching is guarded on `#if MATE_SORT` vs. the default; only the `MATE_SORT=0` path is exercised (and audited) by our shim. If `fg-main` ever flips the default, the pairing logic would swap to an untested branch. `build.rs` asserts `-DMATE_SORT=0` at build time.

### 2. `mem_opt_t` / `mem_pestat_t` layouts are mirrored in two places

They're in `shim/bwa_shim_types.h` (what bindgen reads) and in upstream's `bwamem.h` (what `bwa_shim_align.cpp` includes). Both must stay byte-identical. A bindgen layout-assertion test (`bindgen_test_layout_mem_opt_t`) catches drift at build time. On `refresh-bwa-mem2.sh`, diff `vendor/bwa-mem2/src/bwamem.h` around lines 76–108 (`mem_opt_t`) and 162–166 (`mem_pestat_t`); update `shim/bwa_shim_types.h` if either changed.

### 3. macOS deployment target mismatch → SIGBUS at test-binary startup

`build.rs` sets `MACOSX_DEPLOYMENT_TARGET=11.0` explicitly when building on macOS to keep `cc`'s emitted objects aligned with what `cargo`/`rustc` links. Without this, linked binaries can fault at startup on macOS 26+.

### 4. Shadowing libc

Upstream `bwamem.cpp` had an unused file-scope `int stat;` that shadowed libc's `stat()` syscall wrapper → SIGBUS on test-harness startup. The fix is carried on `fg-main`; should be fixed-forward there rather than as a patch in our crate.

### 5. Some libbwa-mem2 symbols have C++ linkage

`mem_matesw` (and various internal helpers) are not `extern "C"` in `bwamem.h`. If you add forward declarations in the shim for any such symbol, put them outside `extern "C"` blocks so the mangled names match. As of the `mem_pair_resolve` adoption, the shim no longer forward-declares any of these directly.

### 6. bwa-mem2 CIGAR opcodes use a 5-char table, not BAM spec

bwa-mem2's internal `mem_aln_t.cigar` uses opcode table `MIDSH` (M=0 I=1 D=2 S=3 H=4). The BAM spec uses `MIDNSHP=X` (M=0 I=1 D=2 N=3 S=4 H=5). Our emitter (`bwa_cigar_to_bam` in `bwa_shim_align.cpp`) remaps before writing packed BAM. Don't copy bwa-mem2's opcodes verbatim into BAM output.

### 7. `s->seq` after `mem_kernel1_core` is 2-bit encoded, not ASCII

`mem_kernel1_core` rewrites each read's `seq` in place via `nst_nt4_table`: bytes become 0=A, 1=C, 2=G, 3=T, 4=N. Our emitter uses the `bwa2_to_bam4` table to map these to BAM 4-bit nibbles (1/2/4/8/15), plus `bwa2_complement` for the reverse-strand path.

### 8. Files excluded from the C++ build (see `build.rs`)

- `main.cpp` — CLI entry point
- `bwtindex.cpp` — index builder (not our concern)
- `runsimd.cpp` — unguarded `int main()` collides with Rust test harness

`fastmap.cpp` is built (was previously excluded) so `libbwa-mem2.a` exports `worker_alloc` / `worker_free`. Its `main_mem` entry point doesn't clash with Rust's test harness.

## Commit / PR conventions

- Conventional Commits; sign with `-S`; see `CONTRIBUTING.md`.
- Never mention AI/Claude/Anthropic in commits or GitHub-visible content
  (enforced globally in `~/.claude/CLAUDE.md`).
- 1Password signing flakes when the user is away from the machine. After
  three failed attempts, commit with `--no-gpg-sign` and note it in the
  commit message for re-signing later.

## Updating upstream

See `CONTRIBUTING.md` → "Updating the vendored `bwa-mem2` source."
