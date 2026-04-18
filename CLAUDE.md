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

Our shim composes `mem_pair` + `mem_matesw` + `mem_mark_primary_se` + `mem_reg2aln` as a replacement for `mem_sam_pe`. This is only safe because MATE_SORT=0 (upstream default) leaves the MATE_SORT-specific helpers as dead code. If `fg-main` ever flips the default, the shim pairing would diverge. `build.rs` asserts this at build time.

### 2. bwa-mem2 header declarations are stale in places

`bwamem.h` declares `mem_kernel1_core` with 7 parameters; the actual definition in `bwamem.cpp` has 9 (includes `seedBuf` + `seedBufSize`). `mem_kernel2_core` is not declared at all. The shim (`bwa_shim_align.cpp`) carries correct forward declarations; don't rely on `bwamem.h` for these signatures.

### 3. `mem_opt_t` / `mem_pestat_t` layouts are mirrored in two places

They're in `shim/bwa_shim_types.h` (what bindgen reads) and in upstream's `bwamem.h` (what `bwa_shim_align.cpp` includes). Both must stay byte-identical. A bindgen layout-assertion test (`bindgen_test_layout_mem_opt_t`) catches drift at build time. On `refresh-bwa-mem2.sh`, diff `vendor/bwa-mem2/src/bwamem.h` around lines 76–108 (`mem_opt_t`) and 162–166 (`mem_pestat_t`); update `shim/bwa_shim_types.h` if either changed.

### 4. macOS deployment target mismatch → SIGBUS at test-binary startup

`build.rs` sets `MACOSX_DEPLOYMENT_TARGET=11.0` explicitly when building on macOS to keep `cc`'s emitted objects aligned with what `cargo`/`rustc` links. Without this, linked binaries can fault at startup on macOS 26+.

### 5. Shadowing libc

Upstream `bwamem.cpp` had an unused file-scope `int stat;` that shadowed libc's `stat()` syscall wrapper → SIGBUS on test-harness startup. The fix is carried on `fg-main`; should be fixed-forward there rather than as a patch in our crate.

### 6. mem_matesw / mem_kernel1_core have C++ linkage

Not `extern "C"` in `bwamem.h`. Forward declarations in our shim must be outside `extern "C"` blocks so the mangled names match.

### 7. bwa-mem2 CIGAR opcodes use a 5-char table, not BAM spec

bwa-mem2's internal `mem_aln_t.cigar` uses opcode table `MIDSH` (M=0 I=1 D=2 S=3 H=4). The BAM spec uses `MIDNSHP=X` (M=0 I=1 D=2 N=3 S=4 H=5). Our emitter (`bwa_cigar_to_bam` in `bwa_shim_align.cpp`) remaps before writing packed BAM. Don't copy bwa-mem2's opcodes verbatim into BAM output.

### 8. `s->seq` after `mem_kernel1_core` is 2-bit encoded, not ASCII

`mem_kernel1_core` rewrites each read's `seq` in place via `nst_nt4_table`: bytes become 0=A, 1=C, 2=G, 3=T, 4=N. Our emitter uses the `bwa2_to_bam4` table to map these to BAM 4-bit nibbles (1/2/4/8/15), plus `bwa2_complement` for the reverse-strand path.

### 9. Files excluded from the C++ build (see `build.rs`)

- `main.cpp` — CLI entry point
- `fastmap.cpp` — CLI-side batch driver
- `bwtindex.cpp` — index builder (not our concern)
- `runsimd.cpp` — unguarded `int main()` collides with Rust test harness

## Commit / PR conventions

- Conventional Commits; sign with `-S`; see `CONTRIBUTING.md`.
- Never mention AI/Claude/Anthropic in commits or GitHub-visible content
  (enforced globally in `~/.claude/CLAUDE.md`).
- 1Password signing flakes when the user is away from the machine. After
  three failed attempts, commit with `--no-gpg-sign` and note it in the
  commit message for re-signing later.

## Updating upstream

See `CONTRIBUTING.md` → "Updating the vendored `bwa-mem2` source."
