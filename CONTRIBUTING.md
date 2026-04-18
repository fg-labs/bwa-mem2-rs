# Contributing to bwa-mem2-rs

## Development setup

### Prerequisites

- Rust (stable toolchain — see `rust-toolchain.toml` for the pinned version).
- A C++17 compiler (clang on macOS, gcc on Linux).
- `zlib` development headers (`zlib1g-dev` on Debian/Ubuntu; included on macOS).
- For integration tests: a prebuilt bwa-mem2 index — set `BWA_MEM2_RS_TEST_REF`
  to its prefix (the path such that `<prefix>.bwt.2bit.64` exists). Build one
  with `bwa-mem2 index <ref.fa>`.

### Install git hooks

Pre-commit hooks run `cargo ci-fmt` and `cargo ci-lint` before each commit:

```bash
./scripts/install-hooks.sh
```

The hook is skippable with `git commit --no-verify` when you know what you're
doing.

### Run checks manually

```bash
cargo ci-fmt      # format check
cargo ci-lint     # clippy with -D warnings
cargo ci-build    # build everything
cargo ci-test     # unit + integration tests

# or everything at once:
cargo ci-fmt && cargo ci-lint && cargo ci-build && cargo ci-test
```

### Run integration tests against a real index

```bash
BWA_MEM2_RS_TEST_REF=/path/to/hg38.fa cargo test --workspace
```

Without `BWA_MEM2_RS_TEST_REF` the integration tests skip gracefully.

### End-to-end regression

`bwa-mem2-rs-cli` ships a `tests/e2e.rs` test that embeds PhiX174, builds an
index via `bwa-mem2 index`, simulates 1,000 paired reads, runs our `bwa-rs`
CLI, and verifies the output via `samtools quickcheck` + `samtools view -c`.
Runs in CI when `bwa-mem2` + `samtools` are available; locally set
`BWA_MEM2_BIN=/path/to/bwa-mem2` if the CLI is not on `PATH`.

### CLI-parity test (optional)

The `cli_parity` test compares our packed-BAM output against `bwa-mem2 mem`
CLI output. It requires:

- `BWA_MEM2_BIN=/path/to/bwa-mem2`
- `BWA_MEM2_RS_TEST_REF=/path/to/index/prefix`
- `samtools` on `PATH`

## Code style

- Follow the Rust API Guidelines; prefer `impl Trait` over `Box<dyn Trait>` on
  returns; `&[T]` over `&Vec<T>` on borrows.
- `#![forbid(unsafe_code)]` would be ideal but is infeasible for an FFI crate;
  every `unsafe impl Send/Sync` and unsafe block must carry a safety comment
  naming the invariant it relies on (e.g. "BwaIndex has no mutable state after
  load" or "Seeds owns all memory; no aliasing back into shared state").
- Errors return `Result<_, bwa_mem2_rs::Error>`; no panics across the FFI
  boundary.
- Clippy runs with `-D warnings`.

## Commit conventions

- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, etc.
  per <https://conventionalcommits.org>).
- Sign commits (`-S`). Session commits made without signing should be re-signed
  via `git rebase --exec 'git commit --amend --no-edit -S'` before merge.
- Branch names: `JIRA-1234/initials_short-description` or similar.
- Never commit directly to `main`.

## Updating the vendored `bwa-mem2` source

The crate vendors a snapshot of [`fg-labs/bwa-mem2`](https://github.com/fg-labs/bwa-mem2)
at the `fg-main` branch under `bwa-mem2-sys/vendor/bwa-mem2/`. `fg-main` is
our integration branch tracking upstream `bwa-mem2/bwa-mem2` with:

- Apple Silicon / NEON hot-path kernels (PR #288 equivalent).
- `sse2neon` bridge + Linux/macOS aarch64 Makefile targets (PR #271
  equivalent).
- `drop-unused-global-stat` fix (the `int stat;` in `bwamem.cpp` that
  collides with libc on macOS).
- Future: bwa-meth support, XB tag, etc.

To refresh the vendored snapshot to a new `fg-main` tip:

```bash
scripts/refresh-bwa-mem2.sh <commit-hash> [local-source-path]
```

The script:

1. Clones (or rsyncs from a local tree) the target commit into `vendor/`.
2. Writes the commit hash into `vendor/COMMIT`.
3. Verifies the Makefile still has `-DMATE_SORT=0` (the shim's pairing logic
   depends on that default).

After refreshing, run the full test suite against a real index to catch any
upstream-API drift (`mem_kernel1_core`, `mem_kernel2_core`, `mem_matesw`,
etc.).

## Releases

Releases are managed by [release-plz](https://release-plz.dev). Merging to
`main` triggers a release PR that bumps versions and generates `CHANGELOG.md`
entries from conventional commits. Merging the release PR triggers the
publish workflow.
