# bwa-mem3-sys

Low-level FFI bindings to [bwa-mem3](https://github.com/fg-labs/bwa-mem3), vendored from [`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3) at the `main` branch.

**Most users want [`bwa-mem3-rs`](https://crates.io/crates/bwa-mem3-rs) instead** — this crate exists so that crate (and other sibling crates) can share a single vendored copy of bwa-mem3 plus the custom C/C++ shim. Depending on `bwa-mem3-sys` directly means writing `unsafe` FFI yourself.

## What you get

- The full bwa-mem3 source (Apple Silicon / NEON-enabled), compiled as a static library at build time.
- `mem_opt_t` and `mem_pestat_t` POD structs, with bindgen-verified byte-identical layouts.
- A small C shim that exposes opaque `BwaIndex`, `BwaSeeds`, and `BwaBatch` handles plus per-batch alignment entry points that return packed BAM records directly from `mem_aln_t` (no SAM intermediate).
- Runs on x86_64 (AVX2 default, optional `sse42` / `avx512` / `native` cargo features) and aarch64 (Apple Silicon + AWS Graviton).

## Build requirements

- A C++17 compiler (clang on macOS, gcc on Linux).
- `zlib` development headers.
- On macOS, the build pins `MACOSX_DEPLOYMENT_TARGET=11.0` to keep the cc-produced objects aligned with rustc's link target.

## See also

- [Workspace README](https://github.com/fg-labs/bwa-mem3-rs#readme) — architecture, caveats, contribution.
- [`bwa-mem3-rs`](https://crates.io/crates/bwa-mem3-rs) — the safe Rust wrapper.
- [`bwa-mem3-rs-cli`](https://crates.io/crates/bwa-mem3-rs-cli) — a thin CLI over the safe API.

## License

MIT.
