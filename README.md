# bwa-mem2-rs

Rust FFI crate for [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) with packed-BAM output and caller-owned parallelism.

Hosted under [`fg-labs/bwa-mem2-rs`](https://github.com/fg-labs/bwa-mem2-rs); it vendors [`fg-labs/bwa-mem2`](https://github.com/fg-labs/bwa-mem2) at the `fg-main` branch — our integration branch over upstream with Apple Silicon / NEON + Linux-arm64 / CI fixes applied.

## What it does

- Aligns read pairs via bwa-mem2, returning **packed BAM records** per the BAM spec (ready to stream to any BGZF writer).
- Blocking, reentrant, single-threaded per call — **caller drives all parallelism** via rayon, tokio `spawn_blocking`, or any work-stealing pool.
- Phase-split API: `seed_batch` → `extend_batch` lets work-stealing pools balance seeding and extension independently. `align_batch` is a thin convenience wrapper.
- Runs on **x86_64** (AVX2 by default) and **aarch64** (Apple Silicon + AWS Graviton).

## Quick start

```rust
use bwa_mem2_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

let idx = BwaIndex::load("hg38.fa")?;                 // bwa-mem2 index prefix
let mut opts = MemOpts::new()?;
opts.set_pe(true);

let pairs = vec![ReadPair {
    name_r1: b"read1",
    seq_r1:  b"ACGT...",
    qual_r1: Some(b"IIII..."),
    name_r2: b"read1",
    seq_r2:  b"TTTT...",
    qual_r2: Some(b"IIII..."),
}];

let (aln, _pestat) = align_batch(&idx, &opts, &pairs, None)?;
for record in aln.iter() {
    // record.bytes is a packed BAM record: [u32 le block_size][data].
    // Feed to noodles_bam, rust-htslib, or your own BGZF writer.
}
```

## Work-stealing parallelism

```rust
use bwa_mem2_rs::*;
use rayon::prelude::*;
use std::sync::Arc;

let idx = Arc::new(BwaIndex::load("hg38.fa")?);
let mut opts = MemOpts::new()?;
opts.set_pe(true);

let results: Vec<_> = batches
    .par_iter()
    .map(|batch| align_batch(&idx, &opts, batch, None))
    .collect::<Result<_, _>>()?;
```

The crate itself never spawns threads. Every public function blocks its calling thread and is safe to call concurrently from multiple threads sharing the same `BwaIndex`.

## CLI

The `bwa-mem2-rs-cli` crate ships a minimal `bwa-rs` binary that wraps the library:

```
bwa-rs mem <ref-prefix> <r1.fq[.gz]> <r2.fq[.gz]> [-o out.bam]
```

## Index

Use the upstream CLI to build an index:

```
bwa-mem2 index <prefix>.fa
```

Index building is intentionally out-of-scope for this crate.

## Architecture

- `bwa-mem2-sys` vendors upstream bwa-mem2 sources and compiles them alongside
  our custom C++ shim (`shim/bwa_shim*.cpp`). The shim bridges to upstream's
  public API via opaque pointers and emits packed BAM directly from bwa-mem2's
  `mem_aln_t` (no SAM round-trip).
- `bwa-mem2-rs` wraps the FFI crate in a safe Rust API with `Send`/`Sync`
  semantics that let callers drive parallelism externally.
- `bwa-mem2-rs-cli` is a thin `bwa-rs mem` binary.

See `CLAUDE.md` for the project-specific gotchas (CIGAR opcode remap, 2-bit
seq encoding, vendored `MATE_SORT=0` invariant, etc.).

## Testing

```bash
# Unit tests (no index required)
cargo test --workspace

# Integration tests (require a bwa-mem2 index prefix)
BWA_MEM2_RS_TEST_REF=/path/to/hg38.fa cargo test --workspace
```

The end-to-end regression test (PhiX reference + simulated reads + CLI round-trip) runs automatically in CI once `bwa-mem2` and `samtools` are available; locally set `BWA_MEM2_BIN=/path/to/bwa-mem2` if the CLI is not on `PATH`.

## License

MIT.
