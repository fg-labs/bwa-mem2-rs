# bwa-mem3-rs

Safe Rust API for [bwa-mem3](https://github.com/fg-labs/bwa-mem3) alignment with packed-BAM output and caller-owned parallelism.

## What it does

- Aligns read pairs via bwa-mem3, returning **packed BAM records** per the BAM spec — ready to stream to any BGZF writer.
- Blocking, reentrant, single-threaded per call. **The caller drives all parallelism** via rayon, tokio `spawn_blocking`, or any work-stealing pool.
- Phase-split API: `seed_batch` → `extend_batch` lets pools balance seeding and extension independently. `align_batch` is a thin wrapper over the pair.
- Runs on x86_64 (AVX2 default) and aarch64 (Apple Silicon + AWS Graviton).

## Quick start

```rust
use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

let idx = BwaIndex::load("hg38.fa")?;                 // bwa-mem3 index prefix
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

## Parallelism

```rust
use bwa_mem3_rs::*;
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

The crate never spawns threads. Every public function blocks its calling thread and is safe to call concurrently on the same `BwaIndex`.

## Index

Build it with the upstream CLI:

```
bwa-mem3 index <prefix>.fa
```

Index building is intentionally out of scope here.

## See also

- [Workspace README](https://github.com/fg-labs/bwa-mem3-rs#readme) — architecture, tests, contribution.
- [`bwa-mem3-rs-cli`](https://crates.io/crates/bwa-mem3-rs-cli) — thin `bwa-rs` CLI binary over this crate.
- [`bwa-mem3-sys`](https://crates.io/crates/bwa-mem3-sys) — the unsafe FFI crate this one wraps.

## License

MIT.
