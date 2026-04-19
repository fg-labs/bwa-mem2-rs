# bwa-mem2-rs-cli

Minimal `bwa-rs` CLI wrapping [`bwa-mem2-rs`](https://crates.io/crates/bwa-mem2-rs): aligns paired FASTQ to BGZF-compressed BAM.

## Install

```bash
cargo install bwa-mem2-rs-cli
```

## Usage

```bash
bwa-rs mem <ref-prefix> <r1.fq[.gz]> <r2.fq[.gz]> -o out.bam
```

- `<ref-prefix>` — a prebuilt bwa-mem2 index (build one with `bwa-mem2 index <ref.fa>`).
- `<r1.fq[.gz]>` / `<r2.fq[.gz]>` — paired FASTQ; `.gz` is decompressed transparently.
- `-o` / `--output` — output BAM path; stdout if omitted.
- `-k` / `--min-seed-len` — seed length override (same as upstream `-k`).
- `--batch-size` — pairs per alignment batch (default 1024).

Output is a proper BAM stream (magic + header + concatenated packed records + BGZF EOF). Pipe it to `samtools view` / `samtools sort` directly.

## See also

- [Workspace README](https://github.com/fg-labs/bwa-mem2-rs#readme) — architecture, caveats.
- [`bwa-mem2-rs`](https://crates.io/crates/bwa-mem2-rs) — the library behind this binary.

## License

MIT.
