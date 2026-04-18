//! End-to-end regression test: run the `bwa-rs` CLI against a PhiX reference
//! with simulated reads, verify the output parses as valid BAM with the
//! expected record count.
//!
//! Requires `bwa-mem2` (for `index`) and `samtools` on PATH. Set
//! `BWA_MEM2_BIN` to point at the CLI if it's not on PATH. The test skips
//! gracefully when either is missing; CI installs them via conda.

mod phix_seq;

use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

fn find_bwa_mem2() -> Option<String> {
    if let Ok(p) = std::env::var("BWA_MEM2_BIN") {
        if std::path::Path::new(&p).exists() {
            return Some(p);
        }
    }
    let out = Command::new("which").arg("bwa-mem2").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let p = String::from_utf8(out.stdout).ok()?.trim().to_string();
    if p.is_empty() {
        None
    } else {
        Some(p)
    }
}

fn have_samtools() -> bool {
    Command::new("samtools").arg("--version").output().is_ok()
}

fn cli_bin() -> std::path::PathBuf {
    // Cargo sets CARGO_BIN_EXE_<name> for bin targets.
    std::env::var("CARGO_BIN_EXE_bwa-rs").map_or_else(
        |_| std::path::PathBuf::from(env!("CARGO_BIN_EXE_bwa-rs")),
        std::path::PathBuf::from,
    )
}

/// Deterministic PRNG + paired read simulator.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn simulate_pairs(
    ref_seq: &[u8],
    n: usize,
    read_len: usize,
    insert: usize,
    seed: u64,
) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    let mut rng = Rng(seed);
    let max_start = ref_seq.len().saturating_sub(insert + 1);
    (0..n)
        .map(|i| {
            let start = (rng.next() as usize) % max_start;
            let r1 = ref_seq[start..start + read_len].to_vec();
            let span = &ref_seq[start + insert - read_len..start + insert];
            let r2 = revcomp(span);
            (format!("r{i}"), r1, r2)
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

fn write_fastq(path: &std::path::Path, reads: &[(String, Vec<u8>)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in reads {
        let qual = vec![b'I'; seq.len()];
        writeln!(f, "@{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f, "\n+").unwrap();
        f.write_all(&qual).unwrap();
        writeln!(f).unwrap();
    }
}

#[test]
fn phix_e2e_roundtrip() {
    let Some(bwa) = find_bwa_mem2() else {
        eprintln!("skip: bwa-mem2 not on PATH (set BWA_MEM2_BIN)");
        return;
    };
    if !have_samtools() {
        eprintln!("skip: samtools not on PATH");
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // 1. Write PhiX FASTA.
    let ref_fa = dir.join("phix.fa");
    {
        let mut f = fs::File::create(&ref_fa).unwrap();
        writeln!(f, ">phix").unwrap();
        for chunk in phix_seq::PHIX_SEQ.as_bytes().chunks(72) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
    }

    // 2. Build bwa-mem2 index.
    let status = Command::new(&bwa)
        .args(["index"])
        .arg(&ref_fa)
        .status()
        .expect("run bwa-mem2 index");
    assert!(status.success(), "bwa-mem2 index failed");
    assert!(
        ref_fa.with_extension("fa.bwt.2bit.64").exists(),
        "bwa-mem2 index did not produce .bwt.2bit.64"
    );

    // 3. Simulate 1000 paired reads.
    let n_pairs = 1000;
    let read_len = 100;
    let insert = 250;
    let pairs = simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), n_pairs, read_len, insert, 42);
    let r1_fq = dir.join("r1.fq");
    let r2_fq = dir.join("r2.fq");
    let r1_reads: Vec<_> = pairs
        .iter()
        .map(|(n, s, _)| (n.clone(), s.clone()))
        .collect();
    let r2_reads: Vec<_> = pairs
        .iter()
        .map(|(n, _, s)| (n.clone(), s.clone()))
        .collect();
    write_fastq(&r1_fq, &r1_reads);
    write_fastq(&r2_fq, &r2_reads);

    // 4. Run our CLI.
    let out_bam = dir.join("out.bam");
    let status = Command::new(cli_bin())
        .args(["mem"])
        .arg(ref_fa.as_os_str())
        .arg(&r1_fq)
        .arg(&r2_fq)
        .arg("-o")
        .arg(&out_bam)
        .status()
        .expect("run bwa-rs");
    assert!(status.success(), "bwa-rs mem failed");
    assert!(out_bam.exists(), "bwa-rs did not produce output BAM");

    // 5. samtools quickcheck: verifies BGZF framing + BAM magic + EOF block.
    let qc = Command::new("samtools")
        .args(["quickcheck"])
        .arg(&out_bam)
        .output()
        .unwrap();
    assert!(
        qc.status.success(),
        "samtools quickcheck failed: {}",
        String::from_utf8_lossy(&qc.stderr)
    );

    // 6. samtools view -c: record count >= 2 * n_pairs (primaries; may be
    //    higher due to secondary/supplementary).
    let view = Command::new("samtools")
        .args(["view", "-c"])
        .arg(&out_bam)
        .stdout(Stdio::piped())
        .output()
        .unwrap();
    assert!(view.status.success(), "samtools view -c failed");
    let n_recs: usize = String::from_utf8_lossy(&view.stdout)
        .trim()
        .parse()
        .expect("parse record count");
    let expected_min = 2 * n_pairs;
    assert!(
        n_recs >= expected_min,
        "too few records: got {n_recs}, expected at least {expected_min}"
    );
    eprintln!("e2e: aligned {n_pairs} pairs → {n_recs} records");
}
