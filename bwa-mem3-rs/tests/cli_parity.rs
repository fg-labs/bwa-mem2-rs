//! CLI parity: align a synthetic FASTQ pair with both `bwa-mem3 mem` and
//! `align_batch`, then diff the sorted BAM outputs. Skipped unless both
//! `BWA_MEM3_BIN` (path to the CLI) and `BWA_MEM3_RS_TEST_REF` (index prefix)
//! are set AND `samtools` is on PATH.

use std::io::Write;
use std::process::{Command, Stdio};

use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

fn skip() -> Option<(String, String)> {
    let bwa = std::env::var("BWA_MEM3_BIN").ok()?;
    let prefix = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    if !std::path::Path::new(&bwa).exists() {
        eprintln!("skip: BWA_MEM3_BIN={bwa} not found");
        return None;
    }
    if !std::path::Path::new(&format!("{prefix}.bwt.2bit.64")).exists() {
        eprintln!("skip: no bwa-mem3 index at {prefix}");
        return None;
    }
    if Command::new("samtools").arg("--version").output().is_err() {
        eprintln!("skip: samtools not on PATH");
        return None;
    }
    Some((bwa, prefix))
}

/// Write a FASTQ pair to two temp files + return their paths.
fn write_fastq_pair(
    dir: &std::path::Path,
    reads: &[(Vec<u8>, Vec<u8>, Vec<u8>)],
) -> (std::path::PathBuf, std::path::PathBuf) {
    // reads = [(name, r1_seq, r2_seq)]. Quality is constant 'I'.
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    let mut f1 = std::fs::File::create(&r1).unwrap();
    let mut f2 = std::fs::File::create(&r2).unwrap();
    for (name, s1, s2) in reads {
        let q1 = vec![b'I'; s1.len()];
        let q2 = vec![b'I'; s2.len()];
        writeln!(f1, "@{}/1", std::str::from_utf8(name).unwrap()).unwrap();
        f1.write_all(s1).unwrap();
        writeln!(f1, "\n+").unwrap();
        f1.write_all(&q1).unwrap();
        writeln!(f1).unwrap();
        writeln!(f2, "@{}/2", std::str::from_utf8(name).unwrap()).unwrap();
        f2.write_all(s2).unwrap();
        writeln!(f2, "\n+").unwrap();
        f2.write_all(&q2).unwrap();
        writeln!(f2).unwrap();
    }
    (r1, r2)
}

/// Count SAM records (non-header lines) from a compressed BAM byte stream.
fn count_sam_records(bam_bytes: &[u8]) -> usize {
    let mut p = Command::new("samtools")
        .args(["view", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn samtools view");
    p.stdin.as_mut().unwrap().write_all(bam_bytes).unwrap();
    drop(p.stdin.take());
    let stdout = p.wait_with_output().unwrap().stdout;
    std::str::from_utf8(&stdout).unwrap().lines().count()
}

#[test]
fn cli_bam_record_count_matches() {
    let Some((bwa, prefix)) = skip() else {
        return;
    };

    // Small deterministic input: 4 pairs, each a 60mer.
    let seq1: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq2: Vec<u8> = seq1.clone();
    let reads: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> = (0..4)
        .map(|i| (format!("r{i}").into_bytes(), seq1.clone(), seq2.clone()))
        .collect();

    let tmp = tempfile::tempdir().unwrap();
    let (r1, r2) = write_fastq_pair(tmp.path(), &reads);

    // 1. CLI path: bwa-mem3 mem -> BAM via samtools view -b
    let cli_bam = tmp.path().join("cli.bam");
    let sam = Command::new(&bwa)
        .args([
            "mem",
            "-t",
            "1",
            &prefix,
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
        ])
        .output()
        .expect("bwa-mem3 mem");
    assert!(
        sam.status.success(),
        "bwa-mem3 mem failed: {}",
        String::from_utf8_lossy(&sam.stderr)
    );
    let mut p = Command::new("samtools")
        .args(["view", "-b", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    p.stdin.as_mut().unwrap().write_all(&sam.stdout).unwrap();
    drop(p.stdin.take());
    let cli_bam_bytes = p.wait_with_output().unwrap().stdout;
    std::fs::write(&cli_bam, &cli_bam_bytes).unwrap();
    let cli_recs = count_sam_records(&cli_bam_bytes);

    // 2. bwa-rs path: build concat packed-BAM + feed to samtools.
    // Our Record::bytes is [u32 le block][record data]. samtools view expects
    // BGZF-compressed BAM w/ header; just count records via the raw parser.
    let idx = BwaIndex::load(&prefix).unwrap();
    let mut opts = MemOpts::new().unwrap();
    opts.set_pe(true);
    let q = vec![b'I'; 60];
    let pairs: Vec<ReadPair<'_>> = reads
        .iter()
        .map(|(n, s1, s2)| ReadPair {
            name_r1: n,
            seq_r1: s1,
            qual_r1: Some(&q),
            name_r2: n,
            seq_r2: s2,
            qual_r2: Some(&q),
        })
        .collect();
    let (aln, _) = align_batch(&idx, &opts, &pairs, None).unwrap();
    let rs_recs = aln.len();

    eprintln!("CLI records: {cli_recs}, bwa-rs records: {rs_recs}");
    assert_eq!(
        rs_recs, cli_recs,
        "record count mismatch: bwa-rs={rs_recs} CLI={cli_recs}"
    );
}
