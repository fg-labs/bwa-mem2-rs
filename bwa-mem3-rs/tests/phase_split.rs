//! Exercises the separated seed / extend / estimate_pestat APIs.

use bwa_mem3_rs::{
    align_batch, estimate_pestat, extend_batch, seed_batch, BwaIndex, MemOpts, ReadPair,
};

fn ref_prefix() -> Option<String> {
    let p = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    if std::path::Path::new(&format!("{p}.bwt.2bit.64")).exists() {
        Some(p)
    } else {
        None
    }
}

fn make_pairs(seq: &'static [u8], qual: &'static [u8], n: usize) -> Vec<ReadPair<'static>> {
    (0..n)
        .map(|i| {
            let name = Box::leak(format!("r{i}").into_boxed_str()).as_bytes();
            ReadPair {
                name_r1: name,
                seq_r1: seq,
                qual_r1: Some(qual),
                name_r2: name,
                seq_r2: seq,
                qual_r2: Some(qual),
            }
        })
        .collect()
}

#[test]
fn seed_then_extend_matches_align_batch() {
    let Some(prefix) = ref_prefix() else {
        eprintln!("skip: set BWA_MEM3_RS_TEST_REF");
        return;
    };
    let idx = BwaIndex::load(&prefix).expect("load");
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    let seq: &'static [u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual: &'static [u8] = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let pairs = make_pairs(seq, qual, 2);

    // Path A: align_batch (single call).
    let (aln_a, _pes_a) = align_batch(&idx, &opts, &pairs, None).expect("align");
    let n_a = aln_a.len();

    // Path B: seed_batch + extend_batch (separated phases).
    let seeds = seed_batch(&idx, &opts, &pairs).expect("seed");
    let (aln_b, _pes_b) = extend_batch(&idx, &opts, seeds, &pairs, None).expect("extend");
    let n_b = aln_b.len();

    assert_eq!(n_a, n_b, "record count should match: A={n_a} B={n_b}");
    eprintln!("align_batch = {n_a}, seed+extend = {n_b}");
}

#[test]
fn estimate_pestat_runs_without_emitting() {
    let Some(prefix) = ref_prefix() else {
        return;
    };
    let idx = BwaIndex::load(&prefix).expect("load");
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    let seq: &'static [u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual: &'static [u8] = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    let pairs = make_pairs(seq, qual, 4);

    let _pestat = estimate_pestat(&idx, &opts, &pairs).expect("estimate");
    // Not enough well-aligning pairs in this synthetic set to populate the
    // FR orientation, so we just assert the call succeeds and returns a
    // well-formed MemPeStat.
}
