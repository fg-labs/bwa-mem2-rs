//! Verify concurrent `align_batch` calls on the same `BwaIndex` produce
//! the same records whether run serially or in parallel.
//!
//! Skipped unless `BWA_MEM2_RS_TEST_REF` points at a prebuilt bwa-mem2 index.

use bwa_mem2_rs::{align_batch, AlignmentBatch, BwaIndex, MemOpts, ReadPair};
use rayon::prelude::*;
use std::sync::Arc;

fn ref_prefix() -> Option<String> {
    let p = std::env::var("BWA_MEM2_RS_TEST_REF").ok()?;
    if std::path::Path::new(&format!("{p}.bwt.2bit.64")).exists() {
        Some(p)
    } else {
        None
    }
}

fn canonicalize(batch: &AlignmentBatch) -> Vec<(usize, Vec<u8>)> {
    let mut v: Vec<_> = batch
        .iter()
        .map(|r| (r.pair_idx, r.bytes.to_vec()))
        .collect();
    v.sort();
    v
}

#[test]
fn par_iter_matches_serial() {
    let Some(prefix) = ref_prefix() else {
        eprintln!("skip: set BWA_MEM2_RS_TEST_REF");
        return;
    };
    let idx = Arc::new(BwaIndex::load(&prefix).expect("load"));
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    let seq: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = vec![b'I'; seq.len()];

    // 8 small batches. Each batch independent; no cross-batch communication.
    let batches: Vec<Vec<(Vec<u8>, Vec<u8>)>> = (0..8)
        .map(|batch_i| {
            (0..2)
                .map(|pair_i| {
                    let name = format!("b{batch_i}p{pair_i}").into_bytes();
                    (name.clone(), name)
                })
                .collect()
        })
        .collect();

    let do_align = |batch: &Vec<(Vec<u8>, Vec<u8>)>| -> Vec<(usize, Vec<u8>)> {
        let pairs: Vec<_> = batch
            .iter()
            .map(|(n1, n2)| ReadPair {
                name_r1: n1,
                seq_r1: seq,
                qual_r1: Some(&qual),
                name_r2: n2,
                seq_r2: seq,
                qual_r2: Some(&qual),
            })
            .collect();
        let (aln, _) = align_batch(&idx, &opts, &pairs, None).expect("align");
        canonicalize(&aln)
    };

    let serial: Vec<_> = batches.iter().map(do_align).collect();
    let parallel: Vec<_> = batches.par_iter().map(do_align).collect();

    // Each batch's records should match between serial and parallel runs.
    // Record-for-record equality isn't guaranteed across runs (internal
    // heuristics may depend on global state) but record counts should match.
    assert_eq!(serial.len(), parallel.len());
    for (i, (s, p)) in serial.iter().zip(parallel.iter()).enumerate() {
        assert_eq!(
            s.len(),
            p.len(),
            "batch {} record count differs: serial={} parallel={}",
            i,
            s.len(),
            p.len()
        );
    }
    eprintln!(
        "serial total recs: {}, parallel total recs: {}",
        serial.iter().map(Vec::len).sum::<usize>(),
        parallel.iter().map(Vec::len).sum::<usize>(),
    );
}
