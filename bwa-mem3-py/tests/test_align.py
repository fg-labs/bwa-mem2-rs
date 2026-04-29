"""End-to-end alignment test using a real bwa-mem3 index.

Skipped unless ``BWA_MEM3_RS_TEST_REF`` points at a prebuilt index prefix
(such that ``<prefix>.bwt.2bit.64`` exists).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from bwa_mem3 import (
    BwaIndex,
    MemOpts,
    ReadPair,
    align_batch,
    estimate_pestat,
    extend_batch,
    seed_batch,
)


def _ref_prefix() -> str | None:
    p = os.environ.get("BWA_MEM3_RS_TEST_REF")
    if p is None:
        return None
    if not Path(f"{p}.bwt.2bit.64").exists():
        return None
    return p


pytestmark = pytest.mark.skipif(
    _ref_prefix() is None,
    reason="BWA_MEM3_RS_TEST_REF not set or index missing",
)


def _synthetic_pair(name: bytes) -> ReadPair:
    seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    qual = b"I" * len(seq)
    return ReadPair(name_r1=name, seq_r1=seq, qual_r1=qual,
                    name_r2=name, seq_r2=seq, qual_r2=qual)


def test_load_index_and_align() -> None:
    prefix = _ref_prefix()
    assert prefix is not None
    idx = BwaIndex(prefix)
    assert idx.n_contigs() > 0
    contigs = idx.contigs()
    assert all(isinstance(c, tuple) and len(c) == 2 for c in contigs)
    assert all(isinstance(c[0], str) and isinstance(c[1], int) for c in contigs)

    opts = MemOpts()
    opts.set_pe(True)

    pairs = [_synthetic_pair(b"r0"), _synthetic_pair(b"r1")]
    records, _pestat = align_batch(idx, opts, pairs)
    assert len(records) >= 2 * len(pairs)
    for rec in records:
        assert isinstance(rec.bytes, bytes)
        assert rec.pair_idx in (0, 1)


def test_phase_split_matches_align_batch() -> None:
    prefix = _ref_prefix()
    assert prefix is not None
    idx = BwaIndex(prefix)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(b"x")]

    seeds = seed_batch(idx, opts, pairs)
    rec_phased, _ = extend_batch(idx, opts, seeds, pairs)

    rec_one_shot, _ = align_batch(idx, opts, pairs)

    # Same record count and same record bytes per pair (orderings match
    # because both paths run a single thread serially).
    assert len(rec_phased) == len(rec_one_shot)
    for a, b in zip(rec_phased, rec_one_shot):
        assert a.pair_idx == b.pair_idx
        assert a.bytes == b.bytes


def test_seeds_consumed_only_once() -> None:
    """extend_batch consumes the Seeds; reusing it raises ValueError."""
    prefix = _ref_prefix()
    assert prefix is not None
    idx = BwaIndex(prefix)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(b"x")]

    seeds = seed_batch(idx, opts, pairs)
    extend_batch(idx, opts, seeds, pairs)
    with pytest.raises(ValueError):
        extend_batch(idx, opts, seeds, pairs)


def test_estimate_pestat_then_align_with_pestat_in() -> None:
    """estimate_pestat -> align_batch with pestat_in round-trip."""
    prefix = _ref_prefix()
    assert prefix is not None
    idx = BwaIndex(prefix)
    opts = MemOpts()
    opts.set_pe(True)
    pairs = [_synthetic_pair(f"r{i}".encode()) for i in range(4)]

    pestat = estimate_pestat(idx, opts, pairs)
    records, _out_pestat = align_batch(idx, opts, pairs, pestat_in=pestat)
    assert len(records) >= 2 * len(pairs)
