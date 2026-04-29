//! Python bindings for `bwa-mem3-rs`.
//!
//! Exposes the safe Rust API as a CPython extension module
//! (`bwa_mem3._bwa_mem3`). The thin facade in `python/bwa_mem3/__init__.py`
//! re-exports everything for the canonical `import bwa_mem3` namespace.

// pyo3 0.22's #[pymethods] macro injects an Into::into on the return
// type that newer clippy flags as `useless_conversion`. Allow at crate
// level rather than dotting #[allow] across every method.
#![allow(clippy::useless_conversion)]

use std::path::PathBuf;
use std::sync::Arc;

use bwa_mem3_rs as bwa;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList, PyTuple};

fn map_err(e: bwa::Error) -> PyErr {
    match e {
        bwa::Error::InvalidInput(m) => PyValueError::new_err(m),
        other => PyRuntimeError::new_err(other.to_string()),
    }
}

// ---------- BwaIndex ----------

/// Reference index handle.
///
/// `BwaIndex(prefix)` loads from disk; if a shared-memory segment for
/// `prefix` is staged, bwa-mem3 transparently attaches to it. Callers
/// that want fail-fast behavior should probe with `shm.is_staged`
/// first.
#[pyclass(name = "BwaIndex", module = "bwa_mem3._bwa_mem3", frozen)]
struct PyBwaIndex {
    inner: Arc<bwa::BwaIndex>,
}

#[pymethods]
impl PyBwaIndex {
    #[new]
    fn new(prefix: PathBuf) -> PyResult<Self> {
        let inner = bwa::BwaIndex::load(&prefix).map_err(map_err)?;
        Ok(Self {
            inner: Arc::new(inner),
        })
    }

    fn n_contigs(&self) -> usize {
        self.inner.n_contigs()
    }

    /// Returns `[(name, length), ...]` for every contig in the index.
    fn contigs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty_bound(py);
        for (name, len) in self.inner.contigs() {
            let t = PyTuple::new_bound(py, [name.into_py(py), len.into_py(py)]);
            list.append(t)?;
        }
        Ok(list)
    }

    fn __repr__(&self) -> String {
        format!("BwaIndex(n_contigs={})", self.inner.n_contigs())
    }
}

// ---------- MemOpts ----------

/// Bwa-mem3 alignment options.
///
/// Constructed with bwa-mem3 defaults; mutate via setters. Setters return
/// `None` (in-place mutation) to follow Python conventions; chain via the
/// constructor + setattr pattern if needed.
#[pyclass(name = "MemOpts", module = "bwa_mem3._bwa_mem3")]
struct PyMemOpts {
    inner: bwa::MemOpts,
}

#[pymethods]
impl PyMemOpts {
    #[new]
    fn new() -> PyResult<Self> {
        let inner = bwa::MemOpts::new().map_err(map_err)?;
        Ok(Self { inner })
    }

    /// Apply a `-x` preset on top of current values.
    /// `mode` is one of `"pacbio"`, `"ont2d"`, `"intractg"`.
    fn apply_mode(&mut self, mode: &str) -> PyResult<()> {
        let m = match mode {
            "pacbio" => bwa::Mode::Pacbio,
            "ont2d" => bwa::Mode::Ont2d,
            "intractg" => bwa::Mode::Intractg,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown mode {other:?}; expected one of 'pacbio', 'ont2d', 'intractg'"
                )))
            }
        };
        // Drop-and-replace: with_mode consumes self in Rust.
        let inner = std::mem::replace(&mut self.inner, bwa::MemOpts::new().map_err(map_err)?);
        self.inner = inner.with_mode(m);
        Ok(())
    }

    fn set_pe(&mut self, is_pe: bool) {
        self.inner.set_pe(is_pe);
    }

    fn set_soft_clip_supplementary(&mut self, v: bool) {
        self.inner.set_soft_clip_supplementary(v);
    }

    #[getter]
    fn min_seed_len(&self) -> i32 {
        self.inner.min_seed_len()
    }
    #[setter]
    fn set_min_seed_len(&mut self, v: i32) {
        self.inner.set_min_seed_len(v);
    }

    #[getter]
    fn band_width(&self) -> i32 {
        self.inner.band_width()
    }
    #[setter]
    fn set_band_width(&mut self, v: i32) {
        self.inner.set_band_width(v);
    }

    #[getter]
    fn match_score(&self) -> i32 {
        self.inner.match_score()
    }
    #[setter]
    fn set_match_score(&mut self, v: i32) {
        self.inner.set_match_score(v);
    }

    #[getter]
    fn mismatch_penalty(&self) -> i32 {
        self.inner.mismatch_penalty()
    }
    #[setter]
    fn set_mismatch_penalty(&mut self, v: i32) {
        self.inner.set_mismatch_penalty(v);
    }

    fn set_gap_open(&mut self, del: i32, ins: i32) {
        self.inner.set_gap_open(del, ins);
    }
    fn set_gap_extend(&mut self, del: i32, ins: i32) {
        self.inner.set_gap_extend(del, ins);
    }
    fn set_clip_penalty(&mut self, five: i32, three: i32) {
        self.inner.set_clip_penalty(five, three);
    }

    #[getter]
    fn minimum_score(&self) -> i32 {
        self.inner.minimum_score()
    }
    #[setter]
    fn set_minimum_score(&mut self, v: i32) {
        self.inner.set_minimum_score(v);
    }

    #[getter]
    fn max_occurrences(&self) -> i32 {
        self.inner.max_occurrences()
    }
    #[setter]
    fn set_max_occurrences(&mut self, v: i32) {
        self.inner.set_max_occurrences(v);
    }

    fn set_xa_max_hits(&mut self, primary: i32, alt: i32) {
        self.inner.set_xa_max_hits(primary, alt);
    }
    fn set_xa_drop_ratio(&mut self, v: f32) {
        self.inner.set_xa_drop_ratio(v);
    }
    fn set_unpaired_penalty(&mut self, v: i32) {
        self.inner.set_unpaired_penalty(v);
    }

    /// Set the `@RG` ID emitted as `RG:Z:` on every record. `None` clears.
    /// Note: writes to a process-wide global in bwa-mem3.
    #[pyo3(signature = (id=None))]
    fn set_read_group_id(&mut self, id: Option<&str>) -> PyResult<()> {
        self.inner.set_read_group_id(id).map_err(map_err)?;
        Ok(())
    }
}

// ---------- PeOrient / MemPeStat ----------

/// Insert-size statistics for a single orientation.
#[pyclass(name = "PeOrient", module = "bwa_mem3._bwa_mem3", get_all, set_all)]
#[derive(Clone)]
struct PyPeOrient {
    low: i32,
    high: i32,
    failed: bool,
    avg: f64,
    std: f64,
}

#[pymethods]
impl PyPeOrient {
    #[new]
    #[pyo3(signature = (low=0, high=0, failed=false, avg=0.0, std=0.0))]
    fn new(low: i32, high: i32, failed: bool, avg: f64, std: f64) -> Self {
        Self {
            low,
            high,
            failed,
            avg,
            std,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PeOrient(low={}, high={}, failed={}, avg={}, std={})",
            self.low, self.high, self.failed, self.avg, self.std
        )
    }
}

impl From<bwa::PeOrient> for PyPeOrient {
    fn from(p: bwa::PeOrient) -> Self {
        Self {
            low: p.low,
            high: p.high,
            failed: p.failed,
            avg: p.avg,
            std: p.std,
        }
    }
}

impl From<PyPeOrient> for bwa::PeOrient {
    fn from(p: PyPeOrient) -> Self {
        Self {
            low: p.low,
            high: p.high,
            failed: p.failed,
            avg: p.avg,
            std: p.std,
        }
    }
}

fn orientation_from_str(s: &str) -> PyResult<bwa::PeOrientation> {
    match s.to_ascii_uppercase().as_str() {
        "FF" => Ok(bwa::PeOrientation::Ff),
        "FR" => Ok(bwa::PeOrientation::Fr),
        "RF" => Ok(bwa::PeOrientation::Rf),
        "RR" => Ok(bwa::PeOrientation::Rr),
        other => Err(PyValueError::new_err(format!(
            "unknown PE orientation {other:?}; expected one of 'FF', 'FR', 'RF', 'RR'"
        ))),
    }
}

/// 4-orientation paired-end insert-size model (`mem_pestat_t[4]`).
#[pyclass(name = "MemPeStat", module = "bwa_mem3._bwa_mem3")]
struct PyMemPeStat {
    inner: bwa::MemPeStat,
}

#[pymethods]
impl PyMemPeStat {
    #[new]
    fn zero() -> PyResult<Self> {
        Ok(Self {
            inner: bwa::MemPeStat::zero().map_err(map_err)?,
        })
    }

    /// `o` is one of `"FF"`, `"FR"`, `"RF"`, `"RR"`.
    fn orientation(&self, o: &str) -> PyResult<PyPeOrient> {
        let p = self.inner.orientation(orientation_from_str(o)?);
        Ok(p.into())
    }

    fn set_orientation(&mut self, o: &str, v: PyPeOrient) -> PyResult<()> {
        self.inner
            .set_orientation(orientation_from_str(o)?, v.into());
        Ok(())
    }
}

// ---------- ReadPair / Record ----------

/// One paired-end read. Holds owning references to the input bytes so
/// they remain valid while align_batch borrows them.
#[pyclass(name = "ReadPair", module = "bwa_mem3._bwa_mem3", frozen)]
struct PyReadPair {
    name_r1: Py<PyBytes>,
    seq_r1: Py<PyBytes>,
    qual_r1: Option<Py<PyBytes>>,
    name_r2: Py<PyBytes>,
    seq_r2: Py<PyBytes>,
    qual_r2: Option<Py<PyBytes>>,
}

#[pymethods]
impl PyReadPair {
    #[new]
    #[pyo3(signature = (name_r1, seq_r1, name_r2, seq_r2, qual_r1=None, qual_r2=None))]
    fn new(
        name_r1: Py<PyBytes>,
        seq_r1: Py<PyBytes>,
        name_r2: Py<PyBytes>,
        seq_r2: Py<PyBytes>,
        qual_r1: Option<Py<PyBytes>>,
        qual_r2: Option<Py<PyBytes>>,
    ) -> Self {
        Self {
            name_r1,
            seq_r1,
            qual_r1,
            name_r2,
            seq_r2,
            qual_r2,
        }
    }
}

/// One packed BAM record produced by alignment.
#[pyclass(name = "Record", module = "bwa_mem3._bwa_mem3", frozen, get_all)]
struct PyRecord {
    pair_idx: usize,
    bytes: Py<PyBytes>,
}

#[pymethods]
impl PyRecord {
    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "Record(pair_idx={}, bytes=<{} bytes>)",
            self.pair_idx,
            self.bytes.bind(py).len().unwrap_or(0)
        )
    }
}

// ---------- Alignment functions ----------

/// Borrow byte slices from a list of `PyReadPair` while the GIL is held.
/// The returned `Vec<bwa::ReadPair<'a>>` is only valid for the lifetime
/// of `pairs` and the GIL.
///
/// # Safety invariant
///
/// The returned slices alias `Py<PyBytes>` buffers owned by the
/// `PyReadPair`s in `pairs`. They are sound only while *both* (a) the
/// GIL is held (so no concurrent Python code can mutate / drop the
/// `bytes` objects) and (b) every `PyRef` in `pairs` is still alive
/// (which keeps each `Py<PyBytes>` GC-rooted). Callers must not
/// release the GIL via `py.allow_threads(...)` while holding any of
/// the returned slices, and must not store them past the function
/// scope where `pairs` lives.
fn borrow_pairs<'a>(py: Python<'a>, pairs: &'a [PyRef<'a, PyReadPair>]) -> Vec<bwa::ReadPair<'a>> {
    pairs
        .iter()
        .map(|p| bwa::ReadPair {
            name_r1: p.name_r1.bind(py).as_bytes(),
            seq_r1: p.seq_r1.bind(py).as_bytes(),
            qual_r1: p.qual_r1.as_ref().map(|q| q.bind(py).as_bytes()),
            name_r2: p.name_r2.bind(py).as_bytes(),
            seq_r2: p.seq_r2.bind(py).as_bytes(),
            qual_r2: p.qual_r2.as_ref().map(|q| q.bind(py).as_bytes()),
        })
        .collect()
}

fn batch_to_records<'py>(
    py: Python<'py>,
    batch: &bwa::AlignmentBatch,
) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty_bound(py);
    for rec in batch.iter() {
        let bytes = PyBytes::new_bound(py, rec.bytes).unbind();
        let py_rec = PyRecord {
            pair_idx: rec.pair_idx,
            bytes,
        };
        list.append(Bound::new(py, py_rec)?)?;
    }
    Ok(list)
}

/// Align a batch of read pairs. Returns `(records, pestat_out)`.
///
/// **Holds the GIL for the duration of the alignment**, so threading the
/// same `BwaIndex` from `concurrent.futures.ThreadPoolExecutor` will
/// serialize. Use `multiprocessing` (with [`shm.stage`] to share the
/// index across workers) for true parallelism. A future revision may
/// release the GIL once the input bytes are copied or pinned.
#[pyfunction]
#[pyo3(signature = (index, opts, pairs, pestat_in=None))]
fn align_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
    pestat_in: Option<&PyMemPeStat>,
) -> PyResult<(Bound<'py, PyList>, PyMemPeStat)> {
    let borrowed = borrow_pairs(py, &pairs);
    let pestat_ref = pestat_in.map(|p| &p.inner);
    let (batch, pestat_out) =
        bwa::align_batch(&index.inner, &opts.inner, &borrowed, pestat_ref).map_err(map_err)?;
    let records = batch_to_records(py, &batch)?;
    Ok((records, PyMemPeStat { inner: pestat_out }))
}

/// Phase 1: seed only. Returns an opaque handle consumed by `extend_batch`.
#[pyfunction]
fn seed_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
) -> PyResult<PySeeds> {
    let borrowed = borrow_pairs(py, &pairs);
    let seeds = bwa::seed_batch(&index.inner, &opts.inner, &borrowed).map_err(map_err)?;
    Ok(PySeeds { inner: Some(seeds) })
}

/// Phase 2: extend pre-computed seeds to alignments. Consumes `seeds`.
#[pyfunction]
#[pyo3(signature = (index, opts, seeds, pairs, pestat_in=None))]
fn extend_batch<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    seeds: &mut PySeeds,
    pairs: Vec<PyRef<'py, PyReadPair>>,
    pestat_in: Option<&PyMemPeStat>,
) -> PyResult<(Bound<'py, PyList>, PyMemPeStat)> {
    let inner = seeds
        .inner
        .take()
        .ok_or_else(|| PyValueError::new_err("seeds already consumed"))?;
    let borrowed = borrow_pairs(py, &pairs);
    let pestat_ref = pestat_in.map(|p| &p.inner);
    let (batch, pestat_out) =
        bwa::extend_batch(&index.inner, &opts.inner, inner, &borrowed, pestat_ref)
            .map_err(map_err)?;
    let records = batch_to_records(py, &batch)?;
    Ok((records, PyMemPeStat { inner: pestat_out }))
}

/// Run seeding + SE-extension on a pilot batch and return the inferred
/// 4-orientation insert-size model. Discards the alignments.
#[pyfunction]
fn estimate_pestat<'py>(
    py: Python<'py>,
    index: &PyBwaIndex,
    opts: &PyMemOpts,
    pairs: Vec<PyRef<'py, PyReadPair>>,
) -> PyResult<PyMemPeStat> {
    let borrowed = borrow_pairs(py, &pairs);
    let pestat = bwa::estimate_pestat(&index.inner, &opts.inner, &borrowed).map_err(map_err)?;
    Ok(PyMemPeStat { inner: pestat })
}

/// Opaque seeds handle from `seed_batch`. Single-shot — consumed by `extend_batch`.
#[pyclass(name = "Seeds", module = "bwa_mem3._bwa_mem3")]
struct PySeeds {
    inner: Option<bwa::Seeds>,
}

#[pymethods]
impl PySeeds {
    fn __repr__(&self) -> &'static str {
        if self.inner.is_some() {
            "Seeds(<live>)"
        } else {
            "Seeds(<consumed>)"
        }
    }
}

// ---------- shm submodule ----------

#[pyfunction(name = "is_staged")]
fn shm_is_staged(prefix: PathBuf) -> PyResult<bool> {
    bwa::shm::is_staged(&prefix).map_err(map_err)
}

#[pyfunction(name = "stage")]
fn shm_stage(prefix: PathBuf) -> PyResult<()> {
    bwa::shm::stage(&prefix).map_err(map_err)
}

#[pyfunction(name = "destroy")]
fn shm_destroy() -> PyResult<()> {
    bwa::shm::destroy().map_err(map_err)
}

#[pyfunction(name = "list")]
fn shm_list() -> PyResult<()> {
    bwa::shm::list().map_err(map_err)
}

// ---------- module init ----------

#[pymodule]
fn _bwa_mem3(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBwaIndex>()?;
    m.add_class::<PyMemOpts>()?;
    m.add_class::<PyMemPeStat>()?;
    m.add_class::<PyPeOrient>()?;
    m.add_class::<PyReadPair>()?;
    m.add_class::<PyRecord>()?;
    m.add_class::<PySeeds>()?;
    m.add_function(wrap_pyfunction!(align_batch, m)?)?;
    m.add_function(wrap_pyfunction!(seed_batch, m)?)?;
    m.add_function(wrap_pyfunction!(extend_batch, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_pestat, m)?)?;

    let shm = PyModule::new_bound(m.py(), "shm")?;
    shm.add_function(wrap_pyfunction!(shm_is_staged, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_stage, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_destroy, &shm)?)?;
    shm.add_function(wrap_pyfunction!(shm_list, &shm)?)?;
    m.add_submodule(&shm)?;
    Ok(())
}
