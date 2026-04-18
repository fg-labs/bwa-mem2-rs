//! Reference index handle.

use std::ffi::{CStr, CString};
use std::path::Path;

use crate::error::{shim_err, Error, Result};

/// Handle to a loaded bwa-mem2 reference index.
///
/// The index is immutable after loading. Multiple threads may share a
/// `BwaIndex` via `&BwaIndex` or `Arc<BwaIndex>` for concurrent alignment.
pub struct BwaIndex {
    handle: *mut bwa_mem2_sys::BwaIndex,
}

impl BwaIndex {
    /// Load a prebuilt bwa-mem2 index.
    ///
    /// `prefix` is the path without extension; bwa-mem2 appends its own
    /// suffixes (`.bwt.2bit.64`, `.ann`, `.pac`, etc) internally.
    pub fn load(prefix: impl AsRef<Path>) -> Result<Self> {
        let path = prefix.as_ref();
        let s = path
            .to_str()
            .ok_or_else(|| Error::InvalidInput("prefix must be valid UTF-8".into()))?;
        let c = CString::new(s)?;
        let handle = unsafe { bwa_mem2_sys::bwa_shim_idx_load(c.as_ptr()) };
        if handle.is_null() {
            return Err(Error::IndexLoad {
                path: path.to_owned(),
                msg: shim_err("idx load").to_string(),
            });
        }
        Ok(BwaIndex { handle })
    }

    #[must_use]
    pub fn n_contigs(&self) -> usize {
        unsafe { bwa_mem2_sys::bwa_shim_idx_n_contigs(self.handle) }
    }

    #[must_use]
    pub fn contig_name(&self, i: usize) -> &str {
        unsafe {
            let c = bwa_mem2_sys::bwa_shim_idx_contig_name(self.handle, i);
            if c.is_null() {
                ""
            } else {
                CStr::from_ptr(c).to_str().unwrap_or("")
            }
        }
    }

    #[must_use]
    pub fn contig_len(&self, i: usize) -> i64 {
        unsafe { bwa_mem2_sys::bwa_shim_idx_contig_len(self.handle, i) }
    }

    pub fn contigs(&self) -> impl Iterator<Item = (&str, i64)> + '_ {
        (0..self.n_contigs()).map(move |i| (self.contig_name(i), self.contig_len(i)))
    }

    pub(crate) fn raw(&self) -> *mut bwa_mem2_sys::BwaIndex {
        self.handle
    }
}

impl Drop for BwaIndex {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { bwa_mem2_sys::bwa_shim_idx_free(self.handle) };
        }
    }
}

// SAFETY: the shim guarantees BwaIndex is read-only after construction; no
// mutable state is exposed via any accessor. Multiple threads may share a
// `&BwaIndex` or an `Arc<BwaIndex>` for concurrent alignment.
unsafe impl Send for BwaIndex {}
unsafe impl Sync for BwaIndex {}
