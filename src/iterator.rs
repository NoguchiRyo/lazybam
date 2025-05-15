use noodles::bgzf;
use noodles::{bam, sam};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::fs::File;
use std::sync::Arc;
use std::sync::Mutex;

use crate::record::PyBamRecord;

#[pyclass]
pub struct BamReader {
    reader: Arc<Mutex<bam::io::Reader<bgzf::io::reader::Reader<File>>>>,
    header: sam::Header,
    chunk_size: usize,
}

#[pymethods]
impl BamReader {
    #[new]
    fn new(path: &str, chunk_size: Option<usize>) -> PyResult<Self> {
        let chunk_size = chunk_size.unwrap_or(1);
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let header = reader
            .read_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(Self {
            reader: Arc::new(Mutex::new(reader)),
            header,
            chunk_size,
        })
    }

    #[getter]
    fn _header<'py>(&self, py: Python<'py>) -> PyResult<Py<PyBytes>> {
        // ヘッダを SAM テキスト化
        let mut buf = Vec::new();
        let mut w = sam::io::Writer::new(&mut buf);
        w.write_header(&self.header)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Bound<'py, PyBytes> → Py<PyBytes>
        Ok(PyBytes::new(py, &buf).into()) // ここを .into()
    }

    /// context manager
    fn __enter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }
    fn __exit__(
        _slf: PyRefMut<'_, Self>,
        _exc_type: PyObject,
        _exc_val: PyObject,
        _trace: PyObject,
    ) -> PyResult<()> {
        Ok(())
    }

    /// イテレータ
    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        let chunk_size = slf.chunk_size;

        let reader_clone = Arc::clone(&slf.reader);

        let raw_recs: Vec<bam::Record> = py.allow_threads(move || {
            let mut guard = reader_clone.lock().unwrap();
            let mut v = Vec::with_capacity(chunk_size);

            for _ in 0..chunk_size {
                let mut rec = bam::Record::default();
                match guard.read_record(&mut rec) {
                    Ok(0) => break, // EOF
                    Ok(_) => v.push(rec),
                    Err(e) => {
                        // エラー時は空の Vec を返す
                        eprintln!("Error reading BAM record: {}", e);
                        return vec![];
                    }
                }
            }
            v
        });

        if raw_recs.is_empty() {
            return Ok(None);
        }

        let mut out = Vec::with_capacity(raw_recs.len());
        for rec in raw_recs {
            let obj: Py<PyAny> = Py::new(py, PyBamRecord::from_record(rec))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
                .into();
            out.push(obj);
        }

        Ok(Some(out))
    }
}
