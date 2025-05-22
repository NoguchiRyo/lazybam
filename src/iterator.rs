use noodles::bgzf;
use noodles::core::region::Region;
use noodles::{bam, sam};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::fs::File;
use std::str::FromStr;
use std::sync::{Arc, Mutex};

use crate::record::PyBamRecord;

#[pyclass]
pub struct BamReader {
    header: sam::Header,
    chunk_size: usize,

    /// シーケンシャル読み出し用
    reader: Option<Arc<Mutex<bam::io::reader::Reader<bgzf::io::reader::Reader<File>>>>>,

    /// region モード時に全レコードを保持
    region_records: Option<Arc<Vec<bam::Record>>>,

    /// region モード中の現在位置
    region_pos: usize,
}

#[pymethods]
impl BamReader {
    /// path, chunk_size, region を受け取るように変更
    #[new]
    #[pyo3(signature = (path, chunk_size=None, region=None))]
    fn new(path: &str, chunk_size: Option<usize>, region: Option<&str>) -> PyResult<Self> {
        let chunk_size = chunk_size.unwrap_or(1);

        if let Some(raw_region) = region {
            // ── indexed_reader で開いて領域クエリ
            let mut indexed = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            let header = indexed
                .read_header()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            // "*" は unmapped クエリ
            let records: Vec<_> = if raw_region == "*" {
                indexed
                    .query_unmapped()
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?
                    .map(|r| r.map_err(|e| e).unwrap())
                    .collect()
            } else {
                let region =
                    raw_region
                        .parse::<Region>()
                        .map_err(|e: <Region as FromStr>::Err| {
                            PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
                        })?;
                indexed
                    .query(&header, &region)
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?
                    .map(|r| r.map_err(|e| e).unwrap())
                    .collect()
            };

            Ok(BamReader {
                header,
                chunk_size,
                reader: None,
                region_records: Some(Arc::new(records)),
                region_pos: 0,
            })
        } else {
            // ── 従来のシーケンシャル読み出し
            let mut reader = bam::io::reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            let header = reader
                .read_header()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

            Ok(BamReader {
                header,
                chunk_size,
                reader: Some(Arc::new(Mutex::new(reader))),
                region_records: None,
                region_pos: 0,
            })
        }
    }

    #[getter]
    fn _header<'py>(&self, py: Python<'py>) -> PyResult<Py<PyBytes>> {
        let mut buf = Vec::new();
        let mut w = sam::io::Writer::new(&mut buf);
        w.write_header(&self.header)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyBytes::new(py, &buf).into())
    }

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

    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    /// chunk_size ごとにレコードを返す
    fn __next__(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        // --- region_records を一度だけクローンしてローカルに逃がす
        let region_opt: Option<Arc<Vec<bam::Record>>> = slf.region_records.clone();

        if let Some(records) = region_opt {
            let start = slf.region_pos;
            if start >= records.len() {
                return Ok(None);
            }
            let end = (start + slf.chunk_size).min(records.len());
            slf.region_pos = end;

            let slice = &records[start..end];
            let mut out = Vec::with_capacity(slice.len());
            for rec in slice.iter().cloned() {
                let obj: Py<PyAny> = Py::new(py, PyBamRecord::from_record(rec))
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
                    .into();
                out.push(obj);
            }
            return Ok(Some(out));
        }

        // シーケンシャルモード
        let reader_arc = slf.reader.as_ref().unwrap().clone();
        let chunk = slf.chunk_size;
        let raw_recs: Vec<bam::Record> = py.allow_threads(move || {
            let mut guard = reader_arc.lock().unwrap();
            let mut v = Vec::with_capacity(chunk);
            for _ in 0..chunk {
                let mut rec = bam::Record::default();
                match guard.read_record(&mut rec) {
                    Ok(0) => break,
                    Ok(_) => v.push(rec),
                    Err(e) => {
                        eprintln!("Error reading BAM record: {}", e);
                        break;
                    }
                }
            }
            v
        });

        if raw_recs.is_empty() {
            Ok(None)
        } else {
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
}
