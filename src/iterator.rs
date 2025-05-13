use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::alignment::record::data::field::value::Array;
use numpy::{ndarray::Array2, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::IntoPyObjectExt;
use sam::alignment::record::data::field::Value;
use std::fs::File;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Instant;

// 内部で扱う Rust レコード表現
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum PyKind {
    Match = 0,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

impl From<sam::alignment::record::cigar::op::Kind> for PyKind {
    fn from(k: sam::alignment::record::cigar::op::Kind) -> Self {
        use sam::alignment::record::cigar::op::Kind::*;
        match k {
            Match => PyKind::Match,
            Insertion => PyKind::Insertion,
            Deletion => PyKind::Deletion,
            Skip => PyKind::Skip,
            SoftClip => PyKind::SoftClip,
            HardClip => PyKind::HardClip,
            Pad => PyKind::Pad,
            SequenceMatch => PyKind::SequenceMatch,
            SequenceMismatch => PyKind::SequenceMismatch,
        }
    }
}
/// Python 側でラップする BAM レコード
/// Python 側でラップする BAM レコード
#[pyclass]
pub struct PyBamRecord {
    record: bam::Record,
}

/// internal constructor: wrap bam::Record without exposing to Python
impl PyBamRecord {
    fn from_record(record: bam::Record) -> Self {
        PyBamRecord { record }
    }
}
#[pymethods]
impl PyBamRecord {
    #[getter]
    fn qname(&self) -> String {
        let qname = self
            .record
            .name()
            .map(|b| b.to_string())
            .unwrap_or_default();
        return qname;
    }

    #[getter]
    fn flag(&self) -> u16 {
        u16::from(self.record.flags())
    }

    #[getter]
    fn pos(&self) -> i64 {
        self.record
            .alignment_start()
            .and_then(|r| r.ok())
            .map(|p| usize::from(p) as i64)
            .unwrap_or(-1)
    }

    #[getter]
    fn mapq(&self) -> u8 {
        self.record
            .mapping_quality()
            .map(|mq| u8::from(mq))
            .unwrap_or(255)
    }

    #[getter]
    fn len(&self) -> usize {
        self.record.template_length().abs() as usize
    }

    #[getter]
    fn seq(&self) -> String {
        let seq: String = self.record.sequence().iter().map(|b| b as char).collect();
        return seq;
    }

    #[getter]
    fn qual<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        let qual = self.record.quality_scores().as_ref().to_vec();
        Ok(PyArray1::from_vec(py, qual))
    }

    #[getter]
    fn cigar<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let mut cigar_ops = Vec::new();
        for op in self.record.cigar().iter() {
            let py_op = match op {
                Ok(op) => (PyKind::from(op.kind()), op.len() as usize),
                Err(_) => (PyKind::Match, 0), // エラー時は Match で len=0
            };
            cigar_ops.push(py_op);
        }
        let n = cigar_ops.len();
        // len, code を交互に格納する flat ベクタ
        let mut buf: Vec<u32> = Vec::with_capacity(n * 2);
        for (kind, len) in &cigar_ops {
            buf.push(*kind as u32);
            buf.push(*len as u32);
        }
        // shape = (n, 2)
        Ok(PyArray2::from_owned_array(
            py,
            Array2::from_shape_vec((n, 2), buf).unwrap(),
        ))
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> Vec<(String, PyObject)> {
        let mut tags = Vec::new();
        for field in self.record.data().iter().filter_map(Result::ok) {
            let key = String::from_utf8_lossy(field.0.as_ref()).into_owned();
            let val = match field.1 {
                Value::Int8(n) => (n as i32).into_py_any(py).unwrap(),
                Value::UInt8(n) => (n as u32).into_py_any(py).unwrap(),
                Value::Int16(n) => (n as i32).into_py_any(py).unwrap(),
                Value::UInt16(n) => (n as u32).into_py_any(py).unwrap(),
                Value::Int32(n) => (n as i32).into_py_any(py).unwrap(),
                Value::UInt32(n) => (n as u32).into_py_any(py).unwrap(),
                Value::Float(f) => (f as f64).into_py_any(py).unwrap(),
                Value::Character(c) => c.to_string().into_py_any(py).unwrap(),
                Value::String(bs) => String::from_utf8_lossy(bs)
                    .into_owned()
                    .into_py_any(py)
                    .unwrap(),
                Value::Array(arr) => match arr {
                    Array::UInt8(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Int8(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Int16(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Float(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    _ => py.None().into_py_any(py).unwrap(),
                },
                _ => continue,
            };
            tags.push((key, val));
        }
        tags
    }
}
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
    fn header<'py>(&self, py: Python<'py>) -> PyResult<Py<PyBytes>> {
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
        // Reader は Drop で閉じるので特に不要だが、明示的にファイル flush/close するならここ
        Ok(())
    }

    /// イテレータ
    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    fn __next__(slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        let t0 = Instant::now();

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

        let dur_io = t0.elapsed();

        let t1 = Instant::now();

        let dur_prep = t1.elapsed();
        let t2 = Instant::now();
        // GIL 下で PyBamRecord を生の Record から直接ラップ
        let mut out = Vec::with_capacity(raw_recs.len());
        for rec in raw_recs {
            // Py::new で直接インスタンス化
            let obj: Py<PyAny> = Py::new(py, PyBamRecord::from_record(rec))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
                .into();
            out.push(obj);
        }
        let dur_conv = t2.elapsed();
        eprintln!(
            "__next__ timings: IO={:?}, prep={:?}, conv={:?}, total={:?}",
            dur_io,
            dur_prep,
            dur_conv,
            dur_io + dur_prep + dur_conv
        );

        Ok(Some(out))
    }
}
