use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::alignment::record::data::field::value::Array;
use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::IntoPyObjectExt;
use sam::alignment::record::data::field::Value;
use std::fs::File;

#[pyclass]
pub struct PyBamRecord {
    qname: String,
    flag: u16,
    pos: i64,
    mapq: u8,
    len: usize,

    #[pyo3(get)]
    pub tags: Vec<(String, PyObject)>, // list[(tag, value)]

    seq: Vec<u8>,
    qual: Vec<u8>,
}

#[pymethods]
impl PyBamRecord {
    #[getter]
    fn qname(&self) -> &str {
        &self.qname
    }
    #[getter]
    fn flag(&self) -> u16 {
        self.flag
    }
    #[getter]
    fn pos(&self) -> i64 {
        self.pos
    }
    #[getter]
    fn len(&self) -> usize {
        self.len
    }
    #[getter]
    fn mapq(&self) -> u8 {
        self.mapq
    }

    #[getter]
    fn seq<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &self.seq))
    }

    #[getter]
    fn qual<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &self.qual))
    }
}

#[pyclass]
pub struct BamReader {
    reader: bam::io::Reader<bgzf::io::reader::Reader<File>>,
    header: sam::Header,
}

#[pymethods]
impl BamReader {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;

        // header を最初に読む (noodles の Reader は read_header 必須)
        let header = reader
            .read_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;

        Ok(Self { reader, header })
    }

    #[getter]
    fn header<'py>(&self, py: Python<'py>) -> PyResult<Py<PyBytes>> {
        // ヘッダを SAM テキスト化
        let mut buf = Vec::new();
        {
            let mut w = sam::io::Writer::new(&mut buf);
            w.write_header(&self.header)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        }

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

    fn __next__(mut slf: PyRefMut<'_, Self>, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        // --- レコード読み込み -----------------------------------------------------
        let mut record = bam::Record::default();
        let bytes = slf
            .reader
            .read_record(&mut record)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
        if bytes == 0 {
            return Ok(None); // ← EOF で StopIteration
        }

        // --- 基本フィールド -------------------------------------------------------
        let qname = record.name().map(|b| b.to_string()).unwrap_or_default();

        let flag: u16 = u16::from(record.flags());

        // alignment_start : Option<io::Result<Position>>
        let pos: i64 = match record.alignment_start() {
            Some(Ok(p)) => usize::from(p) as i64, // Position → usize → i64
            _ => -1,                              // 未定義なら -1
        };

        let mapq: u8 = record
            .mapping_quality()
            .map(|mq| u8::from(mq))
            .unwrap_or(255);

        let tlen: i64 = i64::from(record.template_length());

        // シーケンスは表示用に文字列、クオリティは bytes で返す
        let seq = record.sequence().as_ref().to_vec();
        let qual = record.quality_scores().as_ref().to_vec();

        // --- タグ（B 型→NumPy 配列に変換）----------------------------------------
        let mut tags = Vec::new();
        for field in record.data().iter() {
            use noodles::sam::alignment::record::data::field::value::array::Values;
            fn to_vec<T: Copy>(values: &dyn Values<'_, T>) -> Vec<T> {
                // iter() は Iterator<Item = Result<T, _>>
                values.iter().filter_map(|r| r.ok()).collect()
            }
            if let Ok((k, v)) = field {
                let tag = String::from_utf8_lossy(k.as_ref()).into_owned();
                let py_val = match v {
                    // 整数 ─────────────────────────────
                    Value::Int8(n) => (n as i32).into_py_any(py).unwrap(),
                    Value::UInt8(n) => (n as u32).into_py_any(py).unwrap(),
                    Value::Int16(n) => (n as i32).into_py_any(py).unwrap(),
                    Value::UInt16(n) => (n as u32).into_py_any(py).unwrap(),
                    Value::Int32(n) => (n as i32).into_py_any(py).unwrap(),
                    Value::UInt32(n) => (n as u32).into_py_any(py).unwrap(),

                    // 浮動小数点
                    Value::Float(f) => (f as f64).into_py_any(py).unwrap(),

                    // 1 文字
                    Value::Character(c) => (c).to_string().into_py_any(py).unwrap(),

                    // 文字列
                    Value::String(s) => String::from_utf8_lossy(s)
                        .into_owned()
                        .into_py_any(py)
                        .unwrap(),

                    // BAM 配列（B 型）→ NumPy
                    Value::Array(arr) => {
                        match arr {
                            Array::UInt8(a) => {
                                let v = to_vec::<u8>(a.as_ref()); // &dyn Values
                                PyArray1::from_vec(py, v).into_py_any(py).unwrap()
                                // → Py<PyAny>
                            }
                            Array::Int16(a) => {
                                let v = to_vec::<i16>(a.as_ref());
                                PyArray1::from_vec(py, v).into_py_any(py).unwrap()
                            }
                            Array::Float(a) => {
                                let v = to_vec::<f32>(a.as_ref());
                                PyArray1::from_vec(py, v).into_py_any(py).unwrap()
                            }
                            _ => py.None().into_py_any(py).unwrap(), // 未対応型
                        }
                    }

                    // そのほかはスキップ
                    _ => continue,
                };
                tags.push((tag, py_val));
            }
        }

        // --- Python オブジェクトを生成 ------------------------------------------
        let py_rec = PyBamRecord {
            qname,
            flag,
            pos,
            mapq,
            len: tlen as usize,
            tags,
            seq,
            qual,
        }
        .into_py_any(py)
        .unwrap();

        Ok(Some(py_rec))
    }
}
