use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use pyo3::prelude::*;
use std::{path::Path, path::PathBuf};

// 既存の関数をインポート
use crate::merge_bams::merge_chunks;
use crate::record::PyBamRecord;
use crate::record_buf::PyRecordBuf;
use crate::write_bams::write_chunk;

#[pyfunction]
pub fn write_chunk_py(
    py: Python<'_>, // ★ 追加
    header_bytes: Vec<u8>,
    records: Vec<PyRef<PyBamRecord>>, // PyRef → Py<...> に
    out_bam: &str,
    sort: bool,
) -> PyResult<()> {
    // ── 1. ヘッダ復元（GIL 必須） ─────────────────────────────
    let hdr_txt = std::str::from_utf8(&header_bytes)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    let header: sam::Header = hdr_txt.parse().map_err(|e: sam::header::ParseError| {
        pyo3::exceptions::PyValueError::new_err(e.to_string())
    })?;

    // ── 2. PyBamRecord → RecordBuf（GIL 必須） ───────────────
    let mut bufs = Vec::with_capacity(records.len());
    for rec in &records {
        let buf: RecordBuf = rec
            .to_record_buf()
            .expect("failed to convert PyBamRecord to RecordBuf");
        bufs.push(buf);
    }
    drop(records); // PyObject の参照を早めに解放（任意）

    // ── 3. 重い処理を GIL なしで実行 ────────────────────────
    py.allow_threads(|| write_chunk(&header, &mut bufs, out_bam, sort))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
}

#[pyfunction]
pub fn write_recordbuf_chunk_py(
    py: Python<'_>,
    header_bytes: Vec<u8>,
    records: Vec<PyRef<PyRecordBuf>>,
    out_bam: &str,
    sort: bool,
) -> PyResult<()> {
    // ── 1. ヘッダー復元 ───────────────────────────────────────────────
    let hdr_txt = std::str::from_utf8(&header_bytes)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    let header: sam::Header = hdr_txt.parse().map_err(|e: sam::header::ParseError| {
        pyo3::exceptions::PyValueError::new_err(e.to_string())
    })?;

    // ── 2. PyRecordBuf → RecordBuf に変換 ───────────────────────────────
    //    RecordBuf は Clone を実装しているので、.as_record_buf().clone() で所有権を得る
    let mut bufs: Vec<RecordBuf> = Vec::with_capacity(records.len());
    for py_rec in &records {
        // PyRecordBuf.as_record_buf() は &RecordBuf を返す
        let buf: RecordBuf = py_rec.as_record_buf().clone();
        bufs.push(buf);
    }
    // 参照が増えていたら早めに解放しておく
    drop(records);

    // ── 3. 重い I/O 処理を GIL なしで実行 ───────────────────────────────
    py.allow_threads(|| write_chunk(&header, &mut bufs, out_bam, sort))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
}

#[pyfunction]
pub fn merge_chunks_py(
    header_bytes: Vec<u8>,
    chunks: Vec<String>,
    out_bam: &str,
    sort: bool,
) -> PyResult<()> {
    // 1) Header を復元
    let hdr_str = std::str::from_utf8(&header_bytes)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
    let header: sam::Header = hdr_str.parse().map_err(|e: sam::header::ParseError| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })?;

    // 2) PathBuf に変換
    let chunk_paths: Vec<PathBuf> = chunks.into_iter().map(PathBuf::from).collect();

    // 3) merge_chunks を呼ぶ
    merge_chunks(&header, &chunk_paths, Path::new(out_bam), sort)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
}
