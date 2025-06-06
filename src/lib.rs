use pyo3::prelude::*;
mod iterator;
mod merge_bams;
mod record;
mod record_buf;
mod record_override;
mod write;
mod write_bams;

/// A Python module implemented in Rust.
#[pymodule(name = "lazybam")]
fn lazybam(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<iterator::BamReader>()?;
    m.add_class::<record::PyBamRecord>()?;
    m.add_class::<record_override::RecordOverride>()?;
    m.add_class::<record_buf::PyRecordBuf>()?;
    m.add_function(wrap_pyfunction!(write::write_chunk_py, m)?)?;
    m.add_function(wrap_pyfunction!(write::write_recordbuf_chunk_py, m)?)?;
    m.add_function(wrap_pyfunction!(write::merge_chunks_py, m)?)?;

    m.add("__doc__", "Rust powered BAM reader built on noodles + PyO3")?;

    py.import("numpy")?;

    Ok(())
}
