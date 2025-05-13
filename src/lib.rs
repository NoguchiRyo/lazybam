use pyo3::prelude::*;
mod iterator;

/// A Python module implemented in Rust.
#[pymodule(name = "_bampy")]
fn bampy(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<iterator::BamReader>()?;
    m.add_class::<iterator::PyBamRecord>()?;

    m.add("__doc__", "Rust powered BAM reader built on noodles + PyO3")?;

    py.import("numpy")?;

    Ok(())
}
