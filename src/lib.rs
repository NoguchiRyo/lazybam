use pyo3::prelude::*;
mod iterator;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
/// Python 拡張モジュール `bampy`
#[pymodule(name = "_bampy")]
fn bampy(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // クラスを Python 側へ登録
    m.add_class::<iterator::BamReader>()?;
    m.add_class::<iterator::PyBamRecord>()?;

    // サンプル関数（任意）
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // モジュール docstring（任意）
    m.add("__doc__", "Rust powered BAM reader built on noodles + PyO3")?;

    // 例: numpy 必須ならここで import しておくと安心
    py.import("numpy")?;

    Ok(())
}
