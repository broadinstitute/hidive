use skydive;
use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn add_one(x: i32) -> PyResult<i32> {
    Ok(x + 1)
}

/// A Python module implemented in Rust.
#[pymodule]
fn pydive(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add_one, m)?)?;

    Ok(())
}