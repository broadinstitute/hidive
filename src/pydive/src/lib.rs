use pyo3::prelude::*;
use skydive;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn add_one(x: i32) -> PyResult<i32> {
    Ok(x + 1)
}


/// Formats the sum of two numbers as string.
#[pyfunction]
fn skytest(x: i32) -> PyResult<i32> {
    skydive::elog!("Hello, world!");
    Ok(x + 1)
}

/// A Python module implemented in Rust.
#[pymodule]
fn pydive(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add_one, m)?)?;
    m.add_function(wrap_pyfunction!(skytest, m)?)?;

    Ok(())
}