use anyhow::Result;
use chrono::{DateTime, Utc};

use std::fs::metadata;
use std::path::PathBuf;

/// Get the update time of a file in the local filesystem.
///
/// # Arguments
///
/// * `path` - A reference to a `PathBuf` object representing the path to the file.
///
/// # Returns
///
/// A `DateTime<Utc>` object representing the update time of the file.
///
/// # Errors
///
/// This function returns an error if the file metadata cannot be accessed.
///
/// # Panics
///
/// This function panics if the file path is invalid.
pub fn local_get_file_update_time(path: &PathBuf) -> Result<DateTime<Utc>> {
    let metadata = metadata(path)?;
    let modified_time = metadata.modified()?;

    Ok(DateTime::<Utc>::from(modified_time))
}
