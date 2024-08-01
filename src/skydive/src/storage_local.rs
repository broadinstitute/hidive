use anyhow::Result;
use chrono::{DateTime, Utc};

use std::fs::metadata;
use std::path::PathBuf;

pub fn local_get_file_update_time(path: &PathBuf) -> Result<DateTime<Utc>> {
    let metadata = metadata(path)?;
    let modified_time = metadata.modified()?;

    Ok(DateTime::<Utc>::from(modified_time))
}
