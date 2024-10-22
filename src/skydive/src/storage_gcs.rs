use anyhow::Result;

use chrono::{DateTime, Utc};
use cloud_storage::{object::ObjectList, sync::Client, ListRequest};

/// Split a GCS path into a bucket name and a prefix.
/// The GCS path should be in the format `gs://bucket_name/prefix`.
///
/// # Arguments
///
/// * `path` - A string slice that holds the GCS path.
///
/// # Returns
///
/// A tuple with the bucket name and the prefix.
///
/// # Panics
///
/// This function panics if the GCS path is invalid.
#[must_use]
pub fn gcs_split_path(path: &str) -> (String, String) {
    let re = regex::Regex::new(r"^gs://").unwrap();
    let path = re.replace(path, "");
    let split: Vec<&str> = path.split('/').collect();

    let bucket_name = split[0].to_string();
    let prefix = split[1..].join("/");

    (bucket_name, prefix)
}
/// List all files in a GCS path.
/// The GCS path should be in the format `gs://bucket_name/prefix`.
///
/// # Arguments
///
/// * `path` - A string slice that holds the GCS path.
///
/// # Returns
///
/// A vector of `ObjectList` objects representing the files in the GCS path.
///
/// # Errors
///
/// This function returns an error if the GCS client cannot be created or the object list cannot be read.
///
/// # Panics
///
/// This function panics if the GCS path is invalid.
pub fn gcs_list_files(path: &str) -> Result<Vec<ObjectList>> {
    let (bucket_name, prefix) = gcs_split_path(path);

    let client = Client::new()?;
    let file_list = client.object().list(
        &bucket_name,
        ListRequest {
            prefix: Some(prefix),
            ..Default::default()
        },
    )?;

    Ok(file_list)
}

/// Get the update time of a file in GCS.
/// The GCS path should be in the format `gs://bucket_name/prefix`.
///
/// # Arguments
///
/// * `path` - A string slice that holds the GCS path.
///
/// # Returns
///
/// A `DateTime<Utc>` object representing the update time of the file.
///
/// # Errors
///
/// This function returns an error if the GCS client cannot be created or the object cannot be read.
///
/// # Panics
///
/// This function panics if the GCS path is invalid.
pub fn gcs_get_file_update_time(path: &str) -> Result<DateTime<Utc>> {
    let (bucket_name, prefix) = gcs_split_path(path);

    let client = Client::new()?;
    let object = client.object().read(&bucket_name, &prefix)?;

    Ok(object.updated)
}

/// Download a file from GCS and return the local filename.
/// The GCS path should be in the format `gs://bucket_name/prefix`.
///
/// # Arguments
///
/// * `path` - A string slice that holds the GCS path.
///
/// # Returns
///
/// A string with the local filename.
///
/// # Errors
///
/// This function returns an error if the GCS client cannot be created
/// or the object cannot be downloaded.
///
/// # Panics
///
/// This function panics if the GCS path is invalid.
pub fn gcs_download_file(path: &str) -> Result<String> {
    let (bucket_name, prefix) = gcs_split_path(&path);
    let filename = prefix.split('/').last().unwrap_or_default().to_string();

    if !std::path::Path::new(&filename).exists() {
        let client = Client::new().unwrap();
        let bytes = client.object().download(&bucket_name, &prefix).unwrap();

        std::fs::write(&filename, bytes)?;
    }

    Ok(filename)
}

/// List all files in a GCS path with a specific suffix.
///
/// # Arguments
///
/// * `path` - A string slice that holds the GCS path.
/// * `suffix` - A string slice that holds the suffix to filter the files.
///
/// # Returns
///
/// A vector of strings with the names of the files that match the suffix.
///
/// # Errors
///
/// This function returns an error if the GCS client cannot be created or the object cannot be read.
///
/// # Panics
///
/// This function panics if the GCS path is invalid.
pub fn gcs_list_files_of_type(path: &str, suffix: &str) -> Result<Vec<String>> {
    let file_list = gcs_list_files(&path).unwrap();

    let bam_files: Vec<_> = file_list
        .iter()
        .flat_map(|fs| {
            fs.items
                .iter()
                .filter_map(|f| {
                    if f.name.ends_with(suffix) {
                        Some(f.name.clone())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

    Ok(bam_files)
}
