use anyhow::Result;

use chrono::{DateTime, Utc};
use cloud_storage::{object::ObjectList, sync::*, ListRequest};

#[must_use]
pub fn gcs_split_path(path: &str) -> (String, String) {
    let re = regex::Regex::new(r"^gs://").unwrap();
    let path = re.replace(path, "");
    let split: Vec<&str> = path.split('/').collect();

    let bucket_name = split[0].to_string();
    let prefix = split[1..].join("/");

    (bucket_name, prefix)
}

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

pub fn gcs_get_file_update_time(path: &str) -> Result<DateTime<Utc>> {
    let (bucket_name, prefix) = gcs_split_path(path);

    let client = Client::new()?;
    let object = client.object().read(&bucket_name, &prefix)?;

    Ok(object.updated)
}

pub fn gcs_download_file(path: String) -> Result<String> {
    let (bucket_name, prefix) = gcs_split_path(&path);
    let filename = prefix.split('/').last().unwrap_or_default().to_string();

    if !std::path::Path::new(&filename).exists() {
        let client = Client::new().unwrap();
        let bytes = client.object().download(&bucket_name, &prefix).unwrap();

        std::fs::write(&filename, bytes)?;
    }

    Ok(filename)
}

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
