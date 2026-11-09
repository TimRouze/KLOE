//! Shared archive metadata and output protection.
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read};
use std::path::{Path, PathBuf};

pub(crate) fn output_root(path: &str) -> &Path {
    Path::new(if path.is_empty() { "." } else { path })
}

/// Hold through input preparation, archive writing and failure-report publication.
/// This also prevents a second writer from passing the initial existence check
/// while the first writer is still downloading or building the graph.
pub(crate) struct OutputLock {
    path: PathBuf,
    _file: File,
}

impl OutputLock {
    pub(crate) fn acquire(root: &Path) -> io::Result<Self> {
        ensure_new_output(root)?;
        fs::create_dir_all(root)?;
        let path = root.join(".kloe-write.lock");
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|error| {
                io::Error::new(
                    error.kind(),
                    format!(
                        "cannot reserve output {}: {error}; choose a fresh output directory",
                        root.display()
                    ),
                )
            })?;
        let lock = Self { path, _file: file };
        ensure_new_output(root)?;
        Ok(lock)
    }
}

impl Drop for OutputLock {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

/// Never truncate a previous archive, including a legacy archive without a manifest.
pub(crate) fn ensure_new_output(root: &Path) -> io::Result<()> {
    for name in [
        crate::compress::ARCHIVE_MANIFEST_FILE,
        "filenames_id.txt",
        "tigs_kloe.fa",
        "positions_kloe.bin",
        "bucket_sizes.txt",
        crate::compress::DATASET_TO_CID_FILE,
        crate::compress::CID_TO_DATASET_FILE,
        "id_to_color_id.txt.zst",
    ] {
        match fs::symlink_metadata(root.join(name)) {
            Ok(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    format!(
                        "output {} already contains archive files; choose a fresh directory",
                        root.display()
                    ),
                ))
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error),
        }
    }
    Ok(())
}

/// Legacy archives do not record k. Invalid modern metadata must not silently fall back.
pub(crate) fn recorded_k(root: &Path) -> io::Result<Option<usize>> {
    let mut file = match File::open(root.join(crate::compress::ARCHIVE_MANIFEST_FILE)) {
        Ok(file) => file,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error),
    };
    let mut header = [0u8; 44];
    file.read_exact(&mut header)?;
    if &header[..4] != b"KLM4" || u32::from_le_bytes(header[4..8].try_into().unwrap()) != 4 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported archive manifest",
        ));
    }
    let k = u32::from_le_bytes(header[8..12].try_into().unwrap()) as usize;
    if k == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "archive k must be positive",
        ));
    }
    Ok(Some(k))
}
