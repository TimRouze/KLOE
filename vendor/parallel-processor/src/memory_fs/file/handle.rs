use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{Seek, SeekFrom, Write},
    path::PathBuf,
    sync::{Arc, Weak},
};

use once_cell::sync::Lazy;
use parking_lot::{ArcMutexGuard, Mutex};

static OPENED_FILES: Lazy<Mutex<HashMap<PathBuf, Arc<Mutex<File>>>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));
static MAX_OPENED_FILES: Lazy<usize> = Lazy::new(|| {
    const DEFAULT_MAX_OPENED_FILES: u32 = 1024;

    let os_limit = match () {
        #[cfg(target_os = "linux")]
        () => limits_rs::get_own_limits()
            .map(|l| {
                l.max_open_files
                    .soft
                    .unwrap_or(l.max_open_files.hard.unwrap_or(DEFAULT_MAX_OPENED_FILES))
            })
            .unwrap_or(DEFAULT_MAX_OPENED_FILES) as usize,
        #[cfg(not(target_os = "linux"))]
        () => DEFAULT_MAX_OPENED_FILES as usize,
    };
    os_limit.saturating_sub(30).min(os_limit / 2).max(4)
});

pub struct FileHandle {
    path: PathBuf,
    file: Weak<Mutex<File>>,
}

impl FileHandle {
    fn open_file(path: &PathBuf, create_new: bool) -> Arc<Mutex<File>> {
        let mut opened_files = OPENED_FILES.lock();

        if opened_files.len() >= *MAX_OPENED_FILES {
            // Check which file to close
            for open_file in opened_files.iter() {
                if Arc::strong_count(open_file.1) == 1 {
                    let to_remove = open_file.0.clone();
                    opened_files.remove(&to_remove);
                    // println!("Closing file {}", to_remove.display());
                    break;
                }
            }
        }

        let file = Arc::new(Mutex::new({
            let mut file = OpenOptions::new()
                .create(create_new)
                .write(true)
                .append(false)
                // .custom_flags(O_DIRECT)
                .open(&path)
                .map_err(|e| format!("Error while opening file {}: {}", path.display(), e))
                .unwrap();
            if !create_new {
                // Seek to the end of the file if the file is reopened
                file.seek(SeekFrom::End(0)).unwrap();
            }
            file
        }));

        opened_files.insert(path.clone(), file.clone());
        file
    }

    pub fn new(path: PathBuf) -> Self {
        let file = Self::open_file(&path, true);
        Self {
            path,
            file: Arc::downgrade(&file),
        }
    }

    pub fn flush(&self) -> std::io::Result<()> {
        // Do not flush if the file is closed
        let Some(file) = self.file.upgrade() else {
            return Ok(());
        };
        let mut file = file.lock();
        file.flush()
    }

    pub fn get_path(&self) -> &PathBuf {
        &self.path
    }

    pub fn get_file(&mut self) -> ArcMutexGuard<parking_lot::RawMutex, File> {
        if let Some(file) = self.file.upgrade() {
            file.lock_arc()
        } else {
            let file = Self::open_file(&self.path, false);
            self.file = Arc::downgrade(&file);
            file.lock_arc()
        }
    }
}

impl Drop for FileHandle {
    fn drop(&mut self) {
        // Close the file writer once the handle is dropped
        OPENED_FILES.lock().remove(&self.path);
    }
}
