use std::sync::atomic::{AtomicU64, Ordering};

static MAX_DISK_USAGE: AtomicU64 = AtomicU64::new(0);
static MAX_FILES_USAGE: AtomicU64 = AtomicU64::new(0);
static CURRENT_DISK_USAGE: AtomicU64 = AtomicU64::new(0);
static CURRENT_FILES_USAGE: AtomicU64 = AtomicU64::new(0);

pub struct MemoryFsStats {
    pub max_disk_usage: u64,
    pub max_files_usage: u64,
    pub current_disk_usage: u64,
    pub current_files_usage: u64,
}

pub fn reset() {
    MAX_DISK_USAGE.store(
        CURRENT_DISK_USAGE.load(Ordering::Relaxed),
        Ordering::Relaxed,
    );
    MAX_FILES_USAGE.store(
        CURRENT_FILES_USAGE.load(Ordering::Relaxed),
        Ordering::Relaxed,
    );
}

pub fn get_stats() -> MemoryFsStats {
    MemoryFsStats {
        max_disk_usage: MAX_DISK_USAGE.load(Ordering::Relaxed),
        max_files_usage: MAX_FILES_USAGE.load(Ordering::Relaxed),
        current_disk_usage: CURRENT_DISK_USAGE.load(Ordering::Relaxed),
        current_files_usage: CURRENT_FILES_USAGE.load(Ordering::Relaxed),
    }
}

pub fn add_files_usage(usage: u64) {
    let current = CURRENT_FILES_USAGE.fetch_add(usage, Ordering::Relaxed) + usage;
    MAX_FILES_USAGE.fetch_max(current, Ordering::Relaxed);
}

pub fn decrease_files_usage(usage: u64) {
    CURRENT_FILES_USAGE.fetch_sub(usage, Ordering::Relaxed);
}

pub fn add_disk_usage(usage: u64) {
    let current = CURRENT_DISK_USAGE.fetch_add(usage, Ordering::Relaxed) + usage;
    MAX_DISK_USAGE.fetch_max(current, Ordering::Relaxed);
}

pub fn decrease_disk_usage(usage: u64) {
    CURRENT_DISK_USAGE.fetch_sub(usage, Ordering::Relaxed);
}
