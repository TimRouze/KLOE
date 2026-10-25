//! Discover the memory available to this process, including nested job cgroups.
use std::path::{Path, PathBuf};

const MIB: u64 = 1024 * 1024;

fn number(path: &Path) -> Option<u64> {
    std::fs::read_to_string(path).ok()?.trim().parse().ok()
}

fn unescape_mount(value: &str) -> PathBuf {
    PathBuf::from(
        value
            .replace("\\040", " ")
            .replace("\\011", "\t")
            .replace("\\012", "\n")
            .replace("\\134", "\\"),
    )
}

fn cgroup_directories(membership: &str, mounts: &str) -> Vec<(PathBuf, PathBuf, bool)> {
    let mut result = Vec::new();
    for line in membership.lines() {
        let fields: Vec<_> = line.splitn(3, ':').collect();
        if fields.len() != 3 {
            continue;
        }
        let v2 = fields[1].is_empty();
        if !v2 && !fields[1].split(',').any(|x| x == "memory") {
            continue;
        }
        for mount in mounts.lines() {
            let Some((before, after)) = mount.split_once(" - ") else {
                continue;
            };
            let left: Vec<_> = before.split_whitespace().collect();
            let right: Vec<_> = after.split_whitespace().collect();
            if left.len() < 5 || right.len() < 3 {
                continue;
            }
            if (v2 && right[0] != "cgroup2")
                || (!v2 && (right[0] != "cgroup" || !right[2].split(',').any(|x| x == "memory")))
            {
                continue;
            }
            let root = unescape_mount(left[3]);
            let mountpoint = unescape_mount(left[4]);
            let path = Path::new(fields[2]);
            // A cgroup namespace may hide the membership path. Its mounted
            // root still exposes a useful limit; never follow '..' outside it.
            if path
                .components()
                .any(|p| p == std::path::Component::ParentDir)
            {
                result.push((mountpoint.clone(), mountpoint, v2));
                continue;
            }
            let relative = path.strip_prefix(&root).or_else(|_| path.strip_prefix("/"));
            if let Ok(relative) = relative {
                result.push((mountpoint.join(relative), mountpoint, v2));
            }
        }
    }
    result
}

fn cgroup_available(dir: &Path, mount: &Path, v2: bool, own_rss: u64) -> Option<u64> {
    let mut available = None;
    for ancestor in dir.ancestors().take_while(|p| p.starts_with(mount)) {
        let usage = number(&ancestor.join(if v2 {
            "memory.current"
        } else {
            "memory.usage_in_bytes"
        }));
        let limits: &[&str] = if v2 {
            &["memory.max", "memory.high"]
        } else {
            &["memory.limit_in_bytes"]
        };
        for limit in limits {
            if let Some(limit) = number(&ancestor.join(limit)).filter(|n| *n < (1u64 << 60)) {
                // Include siblings, page cache and kernel charges, rather than RSS alone.
                let room = limit.saturating_sub(usage.unwrap_or(0).saturating_sub(own_rss));
                available = Some(available.map_or(room, |old: u64| old.min(room)));
            }
        }
    }
    available
}

fn field_kib(text: &str, name: &str) -> Option<u64> {
    text.lines()
        .find_map(|line| {
            let (key, rest) = line.split_once(':')?;
            (key == name).then(|| rest.split_whitespace().next()?.parse::<u64>().ok())?
        })
        .map(|n| n.saturating_mul(1024))
}

fn detect_available() -> Option<u64> {
    let mut bounds = Vec::new();
    let rss = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| field_kib(&s, "VmRSS"))
        .unwrap_or(0);
    if let Some(available) = std::fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|s| field_kib(&s, "MemAvailable"))
    {
        bounds.push(available.saturating_add(rss));
    }
    if let (Ok(membership), Ok(mounts)) = (
        std::fs::read_to_string("/proc/self/cgroup"),
        std::fs::read_to_string("/proc/self/mountinfo"),
    ) {
        for (dir, mount, v2) in cgroup_directories(&membership, &mounts) {
            if let Some(bound) = cgroup_available(&dir, &mount, v2, rss) {
                bounds.push(bound);
            }
        }
    }
    let env_number = |key| std::env::var(key).ok().and_then(|s| s.parse::<u64>().ok());
    if let Some(mb) = env_number("SLURM_MEM_PER_NODE").filter(|n| *n > 0) {
        bounds.push(mb.saturating_mul(MIB));
    }
    if let (Some(mb), Some(cpus)) = (
        env_number("SLURM_MEM_PER_CPU"),
        env_number("SLURM_CPUS_ON_NODE"),
    ) {
        if mb > 0 && cpus > 0 {
            bounds.push(mb.saturating_mul(cpus).saturating_mul(MIB));
        }
    }
    bounds.into_iter().min()
}

/// Cap a working budget to 80% of available job/host memory. This is a startup
/// budget, not a hard allocator limit for arbitrary application data structures.
pub fn effective_memory_budget(requested: usize) -> usize {
    static AVAILABLE: std::sync::OnceLock<Option<u64>> = std::sync::OnceLock::new();
    let available = *AVAILABLE.get_or_init(detect_available);
    available.map_or(requested, |n| {
        requested.min((n / 5 * 4).min(usize::MAX as u64) as usize)
    })
}

/// Leave room for a streaming reader's Zstd window, parser and graph worker
/// buffers. This bounds concurrency, not the size of the global color table.
pub fn streaming_threads(requested: usize, budget: usize) -> usize {
    requested.max(1).min((budget / (256 * 1024 * 1024)).max(1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streaming_concurrency_follows_available_budget() {
        assert_eq!(streaming_threads(20, 800 * 1024 * 1024), 3);
        assert_eq!(streaming_threads(64, 64 * 1024 * 1024 * 1024), 64);
        assert_eq!(streaming_threads(0, 0), 1);
    }

    #[test]
    fn resolves_v1_v2_and_namespaced_mounts() {
        let mounts = "1 0 0:1 / /sys/fs/cgroup rw - cgroup2 cgroup rw\n2 0 0:2 /slurm /sys/fs/cgroup/memory rw - cgroup cgroup rw,memory";
        assert_eq!(
            cgroup_directories("0::/job/step", mounts),
            vec![(
                "/sys/fs/cgroup/job/step".into(),
                "/sys/fs/cgroup".into(),
                true
            )]
        );
        assert_eq!(
            cgroup_directories("4:memory:/slurm/job", mounts),
            vec![(
                "/sys/fs/cgroup/memory/job".into(),
                "/sys/fs/cgroup/memory".into(),
                false
            )]
        );
        assert_eq!(
            cgroup_directories(
                "0::/",
                "1 0 0:1 /host/job /sys/fs/cgroup rw - cgroup2 cgroup rw"
            ),
            vec![("/sys/fs/cgroup/".into(), "/sys/fs/cgroup".into(), true)]
        );
        assert_eq!(
            cgroup_directories("0::/../../job", mounts),
            vec![("/sys/fs/cgroup".into(), "/sys/fs/cgroup".into(), true)]
        );
    }

    #[test]
    fn parent_limits_and_other_job_memory_are_accounted_for() {
        let root = std::env::temp_dir().join(format!("kloe-cgroup-test-{}", std::process::id()));
        let leaf = root.join("step/task");
        std::fs::create_dir_all(&leaf).unwrap();
        std::fs::write(leaf.join("memory.max"), "max").unwrap();
        std::fs::write(root.join("memory.max"), "1000").unwrap();
        std::fs::write(root.join("memory.current"), "400").unwrap();
        assert_eq!(cgroup_available(&leaf, &root, true, 100), Some(700));
        std::fs::write(leaf.parent().unwrap().join("memory.high"), "500").unwrap();
        std::fs::write(leaf.parent().unwrap().join("memory.current"), "250").unwrap();
        assert_eq!(cgroup_available(&leaf, &root, true, 100), Some(350));
        std::fs::write(root.join("memory.limit_in_bytes"), "800").unwrap();
        std::fs::write(root.join("memory.usage_in_bytes"), "300").unwrap();
        assert_eq!(cgroup_available(&leaf, &root, false, 100), Some(600));
        std::fs::remove_dir_all(root).unwrap();
    }
}
