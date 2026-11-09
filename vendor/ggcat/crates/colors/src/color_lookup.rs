use crossbeam::utils::CachePadded;
use hashbrown::HashTable;
use parking_lot::RwLock;

// Keep the fingerprint and CID in one table element. A (u128, u32) map entry
// occupies 32 bytes due to alignment; this layout occupies exactly 16.
#[repr(C)]
struct ColorEntry {
    hash: u64,
    fingerprint: u32,
    cid: u32,
}

const _: () = assert!(std::mem::size_of::<ColorEntry>() == 16);

pub(crate) struct ColorLookup {
    shards: Box<[CachePadded<RwLock<HashTable<ColorEntry>>>]>,
}

impl ColorLookup {
    pub(crate) fn allocation_stats(&self) -> (usize, usize) {
        self.shards.iter().fold((0, 0), |(entries, slots), shard| {
            let table = shard.read();
            let buckets = if table.capacity() == 0 {
                0
            } else {
                table.capacity().next_power_of_two()
            };
            (entries + table.len(), slots + buckets)
        })
    }

    pub(crate) fn new(threads: usize) -> Self {
        let count = (threads.max(1) * 8).next_power_of_two();
        Self {
            shards: (0..count)
                .map(|_| CachePadded::new(RwLock::new(HashTable::new())))
                .collect(),
        }
    }

    pub(crate) fn get_or_insert(&self, hash128: u128, create: impl FnOnce() -> u32) -> u32 {
        let hash = hash128 as u64;
        let fingerprint = (hash128 >> 64) as u32;
        // Use different bits for sharding and the table's low-bit bucket index.
        let shard = &self.shards[((hash >> 32) as usize) & (self.shards.len() - 1)];
        let matches = |entry: &ColorEntry| entry.hash == hash && entry.fingerprint == fingerprint;
        if let Some(entry) = shard.read().find(hash, matches) {
            return entry.cid;
        }
        let mut table = shard.write();
        // Another thread can insert between releasing the read lock and taking
        // the write lock. Allocate/serialize a CID only after this second check.
        if let Some(entry) = table.find(hash, matches) {
            return entry.cid;
        }
        let cid = create();
        // Equality uses all 96 fingerprint bits, so equal 64-bit prefixes can
        // coexist. Rehashing uses the same 64-bit hash as lookup.
        table.insert_unique(
            hash,
            ColorEntry {
                hash,
                fingerprint,
                cid,
            },
            |entry| entry.hash,
        );
        cid
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU32, Ordering};

    #[test]
    fn packed_entry_halves_previous_entry_size() {
        assert_eq!(std::mem::size_of::<ColorEntry>(), 16);
        assert_eq!(std::mem::size_of::<(u128, u32)>(), 32);
    }

    #[test]
    fn equal_hash_prefixes_keep_distinct_fingerprints_and_full_cids() {
        let lookup = ColorLookup::new(1);
        let prefix = 0x1234_5678_90ab_cdefu128;
        let a = prefix | (1u128 << 64);
        let b = prefix | (2u128 << 64);
        assert_eq!(lookup.get_or_insert(a, || u32::MAX - 1), u32::MAX - 1);
        assert_eq!(lookup.get_or_insert(b, || 42), 42);
        assert_eq!(
            lookup.get_or_insert(a, || panic!("duplicate insertion")),
            u32::MAX - 1
        );
        assert_eq!(
            lookup.get_or_insert(b, || panic!("duplicate insertion")),
            42
        );
        // Only the explicitly chosen 96 bits are retained.
        assert_eq!(
            lookup.get_or_insert(a | (7u128 << 96), || panic!("unexpected key")),
            u32::MAX - 1
        );
    }

    #[test]
    fn concurrent_insertions_and_resizing_preserve_ids() {
        let lookup = ColorLookup::new(2);
        let next = AtomicU32::new(0);
        std::thread::scope(|scope| {
            for _ in 0..8 {
                let lookup = &lookup;
                let next = &next;
                scope.spawn(move || {
                    for i in 0..4096u128 {
                        // All keys enter the same shard, with many identical
                        // 64-bit prefixes; exercise equality during resizing.
                        let key = (i << 64) | (i % 32);
                        let cid =
                            lookup.get_or_insert(key, || next.fetch_add(1, Ordering::Relaxed));
                        assert_eq!(lookup.get_or_insert(key, || panic!("lost key")), cid);
                    }
                });
            }
        });
        assert_eq!(next.load(Ordering::Relaxed), 4096);
        let mut ids: Vec<_> = (0..4096u128)
            .map(|i| lookup.get_or_insert((i << 64) | (i % 32), || panic!("missing key")))
            .collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), 4096);
        let (entries, slots) = lookup.allocation_stats();
        assert_eq!(entries, 4096);
        assert!(slots >= entries);
    }
}
