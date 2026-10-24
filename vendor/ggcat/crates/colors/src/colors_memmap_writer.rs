// use crate::storage::roaring::ColorsStorage;
use crate::color_lookup::ColorLookup;
use crate::storage::ColorsSerializerTrait;
use crate::storage::serializer::ColorsSerializer;
use config::ColorIndexType;
use rand::{RngCore, thread_rng};
use siphasher::sip128::{Hasher128, SipHasher13};
use std::hash::Hash;
use std::mem::ManuallyDrop;
use std::path::Path;

pub struct ColorsMemMapWriter<C: ColorsSerializerTrait> {
    colors: ColorLookup,
    colors_storage: ManuallyDrop<ColorsSerializer<C>>,
    hash_keys: (u64, u64),
    print_stats: bool,
}

impl<C: ColorsSerializerTrait> ColorsMemMapWriter<C> {
    pub fn new(
        file: impl AsRef<Path>,
        color_names: &[String],
        implicit_colors_count: Option<u64>,
        threads_count: usize,
        print_stats: bool,
    ) -> anyhow::Result<Self> {
        let mut rng = thread_rng();
        Ok(Self {
            colors: ColorLookup::new(rayon::current_num_threads()),
            colors_storage: ManuallyDrop::new(ColorsSerializer::new(
                file,
                color_names,
                implicit_colors_count,
                threads_count,
                print_stats,
            )?),
            hash_keys: (rng.next_u64(), rng.next_u64()),
            print_stats,
        })
    }

    fn hash_colors(&self, colors: &[ColorIndexType]) -> u128 {
        let mut hasher = SipHasher13::new_with_keys(self.hash_keys.0, self.hash_keys.1);
        colors.hash(&mut hasher);
        hasher.finish128().as_u128()
    }

    #[inline(always)]
    pub fn get_id(&self, colors: &[ColorIndexType]) -> ColorIndexType {
        let hash = self.hash_colors(colors);

        self.colors
            .get_or_insert(hash, || self.colors_storage.serialize_colors(colors))
    }
}

impl<C: ColorsSerializerTrait> Drop for ColorsMemMapWriter<C> {
    fn drop(&mut self) {
        if self.print_stats {
            let (entries, slots) = self.colors.allocation_stats();
            ggcat_logging::info!(
                "Color lookup: fingerprint_bits=96 entry_bytes=16 entries={} slots={} entry_storage_bytes={} previous_layout_entry_storage_bytes={} (excluding control bytes)",
                entries,
                slots,
                slots * 16,
                slots * 32
            );
        }
        unsafe { ManuallyDrop::take(&mut self.colors_storage).finalize() };
    }
}
