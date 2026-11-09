#![cfg_attr(test, feature(type_alias_impl_trait))]
#![allow(warnings)]

#[macro_use]
pub extern crate mt_debug_counters;

pub use mt_debug_counters::logging::enable_counters_logging;

#[macro_use]
pub mod memory_fs;
pub mod buckets;
pub mod debug_allocator;
pub mod execution_manager;
pub mod fast_smart_bucket_sort;
pub mod memory_data_size;
pub mod phase_times_monitor;
pub mod utils;

#[cfg(feature = "process-stats")]
pub mod simple_process_stats;

pub struct Utils {}

impl Utils {
    pub fn multiply_by(val: usize, mult: f64) -> usize {
        ((val as f64) * mult) as usize
    }
}
