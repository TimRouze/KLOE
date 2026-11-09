use std::{
    hint::unreachable_unchecked,
    sync::atomic::{AtomicUsize, Ordering},
};

use parking_lot_core::{DEFAULT_PARK_TOKEN, DEFAULT_UNPARK_TOKEN};

pub struct Notifier {
    status: AtomicUsize,
}

impl Notifier {
    const EMPTY: usize = 0;
    const WAITING: usize = 1;

    pub const fn new() -> Self {
        Self {
            status: AtomicUsize::new(Self::EMPTY),
        }
    }

    pub fn notify_one(&self) {
        if self.status.load(Ordering::SeqCst) == Self::WAITING {
            let key = &self.status as *const AtomicUsize as usize;
            unsafe {
                parking_lot_core::unpark_one(key, |status| {
                    if status.have_more_threads {
                        self.status.store(Self::WAITING, Ordering::SeqCst);
                    } else {
                        self.status.store(Self::EMPTY, Ordering::SeqCst);
                    }
                    DEFAULT_UNPARK_TOKEN
                });
            }
        }
    }

    pub fn notify_all(&self) {
        let key = &self.status as *const AtomicUsize as usize;
        unsafe {
            parking_lot_core::unpark_all(key, DEFAULT_UNPARK_TOKEN);
        }
    }

    #[inline(always)]
    pub fn wait_for_condition(&self, mut checker: impl FnMut() -> bool) {
        while !checker() {
            if self.wait_for_condition_slow(&mut checker) {
                return;
            }
        }
    }

    #[cold]
    #[inline(never)]
    fn wait_for_condition_slow(&self, checker: &mut impl FnMut() -> bool) -> bool {
        let key = &self.status as *const AtomicUsize as usize;

        self.status.store(Self::WAITING, Ordering::SeqCst);

        let mut condition_ok = false;
        match unsafe {
            parking_lot_core::park(
                key,
                || {
                    condition_ok = checker();
                    !condition_ok && self.status.load(Ordering::SeqCst) == Self::WAITING
                },
                || {},
                |_, _| {},
                DEFAULT_PARK_TOKEN,
                None,
            )
        } {
            parking_lot_core::ParkResult::Unparked(_) => {
                // Try to read again
                false
            }
            parking_lot_core::ParkResult::Invalid => {
                // Condition can be satisfied
                condition_ok
            }
            parking_lot_core::ParkResult::TimedOut => unsafe { unreachable_unchecked() },
        }
    }
}
