use std::{
    cell::UnsafeCell,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use crate::execution_manager::notifier::Notifier;

pub struct Scheduler {
    max_threads: usize,
    current_running: AtomicUsize,
    notifier: Notifier,
}

impl Scheduler {
    pub fn new(max_threads: usize) -> Arc<Self> {
        Arc::new(Self {
            max_threads,
            current_running: AtomicUsize::new(0),
            notifier: Notifier::new(),
        })
    }

    pub fn max_threads(&self) -> usize {
        self.max_threads
    }

    pub fn current_running(&self) -> usize {
        self.current_running.load(Ordering::Relaxed)
    }
}

thread_local! {
    static SCHEDULER_GUARD: UnsafeCell<Option<SchedulerGuard>> = const { UnsafeCell::new(None) };
}

struct SchedulerGuard {
    scheduler: Arc<Scheduler>,
}

impl Drop for SchedulerGuard {
    fn drop(&mut self) {
        self.scheduler
            .current_running
            .fetch_sub(1, Ordering::Relaxed);
        self.scheduler.notifier.notify_one();
    }
}

pub fn get_current_scheduler() -> Arc<Scheduler> {
    SCHEDULER_GUARD.with(|guard| {
        let guard = unsafe { &*guard.get() };
        guard.as_ref().map_or_else(
            || panic!("Scheduler is not initialized for the current thread"),
            |g| g.scheduler.clone(),
        )
    })
}

pub fn init_current_thread(scheduler: Arc<Scheduler>) {
    SCHEDULER_GUARD.with(|guard| {
        let guard = unsafe { &mut *guard.get() };
        // assert!(
        //     guard.is_none(),
        //     "Scheduler is already initialized for the current thread"
        // );
        *guard = Some(SchedulerGuard {
            scheduler: scheduler.clone(),
        });
        scheduler.current_running.fetch_add(1, Ordering::Relaxed);
    });
}

pub fn uninit_current_thread() {
    SCHEDULER_GUARD.with(|guard| {
        let guard = unsafe { &mut *guard.get() };
        *guard = None;
    });
}

#[inline(never)]
#[cold]
pub fn run_blocking_op<T>(f: impl FnOnce() -> T) -> T {
    let Some(guard) = SCHEDULER_GUARD.with(|guard| unsafe { &mut *guard.get() }.as_mut()) else {
        return f();
    };

    let active = guard
        .scheduler
        .current_running
        .fetch_sub(1, Ordering::Relaxed);
    guard.scheduler.notifier.notify_one();

    assert!(active > 0);

    let value = f();

    fn reserve_thread(guard: &SchedulerGuard) -> bool {
        guard
            .scheduler
            .current_running
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| {
                if v < guard.scheduler.max_threads {
                    Some(v + 1)
                } else {
                    None
                }
            })
            .is_ok()
    }

    guard
        .scheduler
        .notifier
        .wait_for_condition(|| reserve_thread(guard));
    value
}
