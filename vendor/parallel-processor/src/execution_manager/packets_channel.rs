#![allow(dead_code)]

use std::sync::{
    atomic::{AtomicBool, Ordering},
    Arc,
};

pub(crate) struct PacketsChannelReceiver<Q> {
    queue: Arc<Q>,
}
impl<Q> Clone for PacketsChannelReceiver<Q> {
    fn clone(&self) -> Self {
        Self {
            queue: self.queue.clone(),
        }
    }
}

pub(crate) struct PacketsChannelSender<Q: PacketsQueue> {
    queue: Arc<Q>,
    is_disposed: AtomicBool,
}
impl<Q: PacketsQueue> Clone for PacketsChannelSender<Q> {
    fn clone(&self) -> Self {
        self.queue.incr_senders_count();
        Self {
            queue: self.queue.clone(),
            is_disposed: AtomicBool::new(false),
        }
    }
}

impl<Q: PacketsQueue> Drop for PacketsChannelSender<Q> {
    fn drop(&mut self) {
        if !self.is_disposed.load(Ordering::Relaxed) {
            self.queue.decr_senders_count();
        }
    }
}

pub(crate) trait PacketsQueue {
    type Item;
    fn push(&self, value: Self::Item);
    fn try_pop(&self) -> Option<Self::Item>;
    fn pop(&self) -> Option<Self::Item>;
    fn len(&self) -> usize;
    fn incr_senders_count(&self);
    fn decr_senders_count(&self);
    fn get_senders_count(&self) -> usize;
}

impl<Q: PacketsQueue> PacketsChannelSender<Q> {
    #[inline(always)]
    pub fn send(&self, value: Q::Item) {
        self.queue.push(value);
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.queue.len()
    }

    pub fn dispose(&self) {
        if !self.is_disposed.swap(true, Ordering::Relaxed) {
            self.queue.decr_senders_count();
        }
    }
}

impl<Q: PacketsQueue> PacketsChannelReceiver<Q> {
    #[inline(always)]
    pub fn try_recv(&self) -> Option<Q::Item> {
        self.queue.try_pop()
    }

    #[inline(always)]
    pub fn recv(&self) -> Option<Q::Item> {
        self.queue.pop()
    }

    #[inline(always)]
    pub fn is_active(&self) -> bool {
        self.queue.get_senders_count() > 0 || self.queue.len() > 0
    }

    #[inline(always)]
    pub fn make_sender(&self) -> PacketsChannelSender<Q> {
        self.queue.incr_senders_count();
        PacketsChannelSender {
            queue: self.queue.clone(),
            is_disposed: AtomicBool::new(false),
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.queue.len()
    }
}

pub mod bounded {
    use std::{
        mem::{forget, ManuallyDrop},
        sync::{
            atomic::{AtomicBool, AtomicUsize, Ordering},
            Arc,
        },
    };

    use crossbeam::queue::ArrayQueue;

    use crate::execution_manager::{
        notifier::Notifier,
        objects_pool::PoolObjectTrait,
        packets_channel::{PacketsChannelReceiver, PacketsChannelSender, PacketsQueue},
        scheduler::run_blocking_op,
    };

    pub(crate) struct BoundedQueue<T> {
        queue: ArrayQueue<T>,
        senders_waiting: Notifier,
        receivers_waiting: Notifier,
        senders_count: AtomicUsize,
    }

    impl<T> PacketsQueue for BoundedQueue<T> {
        type Item = T;

        fn push(&self, value: Self::Item) {
            if let Err(value) = self.queue.push(value) {
                let mut value = ManuallyDrop::new(value);
                run_blocking_op(|| {
                    self.senders_waiting.wait_for_condition(|| {
                        match self.queue.push(unsafe { ManuallyDrop::take(&mut value) }) {
                            Ok(_) => true,
                            Err(value) => {
                                // Forget this value as it is only a copy
                                forget(value);

                                // If we failed to push, it means the queue is full
                                // and we need to wait for a receiver to pop an item
                                false
                            }
                        }
                    });
                });
            }
            self.receivers_waiting.notify_one();
        }

        fn try_pop(&self) -> Option<Self::Item> {
            let value = self.queue.pop()?;
            self.senders_waiting.notify_one();
            Some(value)
        }

        fn pop(&self) -> Option<Self::Item> {
            let value = if let Some(value) = self.queue.pop() {
                Some(value)
            } else {
                let mut value = None;

                run_blocking_op(|| {
                    self.receivers_waiting
                        .wait_for_condition(|| match self.queue.pop() {
                            Some(v) => {
                                value = Some(v);
                                true
                            }
                            None => self.senders_count.load(Ordering::Relaxed) == 0,
                        });
                });
                value
            };
            self.senders_waiting.notify_one();
            value
        }

        #[inline(always)]
        fn len(&self) -> usize {
            self.queue.len()
        }

        #[inline(always)]
        fn incr_senders_count(&self) {
            self.senders_count.fetch_add(1, Ordering::Relaxed);
        }

        #[inline(always)]
        fn decr_senders_count(&self) {
            let senders_count = self.senders_count.fetch_sub(1, Ordering::Relaxed);
            if senders_count == 1 {
                self.receivers_waiting.notify_all();
            }
        }

        #[inline(always)]
        fn get_senders_count(&self) -> usize {
            self.senders_count.load(Ordering::Relaxed)
        }
    }

    pub(crate) type PacketsChannelReceiverBounded<T> = PacketsChannelReceiver<BoundedQueue<T>>;
    pub(crate) type PacketsChannelSenderBounded<T> = PacketsChannelSender<BoundedQueue<T>>;

    impl<T> PacketsChannelReceiverBounded<T> {}
    impl<T> PacketsChannelSenderBounded<T> {}

    impl<T: Send + 'static> PoolObjectTrait for PacketsChannelReceiverBounded<T> {
        type InitData = usize; // max_size

        fn allocate_new(init_data: &Self::InitData) -> Self {
            Self {
                queue: Arc::new(BoundedQueue {
                    queue: ArrayQueue::new(*init_data),
                    senders_waiting: Notifier::new(),
                    receivers_waiting: Notifier::new(),
                    senders_count: AtomicUsize::new(0),
                }),
            }
        }

        fn reset(&mut self) {
            assert_eq!(
                self.queue.senders_count.load(Ordering::Relaxed),
                0,
                "Cannot reset PacketsChannelReceiver while senders are active"
            );
            assert_eq!(
                self.queue.len(),
                0,
                "Cannot reset PacketsChannelReceiver while there are items in the queue"
            );
        }
    }

    pub(crate) fn packets_channel_bounded<T: Send + 'static>(
        max_size: usize,
    ) -> (
        PacketsChannelSenderBounded<T>,
        PacketsChannelReceiverBounded<T>,
    ) {
        let internal = Arc::new(BoundedQueue {
            queue: ArrayQueue::new(max_size),
            receivers_waiting: Notifier::new(),
            senders_waiting: Notifier::new(),
            senders_count: AtomicUsize::new(1),
        });
        (
            PacketsChannelSender {
                queue: internal.clone(),
                is_disposed: AtomicBool::new(false),
            },
            PacketsChannelReceiver { queue: internal },
        )
    }
}

pub mod unbounded {
    use crate::execution_manager::{
        packets_channel::{PacketsChannelReceiver, PacketsChannelSender, PacketsQueue},
        scheduler::run_blocking_op,
    };
    use parking_lot::{Condvar, Mutex};
    use std::{
        collections::VecDeque,
        sync::{
            atomic::{AtomicBool, AtomicUsize},
            Arc,
        },
        time::Duration,
    };

    pub(crate) struct UnboundedQueue<T> {
        queue: Mutex<VecDeque<T>>,
        receivers_waiting: Condvar,
        senders_waiting: Condvar,
        senders_count: AtomicUsize,
    }

    impl<T> PacketsQueue for UnboundedQueue<T> {
        type Item = T;

        fn push(&self, value: Self::Item) {
            let mut queue = self.queue.lock();
            queue.push_back(value);
            drop(queue);
            self.receivers_waiting.notify_one();
        }

        fn try_pop(&self) -> Option<Self::Item> {
            self.queue.lock().pop_front()
        }

        fn pop(&self) -> Option<Self::Item> {
            let mut queue = self.queue.lock();
            let value = if let Some(value) = queue.pop_front() {
                drop(queue);
                Some(value)
            } else {
                run_blocking_op(|| loop {
                    match queue.pop_front() {
                        Some(v) => {
                            drop(queue);
                            return Some(v);
                        }
                        None => {
                            if self
                                .senders_count
                                .load(std::sync::atomic::Ordering::Relaxed)
                                == 0
                            {
                                drop(queue);
                                return None;
                            } else {
                                self.receivers_waiting.wait(&mut queue);
                            }
                        }
                    }
                })
            };
            self.senders_waiting.notify_one();
            value
        }

        #[inline(always)]
        fn len(&self) -> usize {
            self.queue.lock().len()
        }

        #[inline(always)]
        fn incr_senders_count(&self) {
            self.senders_count
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }

        #[inline(always)]
        fn decr_senders_count(&self) {
            let senders_count = self
                .senders_count
                .fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
            if senders_count == 1 {
                let _lock = self.queue.lock();
                self.receivers_waiting.notify_all();
            }
        }

        #[inline(always)]
        fn get_senders_count(&self) -> usize {
            self.senders_count
                .load(std::sync::atomic::Ordering::Relaxed)
        }
    }

    impl<T> PacketsChannelSenderUnbounded<T> {
        fn wait_for_space(&self, max_in_queue: usize) {
            run_blocking_op(|| {
                let mut queue = self.queue.queue.lock();
                while queue.len() > max_in_queue {
                    self.queue
                        .senders_waiting
                        .wait_for(&mut queue, Duration::from_millis(50));
                }
            });
        }

        #[inline(always)]
        pub fn send_batch(
            &self,
            values: impl Iterator<Item = T>,
            max_in_queue: Option<usize>,
            high_priority: bool,
        ) {
            if let Some(max_in_queue) = max_in_queue {
                self.wait_for_space(max_in_queue);
            }

            let mut queue = self.queue.queue.lock();
            if high_priority {
                for value in values {
                    queue.push_front(value);
                }
            } else {
                for value in values {
                    queue.push_back(value);
                }
            }
            drop(queue);
            self.queue.receivers_waiting.notify_all();
        }

        pub fn send_with_priority(
            &self,
            value: T,
            high_priority: bool,
            max_in_queue: Option<usize>,
        ) {
            if let Some(max_in_queue) = max_in_queue {
                self.wait_for_space(max_in_queue);
            }

            let mut queue = self.queue.queue.lock();
            if high_priority {
                queue.push_front(value);
            } else {
                queue.push_back(value);
            }
            drop(queue);
            self.queue.receivers_waiting.notify_one();
        }
    }

    pub(crate) type PacketsChannelReceiverUnbounded<T> = PacketsChannelReceiver<UnboundedQueue<T>>;
    pub(crate) type PacketsChannelSenderUnbounded<T> = PacketsChannelSender<UnboundedQueue<T>>;

    pub(crate) fn packets_channel_unbounded<T: Send + 'static>() -> (
        PacketsChannelSenderUnbounded<T>,
        PacketsChannelReceiverUnbounded<T>,
    ) {
        let internal = Arc::new(UnboundedQueue {
            queue: Mutex::new(VecDeque::with_capacity(64)),
            receivers_waiting: Condvar::new(),
            senders_waiting: Condvar::new(),
            senders_count: AtomicUsize::new(1),
        });
        (
            PacketsChannelSender {
                queue: internal.clone(),
                is_disposed: AtomicBool::new(false),
            },
            PacketsChannelReceiver { queue: internal },
        )
    }
}
