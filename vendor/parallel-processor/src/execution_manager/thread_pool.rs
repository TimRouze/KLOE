use crate::execution_manager::executor::{
    AddressConsumer, AddressProducer, AsyncExecutor, ExecutorReceiver,
};
use crate::execution_manager::notifier::Notifier;
use crate::execution_manager::objects_pool::ObjectsPool;
use crate::execution_manager::packet::Packet;
use crate::execution_manager::packets_channel::bounded::{
    packets_channel_bounded, PacketsChannelReceiverBounded, PacketsChannelSenderBounded,
};
use crate::execution_manager::packets_channel::unbounded::{
    packets_channel_unbounded, PacketsChannelSenderUnbounded,
};
use crate::execution_manager::scheduler::{
    init_current_thread, run_blocking_op, uninit_current_thread, Scheduler,
};
use std::marker::PhantomData;
use std::mem::transmute;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Barrier, Weak};
use std::thread::JoinHandle;

pub struct ExecThreadPool<E: AsyncExecutor> {
    name: String,
    threads_count: usize,
    threads: Vec<JoinHandle<()>>,
    scoped_thread_pool: Option<Arc<ScopedThreadPool>>,
    allow_spawning: bool,
    _phantom: PhantomData<E>,
}

pub struct ExecutorsHandle<E: AsyncExecutor> {
    pub(crate) spawner: PacketsChannelSenderUnbounded<AddressConsumer<E>>,
    pub(crate) channels_pool:
        Arc<ObjectsPool<PacketsChannelReceiverBounded<Packet<E::InputPacket>>>>,
    pub(crate) threads_count: usize,
}
impl<E: AsyncExecutor> Clone for ExecutorsHandle<E> {
    fn clone(&self) -> Self {
        ExecutorsHandle {
            spawner: self.spawner.clone(),
            channels_pool: self.channels_pool.clone(),
            threads_count: self.threads_count,
        }
    }
}

impl<E: AsyncExecutor> ExecutorsHandle<E> {
    pub fn add_input_data(
        &self,
        init_data: E::InitData,
        data: impl ExactSizeIterator<Item = E::InputPacket>,
    ) {
        let init_data = Arc::new(init_data);
        if E::ALLOW_PARALLEL_ADDRESS_EXECUTION {
            let address = self.create_new_address(init_data, false);
            for value in data {
                address.send_packet(Packet::new_simple(value));
            }
        } else {
            for value in data {
                let address = self.create_new_address(init_data.clone(), false);
                address.send_packet(Packet::new_simple(value));
            }
        }
    }

    pub fn create_new_address(
        &self,
        data: Arc<E::InitData>,
        high_priority: bool,
    ) -> AddressProducer<E::InputPacket> {
        let channel = self.channels_pool.alloc_object();
        let sender = channel.make_sender();
        self.spawner.send_with_priority(
            AddressConsumer {
                init_data: data,
                packets_queue: Arc::new(channel),
            },
            high_priority,
            None,
        );

        AddressProducer {
            packets_queue: sender,
        }
    }

    pub fn create_new_address_with_limit(
        &self,
        data: Arc<E::InitData>,
        high_priority: bool,
        max_in_queue: usize,
    ) -> AddressProducer<E::InputPacket> {
        let channel = self.channels_pool.alloc_object();
        let sender = channel.make_sender();
        self.spawner.send_with_priority(
            AddressConsumer {
                init_data: data,
                packets_queue: Arc::new(channel),
            },
            high_priority,
            Some(max_in_queue),
        );

        AddressProducer {
            packets_queue: sender,
        }
    }

    pub fn create_new_addresses(
        &self,
        addr_iterator: impl Iterator<Item = Arc<E::InitData>>,
        out_addresses: &mut Vec<AddressProducer<E::InputPacket>>,
        max_in_queue: Option<usize>,
        high_priority: bool,
    ) {
        out_addresses.reserve(addr_iterator.size_hint().0);

        let mut count = 0;
        self.spawner.send_batch(
            addr_iterator.map(|init_data| {
                let channel = self.channels_pool.alloc_object();
                let sender = channel.make_sender();
                count += 1;

                out_addresses.push(AddressProducer {
                    packets_queue: sender,
                });
                AddressConsumer {
                    init_data,
                    packets_queue: Arc::new(channel),
                }
            }),
            max_in_queue,
            high_priority,
        );

        assert!(
            count <= self.threads_count,
            "Cannot create more parallel addresses than the number of threads {} vs {}",
            count,
            self.threads_count
        );
    }

    #[inline(always)]
    pub fn get_pending_executors_count(&self) -> usize {
        self.spawner.len()
    }
}

impl<E: AsyncExecutor> ExecThreadPool<E> {
    pub fn new(threads_count: usize, name: &str, allow_spawning: bool) -> Self {
        Self {
            name: name.to_string(),
            threads_count,
            threads: Vec::with_capacity(threads_count),
            scoped_thread_pool: None,
            allow_spawning,
            _phantom: PhantomData,
        }
    }

    pub fn start(
        &mut self,
        scheduler: Arc<Scheduler>,
        global_params: &Arc<E::GlobalParams>,
    ) -> ExecutorsHandle<E> {
        let (addresses_sender, addresses_receiver) = packets_channel_unbounded();

        let name = self.name.clone();
        let scoped_thread_pool = if self.allow_spawning {
            Some(Arc::new(ScopedThreadPool::new(
                self.threads_count,
                &name,
                &scheduler,
            )))
        } else {
            None
        };
        self.scoped_thread_pool = scoped_thread_pool.clone();

        let handle = ExecutorsHandle {
            spawner: addresses_sender,
            // The channel has unlimited capacity, the actual capacity is limited by the thread count and the pool capacity
            channels_pool: Arc::new(ObjectsPool::new(
                self.threads_count,
                if E::ALLOW_PARALLEL_ADDRESS_EXECUTION {
                    self.threads_count * 2
                } else {
                    2
                },
            )),
            threads_count: self.threads_count,
        };

        let barrier = Arc::new(Barrier::new(self.threads_count));

        for i in 0..self.threads_count {
            let global_params = global_params.clone();
            let addresses_receiver = addresses_receiver.clone();
            let scheduler = scheduler.clone();
            let scoped_thread_pool = scoped_thread_pool.clone();
            let barrier = barrier.clone();

            let thread = std::thread::Builder::new()
                .name(format!("{}-{}", self.name, i))
                .spawn(move || {
                    init_current_thread(scheduler.clone());

                    let mut executor = E::new();
                    executor.executor_main(
                        &global_params,
                        ExecutorReceiver {
                            addresses_receiver,
                            thread_pool: scoped_thread_pool,
                            barrier,
                        },
                    );
                })
                .expect("Failed to spawn thread");
            self.threads.push(thread);
        }

        handle
    }

    pub fn join(mut self) {
        for thread in self.threads {
            if let Err(e) = thread.join() {
                eprintln!("Error joining thread: {:?}", e);
            }
        }
        if let Some(thread_pool) = self.scoped_thread_pool.take() {
            thread_pool.dispose();
        }
    }
}

pub struct ScopedThreadPool {
    threads: Vec<JoinHandle<()>>,
    execution_queue: PacketsChannelSenderBounded<Arc<ScopedTaskData>>,
}

struct ScopedTaskData {
    running_count: AtomicUsize,
    running_notifier: Notifier,
    function: Weak<dyn Fn(usize) + Sync + Send>,
    running_index: AtomicUsize,
}

impl ScopedThreadPool {
    pub fn new(threads_count: usize, name: &str, scheduler: &Arc<Scheduler>) -> Self {
        let (execution_queue, execution_dispatch) =
            packets_channel_bounded::<Arc<ScopedTaskData>>(threads_count * 2);

        Self {
            threads: (0..threads_count)
                .map(|i| {
                    let execution_dispatch = execution_dispatch.clone();
                    let scheduler = scheduler.clone();
                    std::thread::Builder::new()
                        .name(format!("{}-{}-threadpool", name, i))
                        .spawn(move || {
                            while let Some(task) = execution_dispatch.recv() {
                                task.running_count.fetch_add(1, Ordering::SeqCst);
                                if let Some(function) = task.function.upgrade() {
                                    let running_index =
                                        task.running_index.fetch_add(1, Ordering::Relaxed);
                                    init_current_thread(scheduler.clone());
                                    function(running_index);
                                    uninit_current_thread();
                                }
                                let count = task.running_count.fetch_sub(1, Ordering::SeqCst);
                                if count == 1 {
                                    task.running_notifier.notify_all();
                                }
                            }
                        })
                        .unwrap()
                })
                .collect(),
            execution_queue,
        }
    }

    pub fn run_scoped_optional(
        &self,
        concurrency: usize,
        callback: impl Fn(usize) + Sync + Send + Copy,
    ) {
        assert!(concurrency > 0);
        if concurrency == 1 {
            callback(0);
        } else {
            let function: Arc<dyn Fn(usize) + Sync + Send> = Arc::new(callback);
            let function_lt_extended: Arc<dyn Fn(usize) + Sync + Send> =
                unsafe { transmute(function) };

            let task_data = Arc::new(ScopedTaskData {
                running_count: AtomicUsize::new(0),
                running_notifier: Notifier::new(),
                function: Arc::downgrade(&function_lt_extended),
                running_index: AtomicUsize::new(1),
            });

            for _ in 1..concurrency {
                self.execution_queue.send(task_data.clone());
            }

            callback(0);

            drop(function_lt_extended);
            // Wait for all the running tasks to terminate
            if task_data.running_count.load(Ordering::SeqCst) > 0 {
                run_blocking_op(|| {
                    task_data
                        .running_notifier
                        .wait_for_condition(|| task_data.running_count.load(Ordering::SeqCst) == 0);
                });
            }
        }
    }

    pub fn dispose(&self) {
        self.execution_queue.dispose();
    }
}

impl Drop for ScopedThreadPool {
    fn drop(&mut self) {
        for thread in self.threads.drain(..) {
            thread.join().unwrap();
        }
    }
}
