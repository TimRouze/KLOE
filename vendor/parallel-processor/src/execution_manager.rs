pub mod executor;
pub mod notifier;
pub mod objects_pool;
pub mod packet;
pub mod packets_channel;
pub mod scheduler;
pub mod thread_pool;

#[cfg(test)]
mod tests {

    use crate::execution_manager::executor::{AsyncExecutor, ExecutorReceiver};
    use crate::execution_manager::objects_pool::PoolObjectTrait;
    use crate::execution_manager::packet::{PacketTrait, PacketsPool};
    use crate::execution_manager::scheduler::{run_blocking_op, Scheduler};
    use crate::execution_manager::thread_pool::{ExecThreadPool, ExecutorsHandle};
    use crate::set_logger_function;
    use std::num::Wrapping;
    use std::sync::atomic::AtomicUsize;
    use std::sync::Arc;

    struct TestExecutor {}

    impl PoolObjectTrait for usize {
        type InitData = ();

        fn allocate_new(_init_data: &Self::InitData) -> Self {
            0
        }

        fn reset(&mut self) {}
    }

    impl PacketTrait for usize {
        fn get_size(&self) -> usize {
            0
        }
    }

    impl AsyncExecutor for TestExecutor {
        type InputPacket = usize;
        type OutputPacket = usize;
        type GlobalParams = AtomicUsize;
        type InitData = ExecutorsHandle<Self>;
        const ALLOW_PARALLEL_ADDRESS_EXECUTION: bool = false;

        fn new() -> Self {
            Self {}
        }

        fn executor_main<'a>(
            &'a mut self,
            _global_params: &'a Self::GlobalParams,
            mut receiver: ExecutorReceiver<Self>,
            // _memory_tracker: MemoryTracker<Self>,
        ) {
            static INDEX: AtomicUsize = AtomicUsize::new(0);
            let index = INDEX.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let pool = PacketsPool::<usize>::new(1024, ());
            while let Ok(addr) = receiver.obtain_address() {
                // let pool = addr.get_pool();
                while let Some(packet) = addr.receive_packet() {
                    // std::thread::sleep(Duration::from_millis(10));
                    // crate::log_info!("v: {}", x);
                    crate::log_info!("X: {}", *packet);

                    let pvalue = *packet;

                    addr.spawn_executors(3, |_| {
                        run_blocking_op(|| {});
                        // println!("ENABLED count: {}", scheduler.current_running());

                        let address = addr
                            .get_init_data()
                            .create_new_address(Arc::new(addr.get_init_data().clone()), false);

                        let mut x = Wrapping(pvalue);
                        for i in 1..100000 {
                            let i = Wrapping(i as usize);
                            x += i * i + x;
                        }
                        // println!("FINISHED count: {}", scheduler.current_running());

                        if x.0 == usize::MAX {
                            crate::log_info!("Overflow detected in executor {}", index);
                        }

                        let mut npacket = pool.alloc_packet();
                        *npacket = 0 + pvalue * 2;

                        if *npacket > 64 {
                            return;
                        }

                        address.send_packet(npacket);
                    });

                    drop(packet);
                }
            }
            println!("Finished {}!", index);
        }
    }

    #[test]
    #[ignore]
    fn test_executors() {
        set_logger_function(|_level, message| {
            println!("{}", message);
        });

        let scheduler = Scheduler::new(16);

        let mut readers_pool = ExecThreadPool::<TestExecutor>::new(16, "readers-pool", true);

        let running_count = Arc::new(AtomicUsize::new(0));

        let handle = readers_pool.start(scheduler, &running_count);

        let strings = vec![1]; //, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        handle.add_input_data(handle.clone(), strings.into_iter());

        // test_input.set_output_executor::<TestExecutor>(&context, (), 0);

        drop(handle);

        // loop {
        //     // println!(
        //     //     "Running count: {}",
        //     //     running_count.load(std::sync::atomic::Ordering::Relaxed)
        //     // );
        //     std::thread::sleep(Duration::from_millis(1000));
        // }

        // readers_pool
        readers_pool.join();
    }
}
