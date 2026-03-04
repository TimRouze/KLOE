use crate::execution_manager::objects_pool::{PoolObject, PoolObjectTrait};
use crate::execution_manager::packet::{Packet, PacketTrait};
use crate::execution_manager::packets_channel::bounded::{
    PacketsChannelReceiverBounded, PacketsChannelSenderBounded,
};
use crate::execution_manager::packets_channel::unbounded::PacketsChannelReceiverUnbounded;
use crate::execution_manager::thread_pool::ScopedThreadPool;
use std::sync::{Arc, Barrier, BarrierWaitResult};

pub trait AsyncExecutor: Sized + Send + Sync + 'static {
    type InputPacket: PoolObjectTrait;
    type OutputPacket: PacketTrait;
    type GlobalParams: Send + Sync + 'static;
    type InitData: Send + Sync + 'static;
    const ALLOW_PARALLEL_ADDRESS_EXECUTION: bool;

    fn new() -> Self;

    fn executor_main<'a>(
        &'a mut self,
        global_params: &'a Self::GlobalParams,
        receiver: ExecutorReceiver<Self>,
    );
}

pub struct AddressConsumer<E: AsyncExecutor> {
    pub(crate) init_data: Arc<E::InitData>,
    pub(crate) packets_queue:
        Arc<PoolObject<PacketsChannelReceiverBounded<Packet<E::InputPacket>>>>,
}

impl<E: AsyncExecutor> Clone for AddressConsumer<E> {
    fn clone(&self) -> Self {
        AddressConsumer {
            init_data: self.init_data.clone(),
            packets_queue: self.packets_queue.clone(),
        }
    }
}

pub struct AddressProducer<P: PoolObjectTrait> {
    pub(crate) packets_queue: PacketsChannelSenderBounded<Packet<P>>,
}

impl<P: PoolObjectTrait> AddressProducer<P> {
    pub fn send_packet(&self, packet: Packet<P>) {
        self.packets_queue.send(packet);
    }
}

pub struct ExecutorReceiver<E: AsyncExecutor> {
    pub(crate) addresses_receiver: PacketsChannelReceiverUnbounded<AddressConsumer<E>>,
    pub(crate) thread_pool: Option<Arc<ScopedThreadPool>>,
    pub(crate) barrier: Arc<Barrier>,
}

impl<E: AsyncExecutor> ExecutorReceiver<E> {
    pub fn obtain_address(&mut self) -> Result<ExecutorAddressOperations<E>, ()> {
        let address_receiver = match self.addresses_receiver.recv() {
            Some(value) => value,
            None => return Err(()),
        };

        if E::ALLOW_PARALLEL_ADDRESS_EXECUTION && address_receiver.packets_queue.is_active() {
            // Reinsert the current address if it is allowed to be run on multiple threads
            self.addresses_receiver
                .make_sender()
                .send(address_receiver.clone());
        }

        Ok(ExecutorAddressOperations {
            address_data: address_receiver,
            thread_pool: self.thread_pool.clone(),
        })
    }

    pub fn wait_for_executors(&self) -> BarrierWaitResult {
        self.barrier.wait()
    }
}

pub struct ExecutorAddressOperations<E: AsyncExecutor> {
    address_data: AddressConsumer<E>,
    thread_pool: Option<Arc<ScopedThreadPool>>,
}
impl<E: AsyncExecutor> ExecutorAddressOperations<E> {
    pub fn receive_packet(&self) -> Option<Packet<E::InputPacket>> {
        self.address_data.packets_queue.recv()
    }

    pub fn get_init_data(&self) -> &E::InitData {
        &self.address_data.init_data
    }

    pub fn spawn_executors<'a>(
        &'a self,
        count: usize,
        executor: impl Fn(usize) + Send + Sync + 'a,
    ) {
        if count == 1 {
            executor(0);
        } else {
            self.thread_pool
                .as_ref()
                .unwrap()
                .run_scoped_optional(count, |i| executor(i));
        }
    }
}

impl<E: AsyncExecutor> Drop for ExecutorAddressOperations<E> {
    fn drop(&mut self) {
        // Wait for the senders to terminate while ensuring that no more packets are waiting
        assert!(self.address_data.packets_queue.recv().is_none());
    }
}
