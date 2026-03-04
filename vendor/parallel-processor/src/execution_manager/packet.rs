use crate::execution_manager::objects_pool::{ObjectsPool, PoolObject, PoolObjectTrait};
use std::cell::UnsafeCell;
use std::ops::{Deref, DerefMut};

pub trait PacketTrait: PoolObjectTrait + Sync + Send {
    fn get_size(&self) -> usize;
}

pub struct Packet<T: PoolObjectTrait> {
    object: PoolObject<Box<T>>,
    _not_sync: std::marker::PhantomData<UnsafeCell<()>>,
}

pub struct PacketsPool<T: Sync + Send + 'static> {
    objects_pool: ObjectsPool<Box<T>>,
}

// Recursively implement the object trait for the pool, so it can be used recursively
impl<T: PoolObjectTrait> PoolObjectTrait for PacketsPool<T> {
    type InitData = (usize, T::InitData);

    fn allocate_new((cap, init_data): &Self::InitData) -> Self {
        Self::new(*cap, init_data.clone())
    }

    fn reset(&mut self) {}
}

impl<T: PoolObjectTrait> PacketsPool<T> {
    pub fn new(
        cache_size: usize,
        init_data: T::InitData,
        // mem_tracker: &Arc<MemoryTrackerManager>,
    ) -> Self {
        let objects_pool = ObjectsPool::new(cache_size, init_data);
        Self { objects_pool }
    }

    pub fn alloc_packet(&self) -> Packet<T> {
        let object = self.objects_pool.alloc_object();

        Packet {
            object,
            _not_sync: std::marker::PhantomData,
        }
    }
}

impl<T: PoolObjectTrait> Packet<T> {
    pub fn new_simple(data: T) -> Self {
        Packet {
            object: PoolObject::new_simple(Box::new(data)),
            _not_sync: std::marker::PhantomData,
        }
    }
}

impl<T: PoolObjectTrait> Deref for Packet<T> {
    type Target = T;
    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        self.object.deref()
    }
}

impl<T: PoolObjectTrait> DerefMut for Packet<T> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.object.deref_mut()
    }
}

impl PoolObjectTrait for () {
    type InitData = ();
    fn allocate_new(_init_data: &Self::InitData) -> Self {
        panic!("Cannot create () type as object!");
    }

    fn reset(&mut self) {
        panic!("Cannot reset () type as object!");
    }
}
impl PacketTrait for () {
    fn get_size(&self) -> usize {
        0
    }
}
