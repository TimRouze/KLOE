use std::mem::ManuallyDrop;
use std::ops::{Deref, DerefMut};
use std::sync::Arc;

use crossbeam::queue::ArrayQueue;

pub trait PoolObjectTrait: Send + Sync + 'static {
    type InitData: Clone + Sync + Send;

    fn allocate_new(init_data: &Self::InitData) -> Self;
    fn reset(&mut self);
}

impl<T: PoolObjectTrait> PoolObjectTrait for Box<T> {
    type InitData = T::InitData;

    fn allocate_new(init_data: &Self::InitData) -> Self {
        Box::new(T::allocate_new(init_data))
    }
    fn reset(&mut self) {
        T::reset(self);
    }
}

#[derive(Clone)]
pub(crate) struct PoolReturner<T> {
    returner: Arc<ArrayQueue<T>>,
}

// Non blocking pool that caches a fixed number of elements and creates new elements as requested
pub struct ObjectsPool<T: Sync + Send + 'static> {
    queue: Arc<ArrayQueue<T>>,
    allocate_fn: Box<dyn (Fn() -> T) + Sync + Send>,
}

impl<T: PoolObjectTrait> PoolReturner<T> {
    fn return_element(&self, mut el: T) {
        el.reset();
        // Push it into the queue if the capacity allows it, else deallocates
        let _ = self.returner.push(el);
    }
}

impl<T: PoolObjectTrait> ObjectsPool<T> {
    pub fn new(cache_size: usize, init_data: T::InitData) -> Self {
        let queue = Arc::new(ArrayQueue::new(cache_size));

        Self {
            queue,
            allocate_fn: Box::new(move || T::allocate_new(&init_data)),
        }
    }

    pub fn alloc_object(&self) -> PoolObject<T> {
        match self.queue.pop() {
            Some(el) => {
                // Element found in queue, return it
                PoolObject::from_element(el, self)
            }
            None => PoolObject::from_element((self.allocate_fn)(), self),
        }
    }
}

pub struct PoolObject<T: PoolObjectTrait> {
    pub(crate) value: ManuallyDrop<T>,
    pub(crate) returner: Option<PoolReturner<T>>,
}

impl<T: PoolObjectTrait> PoolObject<T> {
    fn from_element(value: T, pool: &ObjectsPool<T>) -> Self {
        Self {
            value: ManuallyDrop::new(value),
            returner: Some(PoolReturner {
                returner: pool.queue.clone(),
            }),
        }
    }
}

impl<T: PoolObjectTrait> PoolObject<T> {
    pub fn new_simple(value: T) -> Self {
        Self {
            value: ManuallyDrop::new(value),
            returner: None,
        }
    }
}

impl<T: PoolObjectTrait> Deref for PoolObject<T> {
    type Target = T;
    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        self.value.deref()
    }
}

impl<T: PoolObjectTrait> DerefMut for PoolObject<T> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.value.deref_mut()
    }
}

impl<T: PoolObjectTrait> Drop for PoolObject<T> {
    fn drop(&mut self) {
        if let Some(returner) = &self.returner {
            returner.return_element(unsafe { ManuallyDrop::take(&mut self.value) });
        } else {
            unsafe { ManuallyDrop::drop(&mut self.value) }
        }
    }
}
