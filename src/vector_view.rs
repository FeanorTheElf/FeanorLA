use super::alg::*;
use std::marker::PhantomData;

pub trait VectorView<T> {
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;
}

pub trait VectorViewMut<T>: VectorView<T> {
    fn at_mut(&mut self, index: usize) -> &mut T;
}
