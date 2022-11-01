use super::*;

use std::marker::PhantomData;

pub struct VectorViewMapped<T, S, V, F>
    where F: for<'a> Fn(&'a T) -> &'a S, V: VectorView<T>
{
    base: V,
    from: PhantomData<T>,
    to: PhantomData<S>,
    f: F
}

impl<T, S, V, F> Clone for VectorViewMapped<T, S, V, F>
    where V: VectorView<T> + Clone, F: Clone + for<'a> Fn(&'a T) -> &'a S, V: VectorView<T>
{
    fn clone(&self) -> Self {
        VectorViewMapped {
            from: PhantomData,
            to: PhantomData,
            f: self.f.clone(),
            base: self.base.clone()
        }
    }
}

impl<T, S, V, F> Copy for VectorViewMapped<T, S, V, F>
    where V: VectorView<T> + Copy, F: Copy + for<'a> Fn(&'a T) -> &'a S, V: VectorView<T>
{}

impl<T, S, V, F> VectorViewMapped<T, S, V, F>
    where F: for<'a> Fn(&'a T) -> &'a S, V: VectorView<T>
{
    pub fn new(base: V, f: F) -> Self {
        VectorViewMapped {
            from: PhantomData,
            to: PhantomData,
            base: base,
            f: f
        }
    }
}

impl<T, S, V, F> VectorView<S> for VectorViewMapped<T, S, V, F>
    where F: for<'a> Fn(&'a T) -> &'a S, V: VectorView<T>
{
    type Subvector = Subvector<VectorViewMapped<T, S, V, F>, S>;

    fn len(&self) -> usize {
        self.base.len()
    }

    fn at(&self, i: usize) -> &S {
        (self.f)(self.base.at(i))
    }

    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        Self::Subvector::new(from, to, self)
    }
}
