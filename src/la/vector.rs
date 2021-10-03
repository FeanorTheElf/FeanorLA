use super::vector_view::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct Subvector<V, T> 
    where V: VectorView<T>
{
    view: V,
    from: usize,
    to: usize,
    element: PhantomData<T>
}

impl<V, T> Copy for Subvector<V, T> 
    where V: VectorView<T> + Copy {}

impl<V, T> Clone for Subvector<V, T>  
    where V: VectorView<T> + Clone
{
    fn clone(&self) -> Self {
        Subvector {
            view: self.view.clone(),
            from: self.from,
            to: self.to,
            element: PhantomData
        }
    }
}

impl<V, T> VectorView<T> for Subvector<V, T> 
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        self.to - self.from
    }

    fn at(&self, index: usize) -> &T {
        self.view.at(index + self.from)
    }
}

impl<V, T> VectorViewMut<T> for Subvector<V, T>
    where V: VectorViewMut<T>
{
    fn at_mut(&mut self, index: usize) -> &mut T {
        self.view.at_mut(index + self.from)
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        self.view.swap(i + self.from, j + self.from);
    }
}

impl<V, T> Subvector<V, T> 
    where V: VectorView<T>
{
    pub fn new(from: usize, to: usize, vector: V) -> Self {
        assert!(from <= to);
        assert!(to <= vector.len());
        Subvector {
            from: from,
            to: to,
            view: vector,
            element: PhantomData
        }
    }
}
