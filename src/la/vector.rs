use super::super::alg::*;
use super::vector_view::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct VectorRef<'a, V, T> 
    where V: VectorView<T>
{
    view: &'a V,
    from: usize,
    to: usize,
    element: PhantomData<T>
}

#[derive(Debug)]
pub struct VectorRefMut<'a, V, T> 
    where V: VectorViewMut<T>
{
    view: &'a mut V,
    from: usize,
    to: usize,
    element: PhantomData<T>
}

impl<'a, V, T> Copy for VectorRef<'a, V, T> 
    where V: VectorView<T> {}

impl<'a, V, T> Clone for VectorRef<'a, V, T>  
    where V: VectorView<T>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, V, T> VectorView<T> for VectorRef<'a, V, T> 
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        self.to - self.from
    }

    fn at(&self, index: usize) -> &T {
        self.view.at(index + self.from)
    }
}

impl<'a, V, T> VectorView<T> for VectorRefMut<'a, V, T>
    where V: VectorViewMut<T>
{

    fn len(&self) -> usize {
        self.to - self.from
    }

    fn at(&self, index: usize) -> &T {
        self.view.at(index + self.from)
    }
}

impl<'a, V, T> VectorViewMut<T> for VectorRefMut<'a, V, T>
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

impl<'a, V, T> VectorRef<'a, V, T> 
    where V: VectorView<T>
{
    pub fn new(from: usize, to: usize, vector: &'a V) -> Self {
        assert!(from < to);
        assert!(to <= vector.len());
        VectorRef {
            from: from,
            to: to,
            view: vector,
            element: PhantomData
        }
    }
}

impl<'a, V, T> VectorRefMut<'a, V, T> 
    where V: VectorViewMut<T>
{
    pub fn new(from: usize, to: usize, vector: &'a mut V) -> Self {
        assert!(from < to);
        assert!(to <= vector.len());
        VectorRefMut {
            from: from,
            to: to,
            view: vector,
            element: PhantomData
        }
    }
}

#[derive(Debug)]
pub struct VectorRestriction<V, T>
    where V: VectorView<T>
{
    base: V,
    from: usize,
    to: usize,
    element: PhantomData<T>
}


impl<V, T> VectorRestriction<V, T>
    where V: VectorView<T>
{
    pub fn restrict(vector: V, from: usize, to: usize) -> Self {
        vector.assert_in_range(from);
        assert!(to <= vector.len());
        assert!(from < to);
        VectorRestriction {
            base: vector,
            from: from,
            to: to,
            element: PhantomData
        }
    }
}

impl<V, T> VectorView<T> for VectorRestriction<V, T>
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        self.to - self.from
    }

    fn at(&self, i: usize) -> &T {
        self.base.at(i + self.from)
    }
}

impl<V, T> VectorViewMut<T> for VectorRestriction<V, T>
    where V: VectorViewMut<T>
{
    fn at_mut(&mut self, i: usize) -> &mut T {
        self.base.at_mut(i + self.from)
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        self.base.swap(i + self.from, j + self.from);
    }
}
