use super::vector_view::*;
use std::marker::PhantomData;
use super::alg::*;

#[derive(Debug, Clone)]
pub struct VectorOwned<T> {
    data: Box<[T]>,
}

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

impl<T> VectorView<T> for VectorOwned<T> {

    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<T> VectorViewMut<T> for VectorOwned<T> {

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }
}

impl<'a, V, T> VectorView<T> for VectorRef<'a, V, T> 
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        self.from - self.to
    }

    fn at(&self, index: usize) -> &T {
        self.view.at(index + self.from)
    }
}

impl<'a, V, T> VectorView<T> for VectorRefMut<'a, V, T>
    where V: VectorViewMut<T>
{

    fn len(&self) -> usize {
        self.from - self.to
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
}

impl<T> VectorOwned<T> {

    pub fn new(data: Box<[T]>) -> VectorOwned<T> {
        VectorOwned { data }
    }
}

impl<T> VectorOwned<T>
    where T: Zero + One 
{
    pub fn unit_vector(index: usize, len: usize) -> VectorOwned<T> {
        assert!(index < len, "Expected index of 1-entry in unit vector to be within 0 and {}, got {}", len, index);
        VectorOwned::new((0..len).map(|i| if i == index { T::one() } else { T::zero() }).collect::<Vec<T>>().into_boxed_slice())
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
}
