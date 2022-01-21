use std::marker::PhantomData;

pub trait VectorView<T>: Sized {
    
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;

    fn assert_in_range(&self, index: usize) {
        assert!(index < self.len(), "Vector index {} out of range 0..{}", index, self.len());
    }

    fn into_owned(self) -> VectorOwned<T>
        where T: Clone
    {
        (0..self.len()).map(|i| self.at(i).clone()).collect()
    }
}

pub trait VectorViewMut<T>: VectorView<T> {
    fn at_mut(&mut self, index: usize) -> &mut T;
    fn swap(&mut self, i: usize, j: usize);
}

pub type VectorOwned<T> = Vec<T>;

impl<T> VectorView<T> for VectorOwned<T> {

    fn len(&self) -> usize {
        self.len()
    }

    fn at(&self, index: usize) -> &T {
        &self[index]
    }
}

impl<T> VectorViewMut<T> for VectorOwned<T> {

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        if i == j {
            return;
        }
        // I got an infinite recursion here before, so I want to be
        // explicit that I don't call `<Vec<T> as VectorViewMut<T>>::swap`
        let self_ref: &mut [T] = &mut self[..];
        self_ref.swap(i, j);
    }
}

impl<'a, T, V> VectorView<T> for &'a V
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        (**self).len()
    }

    fn at(&self, i: usize) -> &T {
        (**self).at(i)
    }
}

impl<'a, T, V> VectorView<T> for &'a mut V
    where V: VectorView<T>
{
    fn len(&self) -> usize {
        (**self).len()
    }

    fn at(&self, i: usize) -> &T {
        (**self).at(i)
    }
}

impl<'a, T, V> VectorViewMut<T> for &'a mut V
    where V: VectorViewMut<T>
{
    fn at_mut(&mut self, i: usize) -> &mut T {
        (**self).at_mut(i)
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        (**self).swap(fst, snd)
    }
}

//
// It would be nice to build this in the same way as `MatrixRowIter`,
// so directly store an element of type V which then may (probably will)
// be a reference in generic implementations. However, this causes a
// problem, as the lifetime of that object is then hidden within V,
// and we cannot declare the correct lifetime of the returned reference.
// In case of `MatrixRowIter`, this is no problem, as it returns `MatrixRow`s`
// instead of references.
//
pub struct VectorIter<'a, T, V> {
    vector: &'a V,
    current: usize,
    element: PhantomData<T>,
}

impl<'a, V, T> VectorIter<'a, T, V> 
    where V: 'a + VectorView<T>, T: 'a
{
    pub fn new(vector: &'a V) -> Self {
        VectorIter {
            vector: vector,
            current: 0,
            element: PhantomData
        }
    }
}

impl<'a, V, T> Iterator for VectorIter<'a, T, V> 
    where V: 'a + VectorView<T>, T: 'a
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.vector.len() {
            self.current += 1;
            return Some(self.vector.at(self.current - 1));
        } else {
            return None;
        }
    }
}

impl<'a, T, V> std::iter::FusedIterator for VectorIter<'a, T, V> 
    where V: 'a + VectorView<T>, T: 'a {}

impl<'a, T, V> Copy for VectorIter<'a, T, V> 
    where V: 'a + VectorView<T>, T: 'a {}

impl<'a, T, V> Clone for VectorIter<'a, T, V> 
    where V: 'a + VectorView<T>, T: 'a 
{
    fn clone(&self) -> Self {
        *self
    }
}

#[test]
fn test_vec_swap() {
    let mut a = vec![1, 2];
    <VectorOwned<i32> as VectorViewMut<i32>>::swap(&mut a, 0, 1);
    assert_eq!(2, *a.at(0));
}