pub mod compile_time_vector;
pub mod matrix_diagonal;
pub mod matrix_row_col;
pub mod constant_value_vector;

use std::marker::PhantomData;

pub trait VectorView<T>: Sized {

    type Subvector: VectorView<T>;
    
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;
    fn subvector(self, from: usize, to: usize) -> Self::Subvector;

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

    type SubvectorMut: VectorView<T>;

    fn at_mut(&mut self, index: usize) -> &mut T;
    fn swap(&mut self, i: usize, j: usize);

    ///
    /// Currently, there is no way to strengthen the associated type bound in
    /// a subtrait, i.e. require `Self::Subvector: VectorViewMut` in VectorViewMut.
    /// "Faking" this by adding a `where Self::Subvector: VectorViewMut` to the
    /// subtrait definition does not work:
    /// Adding such a where-clause then requires adding the corresponding where-clause
    /// to every use of VectorViewMut, and besides the inconvenience, this causes a
    /// recursion. In other words, we would have to add
    /// ```text
    ///     where V: VectorViewMut<T>,
    ///           V::Subvector: VectorViewMut<T>,
    ///           V::Subvector::Subvector: VectorViewMut<T>,
    ///           ...
    /// ```
    /// This is the best workaround I found. Just implement this function as the identity...
    /// 
    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut;
}

pub type VectorOwned<T> = Vec<T>;

impl<T> VectorView<T> for VectorOwned<T> {

    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        self.len()
    }

    fn at(&self, index: usize) -> &T {
        &self[index]
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<T> VectorViewMut<T> for VectorOwned<T> {

    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        if i == j {
            return;
        }
        <[T]>::swap(&mut self[..], i, j);
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
    }
}

impl<'a, T> VectorView<T> for &'a [T] {

    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    fn at(&self, index: usize) -> &T {
        &self[index]
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<'a, T> VectorView<T> for &'a mut [T] {

    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    fn at(&self, index: usize) -> &T {
        &self[index]
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<'a, T> VectorViewMut<T> for &'a mut [T] {

    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        if i == j {
            return;
        }
        <[T]>::swap(self, i, j);
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
    }
}

impl<'a, T, V> VectorView<T> for &'a V
    where V: VectorView<T>
{
    type Subvector = Subvector<&'a V, T>;

    fn len(&self) -> usize {
        (**self).len()
    }

    fn at(&self, i: usize) -> &T {
        (**self).at(i)
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<'a, T, V> VectorView<T> for &'a mut V
    where V: VectorView<T>
{
    type Subvector = Subvector<&'a mut V, T>;

    fn len(&self) -> usize {
        (**self).len()
    }

    fn at(&self, i: usize) -> &T {
        (**self).at(i)
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<'a, T, V> VectorViewMut<T> for &'a mut V
    where V: VectorViewMut<T>
{
    type SubvectorMut = Subvector<&'a mut V, T>;

    fn at_mut(&mut self, i: usize) -> &mut T {
        (**self).at_mut(i)
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        (**self).swap(fst, snd)
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
    }
}

#[derive(Debug, Clone, Copy)]
pub struct VectorArray<T, const N: usize>([T; N]);

impl<T, const N: usize> VectorArray<T, N> {

    pub const fn new(data: [T; N]) -> Self {
        VectorArray(data)
    }

    pub fn destruct(self) -> [T; N] {
        self.0
    }
}

impl<T, const N: usize> VectorView<T> for VectorArray<T, N> {

    type Subvector = Subvector<Self, T>;
    
    fn len(&self) -> usize {
        self.0.len()
    }

    fn at(&self, i: usize) -> &T {
        &self.0[i]
    }

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<T, const N: usize> VectorViewMut<T> for VectorArray<T, N> {
    
    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, i: usize) -> &mut T {
        &mut self.0[i]
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        self.0.swap(fst, snd)
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
    }
}

impl<T, const N: usize> std::iter::FromIterator<T> for VectorArray<T, N> {
    
    fn from_iter<I: IntoIterator<Item = T>>(into_iter: I) -> Self {
        let mut iter = into_iter.into_iter();
        let result = core::array::from_fn::<T, N, _>(|_| iter.next().unwrap());
        assert!(iter.next().is_none());
        return VectorArray(result);
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
    type Subvector = Self;

    fn len(&self) -> usize {
        self.to - self.from
    }

    fn at(&self, index: usize) -> &T {
        self.assert_in_range(index);
        self.view.at(index + self.from)
    }

    fn subvector(mut self, from: usize, to: usize) -> Self::Subvector {
        assert!(to <= self.len());
        self.to = self.from + to;
        self.from = self.from + from;
        return self;
    }
}

impl<V, T> VectorViewMut<T> for Subvector<V, T>
    where V: VectorViewMut<T>
{
    type SubvectorMut = Self;

    fn at_mut(&mut self, index: usize) -> &mut T {
        self.view.at_mut(index + self.from)
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        self.view.swap(i + self.from, j + self.from);
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
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

#[test]
fn test_vec_swap() {
    let mut a = vec![1, 2];
    <VectorOwned<i32> as VectorViewMut<i32>>::swap(&mut a, 0, 1);
    assert_eq!(2, *a.at(0));
}