use super::super::alg::*;
use super::vector_view::*;
use super::matrix_vector::*;
use super::constant::*;
use super::ops::*;
use super::vector::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, Add, SubAssign, Sub, MulAssign, RangeBounds, Bound, Index};

#[derive(Debug)]
pub struct Vector<V, T>
    where V: VectorView<T>
{
    data: V,
    element: PhantomData<T>
}

impl<V, T> Copy for Vector<V, T>
    where V: VectorView<T> + Copy
{}

impl<V, T> Clone for Vector<V, T>
    where V: VectorView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>
{
    pub fn new(vector: V) -> Self {
        Vector {
            data: vector,
            element: PhantomData
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn at(&self, i: usize) -> &T {
        self.data.at(i)
    }

    pub fn as_ref<'a>(&'a self) -> Vector<VectorRef<V, T>, T> {
        self.subvector(..)
    }

    pub fn subvector<'a, R>(&'a self, range: R) -> Vector<VectorRef<V, T>, T>
        where R: RangeBounds<usize>
    {
        let begin = match range.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let end = match range.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.len(),
        };
        return Vector::new(
            VectorRef::new(begin, end, &self.data)
        )
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>, T: Clone
{
    pub fn assign<W>(&mut self, rhs: Vector<W, T>) 
        where W: VectorView<T>
    {
        <T as MatrixAssign<_, _>>::assign_matrix(&mut self.as_mut().as_column_vector(), rhs.as_ref().as_column_vector())
    }
}

impl<V, T> Index<usize> for Vector<V, T>
    where V: VectorView<T>
{
    type Output = T;

    fn index(&self, i: usize) -> &T {
        self.at(i)
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>
{
    pub fn as_mut<'a>(&'a mut self) -> Vector<VectorRefMut<V, T>, T> {
        self.subvector_mut(..)
    }

    pub fn subvector_mut<'a, R>(&'a mut self, range: R) -> Vector<VectorRefMut<V, T>, T>
        where R: RangeBounds<usize>
    {
        let begin = match range.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let end = match range.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.len(),
        };
        return Vector::new(
            VectorRefMut::new(begin, end, &mut self.data)
        )
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>, T: Clone + std::fmt::Debug
{
    pub fn scale<R>(&mut self, rhs: &T, ring: &R) 
        where R: Ring<El = T>
    {
        <R as MatrixScale<_>>::scale_matrix(ring, rhs, &mut self.as_mut().as_column_vector());
    }

    pub fn add_assign<R, W>(&mut self, rhs: Vector<W, T>, ring: &R)
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixAddAssign<_, _>>::add_assign_matrix(ring, &mut self.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
    }

    pub fn sub_assign<R, W>(&mut self, rhs: Vector<W, T>, ring: &R)
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixAddAssign<_, _>>::sub_assign_matrix(ring, &mut self.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
    }
}
impl<V, T> Vector<V, T>
    where V: VectorView<T>, T: Clone + std::fmt::Debug
{
    pub fn add<R, W>(self, rhs: Vector<W, T>, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>, W: VectorView<T>
    {
        let mut result = self.to_owned();
        <R as MatrixAddAssign<_, _>>::add_assign_matrix(ring, &mut result.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
        return result;
    }

    pub fn sub<R, W>(self, rhs: Vector<W, T>, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>, W: VectorView<T>
    {
        let mut result = self.to_owned();
        <R as MatrixAddAssign<_, _>>::sub_assign_matrix(ring, &mut result.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
        return result;
    }
}

impl<V, W, T, U> PartialEq<Vector<W, U>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<U>, T: PartialEq<U>
{
    fn eq(&self, rhs: &Vector<W, U>) -> bool {
        for i in 0..self.len() {
            if self.at(i) != rhs.at(i) {
                return false;
            }
        }
        return true;
    }
}

impl<V, T> Eq for Vector<V, T>
    where V: VectorView<T>, T: Eq
{}

impl<V, T> Vector<V, T>
    where V: VectorView<T>, T: Clone
{
    pub fn to_owned(self) -> Vector<VectorOwned<T>, T> {
        Vector::new(self.data.to_owned())
    }
}

impl<T> Vector<VectorOwned<T>, T> {
    
    pub fn from_array<const L: usize>(data: [T; L]) -> Self {
        Vector::new(
            VectorOwned::from_array(data)
        )
    }

    pub fn from_fn<F>(len: usize, f: F) -> Self
        where F: FnMut(usize) -> T
    {
        Vector::new(VectorOwned::from_fn(len, f))
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>
{
    pub fn as_column_vector(self) -> ColumnVector<V, T> {
        ColumnVector::new(self.data)
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>
{
    pub fn as_row_vector(self) -> RowVector<V, T> {
        RowVector::new(self.data)
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>
{
    pub fn at_mut(&mut self, i: usize) -> &mut T {
        self.data.at_mut(i)
    }
}

impl<V, T, W> AddAssign<Vector<W, T>> for Vector<V, T>
    where V: VectorViewMut<T>, W: VectorView<T>, T: RingEl
{
    fn add_assign(&mut self, rhs: Vector<W, T>) {
        self.add_assign(rhs, &StaticRing::<T>::RING)
    }
}

impl<V, T, W> SubAssign<Vector<W, T>> for Vector<V, T>
    where V: VectorViewMut<T>, W: VectorView<T>, T: RingEl
{
    fn sub_assign(&mut self, rhs: Vector<W, T>) {
        self.sub_assign(rhs, &StaticRing::<T>::RING)
    }
}

impl<V, T, W> Add<Vector<W, T>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<T>, T: RingEl
{
    type Output = Vector<VectorOwned<T>, T>;

    fn add(self, rhs: Vector<W, T>) -> Self::Output {
        self.add(rhs, &StaticRing::<T>::RING)
    }
}

impl<V, T, W> Sub<Vector<W, T>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<T>, T: RingEl
{
    type Output = Vector<VectorOwned<T>, T>;

    fn sub(self, rhs: Vector<W, T>) -> Self::Output {
        self.sub(rhs, &StaticRing::<T>::RING)
    }
}

impl<V, T> MulAssign<T> for Vector<V, T>
    where V: VectorViewMut<T>, T: RingEl
{
    fn mul_assign(&mut self, rhs: T) {
        self.scale(&rhs, &StaticRing::<T>::RING)
    }
}

impl<T> Vector<VectorConstant<T>, T>
    where T: std::fmt::Debug + Clone
{
    pub fn zero_ring<R>(len: usize, ring: &R) -> Self 
        where R: Ring<El = T>
    {
        Vector::new(VectorConstant::new(len, ring.zero()))
    }
}

impl<T> Vector<VectorConstant<T>, T>
    where T: Zero
{
    pub fn zero(len: usize) -> Self {
        Vector::new(VectorConstant::new(len, T::zero()))
    }
}

impl<T> Vector<VectorOwned<T>, T>
    where T: Zero + One
{
    pub fn unit_vector(i: usize, len: usize) -> Self {
        Vector::new(VectorOwned::unit_vector(i, len))
    }
}

impl<V, T> std::fmt::Display for Vector<V, T> 
    where V: VectorView<T>, T: std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        write!(f, "{}", self.at(0))?;
        for i in 1..self.len() {
            write!(f, ", {}", self.at(i))?;
        }
        write!(f, "]")?;
        return Ok(());
    }
}
