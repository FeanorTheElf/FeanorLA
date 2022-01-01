use super::super::alg::*;
use super::vector_view::*;
use super::matrix_vector::*;
use super::diagonal::*;
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
    where V: VectorView<T> + Clone
{
    fn clone(&self) -> Self {
        Vector::new(self.data.clone())
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>
{
    pub const fn new(vector: V) -> Self {
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

    pub fn as_ref<'a>(&'a self) -> Vector<&'a V, T> {
        Vector::new(&self.data)
    }

    pub fn into_subvector<R>(self, range: R) -> Vector<Subvector<V, T>, T>
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
            Subvector::new(begin, end, self.data)
        )
    }

    pub fn subvector<'a, R>(&'a self, range: R) -> Vector<Subvector<&'a V, T>, T>
        where R: RangeBounds<usize>
    {
        Vector::new(&self.data).into_subvector(range)
    }

    pub fn iter<'a>(&'a self) -> VectorIter<'a, T, V> {
        VectorIter::new(&self.data)
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
    pub fn as_mut<'a>(&'a mut self) -> Vector<&'a mut V, T> {
        Vector::new(&mut self.data)
    }

    pub fn subvector_mut<'a, R>(&'a mut self, range: R) -> Vector<Subvector<&'a mut V, T>, T>
        where R: RangeBounds<usize>
    {
        self.as_mut().into_subvector(range)
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
        let mut result = self.into_owned();
        <R as MatrixAddAssign<_, _>>::add_assign_matrix(ring, &mut result.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
        return result;
    }

    pub fn sub<R, W>(self, rhs: Vector<W, T>, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>, W: VectorView<T>
    {
        let mut result = self.into_owned();
        <R as MatrixAddAssign<_, _>>::sub_assign_matrix(ring, &mut result.as_mut().as_column_vector(), rhs.as_ref().as_column_vector());
        return result;
    }

    pub fn neg<R>(self, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>
    {
        let mut result = self.into_owned();
        <R as MatrixScale<_>>::negate_matrix(ring, &mut result.as_mut().as_column_vector());
        return result;
    }

    pub fn eq<R, W>(self, rhs: Vector<W, T>, ring: &R) -> bool
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixEq<_, _>>::eq_matrix(ring, &self.as_ref().as_column_vector(), rhs.as_ref().as_column_vector())
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
    pub fn into_owned(self) -> Vector<VectorOwned<T>, T> {
        Vector::new(self.data.into_owned())
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

    pub fn unit_vector_ring<R>(i: usize, len: usize, ring: &R) -> Self 
        where R: Ring<El = T>, T: std::fmt::Debug + Clone
    {
        let mut result = Vector::zero_ring(len, ring).into_owned();
        *result.at_mut(i) = ring.one();
        return result;
    }

    pub fn raw_data(self) -> Box<[T]> {
        self.data.raw_data()
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>
{
    pub fn as_column_vector(self) -> ColumnVector<V, T> {
        ColumnVector::new(self.data)
    }

    pub fn as_row_vector(self) -> RowVector<V, T> {
        RowVector::new(self.data)
    }

    pub fn as_diag_matrix(self, diag_index: i64, zero: T) -> DiagonalMatrix<V, T> {
        DiagonalMatrix::new(self.data, diag_index, zero)
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
        self.add_assign(rhs, &T::RING)
    }
}

impl<V, T, W> SubAssign<Vector<W, T>> for Vector<V, T>
    where V: VectorViewMut<T>, W: VectorView<T>, T: RingEl
{
    fn sub_assign(&mut self, rhs: Vector<W, T>) {
        self.sub_assign(rhs, &T::RING)
    }
}

impl<V, T, W> Add<Vector<W, T>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<T>, T: RingEl
{
    type Output = Vector<VectorOwned<T>, T>;

    fn add(self, rhs: Vector<W, T>) -> Self::Output {
        self.add(rhs, &T::RING)
    }
}

impl<V, T, W> Sub<Vector<W, T>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<T>, T: RingEl
{
    type Output = Vector<VectorOwned<T>, T>;

    fn sub(self, rhs: Vector<W, T>) -> Self::Output {
        self.sub(rhs, &T::RING)
    }
}

impl<V, T> MulAssign<T> for Vector<V, T>
    where V: VectorViewMut<T>, T: RingEl
{
    fn mul_assign(&mut self, rhs: T) {
        self.scale(&rhs, &T::RING)
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

impl<T> Vector<VectorConstant<T>, T> {

    pub fn constant(len: usize, c: T) -> Self {
        Vector::new(VectorConstant::new(len, c))
    }
}

impl<T> Vector<VectorConstant<T>, T>
    where T: Zero
{
    pub fn zero(len: usize) -> Self {
        Self::constant(len, T::zero())
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
