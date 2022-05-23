use super::super::ring::*;
use super::super::primitive::*;
use super::ops::*;

pub use super::vector_view::*;
pub use super::vector_view::compile_time_vector::*;
pub use super::vector_view::constant_value_vector::*;
pub use super::vector_view::matrix_diagonal::*;
pub use super::vector_view::matrix_row_col::*;
use super::matrix_view::matrix_vector::*;

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

    pub fn into_subvector<R>(self, range: R) -> Vector<V::Subvector, T>
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
        assert!(begin <= end);
        assert!(end <= self.len());
        return Vector::new(
            self.data.subvector(begin, end)
        )
    }

    ///
    /// Computes the subvector containing all elements with index within range 
    /// that exist in this vector. In particular, this behaves like `into_subvector`
    /// if the given range is completely contained in self.
    /// 
    /// # Examples
    /// 
    /// ```
    /// # use feanor_la::la::vec::*;
    /// assert_eq!(Vector::<_, i64>::from_array([]), Vector::from_array([1, 2]).into_subvector_intersect(2..));
    /// assert_eq!(Vector::from_array([2]), Vector::from_array([1, 2, 3]).into_subvector_intersect(1..2));
    /// assert_eq!(Vector::<_, i64>::from_array([]), Vector::from_array([1, 2]).into_subvector_intersect(1..0));
    /// ```
    /// 
    pub fn into_subvector_intersect<R>(self, range: R) -> Vector<V::Subvector, T>
        where R: RangeBounds<usize>
    {
        let begin = std::cmp::max(0, match range.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        });
        let end = std::cmp::min(self.len(), match range.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.len(),
        });
        if begin >= end {
            return Vector::new(
                self.data.subvector(0, 0)
            );
        } else {
            return Vector::new(
                self.data.subvector(begin, end)
            );
        }
    }

    pub fn subvector<'a, R>(&'a self, range: R) -> Vector<<&'a V as VectorView<T>>::Subvector, T>
        where R: RangeBounds<usize>
    {
        Vector::new(&self.data).into_subvector(range)
    }

    ///
    /// Look at `into_subvector_intersect()` to see how this differs from `subvector_mut()`.
    /// 
    pub fn subvector_intersect<'a, R>(&'a self, range: R) -> Vector<<&'a V as VectorView<T>>::Subvector, T>
        where R: RangeBounds<usize>
    {
        Vector::new(&self.data).into_subvector_intersect(range)
    }

    pub fn iter<'a>(&'a self) -> VectorIter<'a, T, V> {
        VectorIter::new(&self.data)
    }

    pub fn raw_data(self) -> V {
        self.data
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>, T: Clone
{
    pub fn assign<W>(&mut self, rhs: Vector<W, T>) 
        where W: VectorView<T>
    {
        <T as MatrixAssign<_, _>>::assign_matrix(&mut ColumnVector::new(self.as_mut().data), ColumnVector::new(rhs.as_ref().data))
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

    pub fn into_subvector_mut<R>(self, range: R) -> Vector<V::SubvectorMut, T>
        where R: RangeBounds<usize>
    {
        Vector::new(V::cast_subvector(self.into_subvector(range).data))
    }

    ///
    /// Look at `into_subvector_intersect()` to see how this differs from `subvector_mut()`.
    /// 
    pub fn into_subvector_mut_intersect<R>(self, range: R) -> Vector<V::SubvectorMut, T>
        where R: RangeBounds<usize>
    {
        Vector::new(V::cast_subvector(self.into_subvector_intersect(range).data))
    }

    pub fn subvector_mut<'a, R>(&'a mut self, range: R) -> Vector<<&'a mut V as VectorViewMut<T>>::SubvectorMut, T>
        where R: RangeBounds<usize>
    {
        self.as_mut().into_subvector(range)
    }

    ///
    /// Look at `into_subvector_intersect()` to see how this differs from `subvector_mut()`.
    /// 
    pub fn subvector_mut_intersect<'a, R>(&'a mut self, range: R) -> Vector<<&'a mut V as VectorViewMut<T>>::SubvectorMut, T>
        where R: RangeBounds<usize>
    {
        self.as_mut().into_subvector_intersect(range)
    }
}

impl<'a, V, T> Vector<&'a mut V, T>
    where V: VectorViewMut<T>
{
    pub fn borrow<'b: 'a>(&'b mut self) -> Vector<&'b mut V, T> {
        Vector::new(&mut *self.data)
    }
}

impl<'a, T> Vector<&'a mut [T], T>
{
    pub fn split(self, mid: usize) -> (Self, Self) {
        let (l, r) = self.data.split_at_mut(mid);
        (Vector::new(l), Vector::new(r))
    }
    
    pub fn borrow<'b>(&'b mut self) -> Vector<&'b mut [T], T>
        where 'a: 'b
    {
        Vector::new(&mut *self.data)
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>, T: Clone + std::fmt::Debug
{
    pub fn scale<R>(&mut self, rhs: &T, ring: &R) 
        where R: Ring<El = T>
    {
        <R as MatrixScale<_>>::scale_matrix(ring, rhs, &mut ColumnVector::new(self.as_mut().data));
    }

    pub fn add_assign<R, W>(&mut self, rhs: Vector<W, T>, ring: &R)
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixAddAssign<_, _>>::add_assign_matrix(ring, &mut ColumnVector::new(self.as_mut().data), ColumnVector::new(rhs.as_ref().data));
    }

    pub fn sub_assign<R, W>(&mut self, rhs: Vector<W, T>, ring: &R)
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixAddAssign<_, _>>::sub_assign_matrix(ring, &mut ColumnVector::new(self.as_mut().data), ColumnVector::new(rhs.as_ref().data));
    }
}

impl<V, T> Vector<V, T>
    where V: VectorView<T>, T: Clone + std::fmt::Debug
{
    pub fn add<R, W>(self, rhs: Vector<W, T>, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>, W: VectorView<T>
    {
        let mut result = self.into_owned();
        <R as MatrixAddAssign<_, _>>::add_assign_matrix(ring, &mut ColumnVector::new(result.as_mut().data), ColumnVector::new(rhs.as_ref().data));
        return result;
    }

    pub fn sub<R, W>(self, rhs: Vector<W, T>, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>, W: VectorView<T>
    {
        let mut result = self.into_owned();
        <R as MatrixAddAssign<_, _>>::sub_assign_matrix(ring, &mut ColumnVector::new(result.as_mut().data), ColumnVector::new(rhs.as_ref().data));
        return result;
    }

    pub fn neg<R>(self, ring: &R) -> Vector<VectorOwned<T>, T>
        where R: Ring<El = T>
    {
        let mut result = self.into_owned();
        <R as MatrixScale<_>>::negate_matrix(ring, &mut ColumnVector::new(result.as_mut().data));
        return result;
    }

    pub fn eq<R, W>(self, rhs: Vector<W, T>, ring: &R) -> bool
        where R: Ring<El = T>, W: VectorView<T>
    {
        <R as MatrixEq<_, _>>::eq_matrix(ring, ColumnVector::new(self.as_ref().data), ColumnVector::new(rhs.as_ref().data))
    }

    pub fn l2_norm_square<R>(self, ring: &R) -> R::El
        where R: Ring<El = T>
    {
        <R as MatrixFrobenius<_>>::calc_matrix_frobenius_norm_square(ring, ColumnVector::new(self.as_ref().data))
    }
}

impl<V, W, T, U> PartialEq<Vector<W, U>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<U>, T: PartialEq<U>
{
    fn eq(&self, rhs: &Vector<W, U>) -> bool {
        assert_eq!(self.len(), rhs.len());
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

impl<T, const N: usize> Vector<VectorArray<T, N>, T> {
    
    pub fn from_array(data: [T; N]) -> Self {
        Vector::new(VectorArray::new(data))
    }
}

impl<T> Vector<VectorOwned<T>, T> {
    
    pub fn from_fn<F>(len: usize, f: F) -> Self
        where F: FnMut(usize) -> T
    {
        Vector::new((0..len).map(f).collect())
    }

    pub fn unit_vector_ring<R>(i: usize, len: usize, ring: &R) -> Self 
        where R: Ring<El = T>, T: std::fmt::Debug + Clone
    {
        let mut result = Vector::zero_ring(len, ring).into_owned();
        *result.at_mut(i) = ring.one();
        return result;
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

impl<T> Vector<ConstantVectorView<T>, T>
    where T: std::fmt::Debug + Clone
{
    pub fn zero_ring<R>(len: usize, ring: &R) -> Self 
        where R: Ring<El = T>
    {
        Vector::new(ConstantVectorView::new(len, ring.zero()))
    }
}

impl<T> Vector<ConstantVectorView<T>, T> {

    pub fn constant(len: usize, c: T) -> Self {
        Vector::new(ConstantVectorView::new(len, c))
    }
}

impl<T> Vector<ConstantVectorView<T>, T>
    where T: Zero
{
    pub fn zero(len: usize) -> Self {
        Self::constant(len, T::zero())
    }
}

impl<T> Vector<UnitVectorView<T>, T>
    where T: Zero + One
{
    pub fn unit_vector(i: usize, len: usize) -> Self {
        Vector::new(UnitVectorView::new(len, i, T::zero(), T::one()))
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

#[test]
fn test_subvector_subvector_same() {
    let v = Vector::from_array([1, 2]);
    let a: Vector<Subvector<&VectorArray<i64, 2>, i64>, i64> = v.subvector(..);
    let b: Vector<Subvector<&VectorArray<i64, 2>, i64>, i64> = a.into_subvector(..);
}