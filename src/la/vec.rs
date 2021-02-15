use super::super::alg::*;
use super::vector_view::*;
use super::vector::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, Add, Mul, MulAssign, RangeBounds, Bound};

#[derive(Debug, Clone, Copy)]
pub struct Vector<V, T>
    where V: VectorView<T>
{
    data: V,
    element: PhantomData<T>
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
    
    pub fn scale<U>(&mut self, rhs: U) 
        where T: MulAssign<U>, U: Clone
    {
        for i in 0..self.len() {
            *self.at_mut(i) *= rhs.clone();
        }
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
    pub fn to_owned(&self) -> Vector<VectorOwned<T>, T> {
        Vector::new(
            VectorOwned::from_fn(self.len(), |i| self.at(i).clone())
        )
    }
}

impl<T> Vector<VectorOwned<T>, T> {
    
    pub fn from_array<const L: usize>(data: [T; L]) -> Self {
        Vector::new(
            VectorOwned::from_array(data)
        )
    }
}

impl<V, T> Vector<V, T>
    where V: VectorViewMut<T>
{
    pub fn at_mut(&mut self, i: usize) -> &mut T {
        self.data.at_mut(i)
    }
}

impl<V, T, W, U> AddAssign<Vector<W, U>> for Vector<V, T>
    where V: VectorViewMut<T>, W: VectorView<U>, U: Clone, T: AddAssign<U>
{
    fn add_assign(&mut self, rhs: Vector<W, U>) {
        assert_eq!(self.len(), rhs.len());
        for i in 0..self.len() {
            *self.at_mut(i) += rhs.at(i).clone();
        }
    }
}

impl<V, T, W, U> Add<Vector<W, U>> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<U>, U: Clone, T: Add<U> + Clone
{
    type Output = Vector<VectorOwned<T::Output>, T::Output>;

    fn add(self, rhs: Vector<W, U>) -> Self::Output {
        // TODO: Implementation that does not require U, T: Clone
        assert_eq!(self.len(), rhs.len());
        Vector::new(
            VectorOwned::from_fn(self.len(), |i| self.at(i).clone() + rhs.at(i).clone())
        )
    }
}

impl<V, T, U> MulAssign<U> for Vector<V, T>
    where V: VectorViewMut<T>, U: Clone, T: MulAssign<U>
{
    fn mul_assign(&mut self, rhs: U) {
        for i in 0..self.len() {
            *self.at_mut(i) *= rhs.clone();
        }
    }
}

impl<V, T, U> Mul<U> for Vector<V, T>
    where V: VectorView<T>, U: Clone, T: Mul<U> + Clone
{
    type Output = Vector<VectorOwned<T::Output>, T::Output>;

    fn mul(self, rhs: U) -> Self::Output {
        // TODO: Implementation that does not require U, T: Clone
        Vector::new(
            VectorOwned::from_fn(self.len(), |i| self.at(i).clone() * rhs.clone())
        )
    }
}

impl<T> Vector<VectorOwned<T>, T>
    where T: Zero
{
    pub fn zero(len: usize) -> Self {
        Vector::new(VectorOwned::zero(len))
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