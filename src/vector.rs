use super::arith::*;
use super::indexed::{Indexed, IndexedMut};
use std::ops::{
    Add, AddAssign, Index, IndexMut, Neg, RangeBounds, Bound,
    Sub, SubAssign,
};

pub trait VectorView<T> {
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;
}

pub trait VectorViewMut<T>: VectorView<T> {
    fn at_mut(&mut self, index: usize) -> &mut T;
}

#[derive(Debug, Clone)]
pub struct Vector<T> {
    data: Box<[T]>,
}

#[derive(Debug, Clone, Copy)]
pub struct VectorRef<'a, T> {
    data: &'a [T],
}

#[derive(Debug)]
pub struct VectorRefMut<'a, T> {
    data: &'a mut [T],
}

impl<T> VectorView<T> for Vector<T> {

    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<T> VectorViewMut<T> for Vector<T> {

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }
}

impl<'a, T> VectorView<T> for VectorRef<'a, T> {

    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<'a, T> VectorView<T> for VectorRefMut<'a, T> {

    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<'a, T> VectorViewMut<T> for VectorRefMut<'a, T> {

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }
}

impl<T> Vector<T> {

    pub fn new(data: Box<[T]>) -> Vector<T> {
        Vector { data }
    }

    pub fn as_ref(&self) -> VectorRef<T> {
        VectorRef {
            data: &self.data
        }
    }

    pub fn as_mut(&mut self) -> VectorRefMut<T> {
        VectorRefMut {
            data: &mut self.data
        }
    }
}

impl<T> Vector<T>
    where T: Zero + One 
{
    pub fn unit_vector(index: usize, len: usize) -> Vector<T> {
        assert!(index < len, "Expected index of 1-entry in unit vector to be within 0 and {}, got {}", len, index);
        Vector::new((0..len).map(|i| if i == index { T::one() } else { T::zero() }).collect::<Vec<T>>().into_boxed_slice())
    }
}

impl<T> Vector<T>
    where T: Zero 
{
    pub fn zero(len: usize) -> Vector<T> {
        Vector::new((0..len).map(|_| T::zero()).collect::<Vec<T>>().into_boxed_slice())
    }
}

impl<'a, T> VectorRef<'a, T> {

    pub fn create(data: &'a [T]) -> VectorRef<T> {
        VectorRef { data }
    }
}

impl<'a, T> VectorRefMut<'a, T> {

    pub fn create(data: &'a mut [T]) -> VectorRefMut<T> {
        VectorRefMut { data }
    }

    pub fn as_const(&self) -> VectorRef<T> {
        VectorRef {
            data: self.data
        }
    }

    pub fn into_const(self) -> VectorRef<'a, T> {
        VectorRef {
            data: self.data
        }
    }

    pub fn into_subrange<R: RangeBounds<usize>>(self, range: R) -> VectorRefMut<'a, T> {
        VectorRefMut {
            data: &mut self.data[get_lower_index(&range, self.len())..get_upper_index(&range, self.len())]
        }
    }

    pub fn into_slice(self) -> &'a mut [T] {
        self.data
    }
}

// ================== Arithmetic operations ==================

impl<T, V> Add<V> for Vector<T> 
    where T: for<'b> AddAssign<&'b T> + Clone, V: VectorView<T>
{
    type Output = Vector<T>;

    fn add(mut self, rhs: V) -> Self::Output {
        self += rhs;
        return self;
    }
}

impl<'a, T, V> Add<V> for VectorRef<'a, T> 
    where T: Add<T, Output = T> + Clone, V: VectorView<T>
{
    type Output = Vector<T>;

    fn add(self, rhs: V) -> Self::Output {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Expected the lengths of summed vectors to be equal, but got {} and {}",
            self.len(),
            rhs.len()
        );
        let values = (0..self.len()).map(|i| self.at(i).clone() + rhs.at(i).clone());
        Vector::new(values.collect::<Vec<T>>().into_boxed_slice())
    }
}

impl<'a, T, V> Add<V> for VectorRefMut<'a, T> 
    where T: Add<T, Output = T> + Clone, V: VectorView<T>
{
    type Output = Vector<T>;

    fn add(self, rhs: V) -> Self::Output {
        self.into_const().add(rhs)
    }
}

impl<T, V> Sub<V> for Vector<T> 
    where T: for<'b> SubAssign<&'b T> + Clone, V: VectorView<T>
{
    type Output = Vector<T>;

    fn sub(mut self, rhs: V) -> Self::Output {
        self -= rhs;
        return self;
    }
}

impl<'a, V, T> Sub<V> for VectorRef<'a, T> 
    where T: Sub<T, Output = T> + Copy, V: VectorView<T>
{
    type Output = Vector<T>;

    fn sub(self, rhs: V) -> Self::Output {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Expected the lengths of subtracted vectors to be equal, but got {} and {}",
            self.len(),
            rhs.len()
        );
        let values = (0..self.len()).map(|i| self.at(i).clone() - rhs.at(i).clone());
        Vector::new(values.collect::<Vec<T>>().into_boxed_slice())
    }
}

impl<'a, V, T> Sub<V> for VectorRefMut<'a, T> 
    where T: Sub<T, Output = T> + Copy, V: VectorView<T>
{
    type Output = Vector<T>;

    fn sub(self, rhs: V) -> Self::Output {
        self.into_const().sub(rhs)
    }
}

impl<T> Neg for Vector<T> 
    where T: Neg + Clone
{
    type Output = Vector<T::Output>;

    fn neg(self) -> Self::Output {
        self.as_ref().neg()
    }
}

impl<'a, T> Neg for VectorRef<'a, T> 
    where T: Neg + Clone
{
    type Output = Vector<T::Output>;

    fn neg(self) -> Self::Output {
        Vector::new((0..self.len()).map(|i| -self.at(i).clone()).collect::<Vec<T::Output>>().into_boxed_slice())
    }
}

impl<'a, T> Neg for VectorRefMut<'a, T> 
    where T: Neg + Clone
{
    type Output = Vector<T::Output>;

    fn neg(self) -> Self::Output {
        self.into_const().neg()
    }
}

impl<'a, T> IntoIterator for VectorRef<'a, T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<'a, T> IntoIterator for VectorRefMut<'a, T> {
    type Item = &'a mut T;
    type IntoIter = std::slice::IterMut<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter_mut()
    }
}

impl<T, V> AddAssign<V> for Vector<T>
where
    T: for<'b> AddAssign<&'b T>, V: VectorView<T>
{
    fn add_assign(&mut self, rhs: V) {
        self.as_mut().add_assign(rhs);
    }
}

impl<'a, T, V> AddAssign<V> for VectorRefMut<'a, T>
where
    T: for<'b> AddAssign<&'b T>, V: VectorView<T>
{
    fn add_assign(&mut self, rhs: V) {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Expected the lengths of summed vectors to be equal, but got {} and {}",
            self.len(),
            rhs.len()
        );
        for i in 0..self.len() {
            (*self.at_mut(i)).add_assign(rhs.at(i));
        }
    }
}

impl<T, V> SubAssign<V> for Vector<T>
where
    T: for<'b> SubAssign<&'b T>, V: VectorView<T>
{
    fn sub_assign(&mut self, rhs: V) {
        self.as_mut().sub_assign(rhs)
    }
}

impl<'a, T, V> SubAssign<V> for VectorRefMut<'a, T>
where
    T: for<'b> SubAssign<&'b T>, V: VectorView<T>
{
    fn sub_assign(&mut self, rhs: V) {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Expected the lengths of subtracted vectors to be equal, but got {} and {}",
            self.len(),
            rhs.len()
        );
        for i in 0..self.len() {
            (*self.at_mut(i)).sub_assign(rhs.at(i));
        }
    }
}

impl<T, V> PartialEq<V> for Vector<T>
where
    T: PartialEq, V: VectorView<T>
{
    fn eq(&self, rhs: &V) -> bool {
        self.as_ref() == *rhs
    }
}

impl<'a, T, V> PartialEq<V> for VectorRef<'a, T>
where
    T: PartialEq, V: VectorView<T>
{
    fn eq(&self, rhs: &V) -> bool {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Expected the lengths of compared vectors to be equal, but got {} and {}",
            self.len(),
            rhs.len()
        );
        for i in 0..self.len() {
            if self.at(i) != rhs.at(i) {
                return false;
            }
        }
        return true;
    }
}

impl<'a, T, V> PartialEq<V> for VectorRefMut<'a, T>
where
    T: PartialEq, V: VectorView<T>
{
    fn eq(&self, rhs: &V) -> bool {
        self.as_const() == *rhs
    }
}

impl<T> Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.at(index)
    }
}

impl<T> IndexMut<usize> for Vector<T> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.at_mut(index)
    }
}

impl<'a, T> Index<usize> for VectorRef<'a, T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.at(index)
    }
}

impl<'a, T> Index<usize> for VectorRefMut<'a, T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.at(index)
    }
}

impl<'a, T> IndexMut<usize> for VectorRefMut<'a, T> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.at_mut(index)
    }
}

fn get_lower_index<R: RangeBounds<usize>>(range: &R, size: usize) -> usize {
    match range.start_bound() {
        Bound::Excluded(i) => *i + 1,
        Bound::Included(i) => *i,
        Bound::Unbounded => 0
    }
}

fn get_upper_index<R: RangeBounds<usize>>(range: &R, size: usize) -> usize {
    match range.end_bound() {
        Bound::Excluded(i) => {
            assert!(*i <= size, "Index {} out of range for size {}", i - 1, size);
            *i
        },
        Bound::Included(i) => {
            assert!(*i < size, "Index {} out of range for size {}", i, size);
            *i + 1
        },
        Bound::Unbounded => size
    }
}

impl<'b, T: 'b, R: RangeBounds<usize>> Indexed<'b, R> for Vector<T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, index: R) -> Self::Output {
        self.as_ref().get(index)
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>> Indexed<'b, R> for VectorRef<'a, T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, index: R) -> Self::Output {
        VectorRef {
            data: &self.data[get_lower_index(&index, self.len())..get_upper_index(&index, self.len())]
        }
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>> Indexed<'b, R> for VectorRefMut<'a, T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, index: R) -> Self::Output {
        self.into_const().get(index)
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>> IndexedMut<'b, R> for Vector<T> {
    type Output = VectorRefMut<'b, T>;

    fn get_mut(&'b mut self, index: R) -> Self::Output {
        self.as_mut().get_mut(index)
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>> IndexedMut<'b, R> for VectorRefMut<'a, T> {
    type Output = VectorRefMut<'b, T>;

    fn get_mut(&'b mut self, index: R) -> Self::Output {
        VectorRefMut {
            data: &mut self.data[get_lower_index(&index, self.len())..get_upper_index(&index, self.len())]
        }
    }
}
