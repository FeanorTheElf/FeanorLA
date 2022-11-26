use super::vec::*;
use super::super::prelude::*;

use std::marker::PhantomData;

pub trait VecFn<T>: Sized {

    fn len(&self) -> usize;
    fn at(&self, i: usize) -> T;

    fn add_to<R, V>(self, dst: &mut Vector<V, T>, ring: &R)
        where V: VectorViewMut<T>, R: Ring<El = T>
    {
        assert_eq!(self.len(), dst.len());
        for i in 0..self.len() {
            ring.add_assign(dst.at_mut(i), self.at(i));
        }
    }

    fn assign_to<V>(self, dst: &mut Vector<V, T>)
        where V: VectorViewMut<T>
    {
        assert_eq!(self.len(), dst.len());
        for i in 0..self.len() {
            *dst.at_mut(i) = self.at(i);
        }
    }

    fn compute(self) -> Vector<Vec<T>, T> {
        Vector::from_fn(self.len(), |i| self.at(i))
    }

    fn scaled_ring<R>(self, coeff: El<R>, ring: R) -> VectorScaled<R, Self> 
        where R: Ring<El = T>
    {
        VectorScaled::new(self, coeff, ring)
    }

    fn scaled(self, coeff: T) -> VectorScaled<T::RingType, Self> 
        where T: RingEl
    {
        self.scaled_ring(coeff, T::RING)
    }

    fn neg_ring<R>(self, ring: R) -> VectorNeg<R, Self> 
        where R: Ring<El = T>
    {
        VectorNeg::new(self, ring)
    }

    fn neg(self) -> VectorNeg<T::RingType, Self> 
        where T: RingEl
    {
        self.neg_ring(T::RING)
    }

    fn add_ring<R, V>(self, rhs: V, ring: R) -> VectorSum<R, Self, V>
        where R: Ring<El = T>, V: VecFn<T>
    {
        assert_eq!(self.len(), rhs.len());
        VectorSum::new(self, rhs, ring)
    }

    fn add<V>(self, rhs: V) -> VectorSum<T::RingType, Self, V>
        where T: RingEl, V: VecFn<T>
    {
        self.add_ring(rhs, T::RING)
    }

    fn sub_ring<R, V>(self, rhs: V, ring: R) -> VectorSum<R, Self, VectorNeg<R, V>>
        where R: Ring<El = T>, V: VecFn<T>
    {
        assert_eq!(self.len(), rhs.len());
        self.add_ring(rhs.neg_ring(ring.clone()), ring)
    }

    fn sub<N>(self, rhs: N) -> VectorSum<T::RingType, Self, VectorNeg<T::RingType, N>>
        where T: RingEl, N: VecFn<T>
    {
        self.sub_ring(rhs, T::RING)
    }

    fn map<F, U>(self, func: F) -> VectorMapped<Self, T, U, F>
        where F: Fn(T) -> U
    {
        VectorMapped::new(self, func)
    }
}

impl<V, T> VecFn<T> for Vector<V, T>
    where T: Clone, V: VectorView<T>
{
    fn len(&self) -> usize {
        Vector::<V, T>::len(self)
    }

    fn at(&self, i: usize) -> T {
        Vector::<V, T>::at(self, i).clone()
    }

    default fn add_to<R, W>(self, dst: &mut Vector<W, T>, ring: &R)
        where W: VectorViewMut<T>, R: Ring<El = T>
    {
        assert_eq!(self.len(), dst.len());
        for i in 0..self.len() {
            ring.add_assign_ref(dst.at_mut(i), self.at(i));
        }
    }

    fn compute(self) -> Vector<Vec<T>, T> {
        self.into_owned()
    }
}

pub struct VectorSum<R, V1, V2>
    where R: Ring, V1: VecFn<El<R>>, V2: VecFn<El<R>>
{
    ring: R,
    lhs: V1,
    rhs: V2
}

impl<R, V1, V2> VectorSum<R, V1, V2>
    where R: Ring, V1: VecFn<El<R>>, V2: VecFn<El<R>>
{
    pub fn new(lhs: V1, rhs: V2, ring: R) -> Self {
        assert_eq!(lhs.len(), rhs.len());
        VectorSum { ring, lhs, rhs }
    }
}

impl<R, V1, V2> VecFn<El<R>> for VectorSum<R, V1, V2>
    where R: Ring, V1: VecFn<El<R>>, V2: VecFn<El<R>>
{
    fn len(&self) -> usize {
        self.lhs.len()
    }

    fn at(&self, i: usize) -> El<R> {
        self.ring.add(self.lhs.at(i), self.rhs.at(i))
    }
    
    fn add_to<R2, V>(self, dst: &mut Vector<V, El<R>>, ring: &R2)
        where V: VectorViewMut<El<R>>, R2: Ring<El = El<R>>
    {
        assert_eq!(self.len(), dst.len());
        self.lhs.add_to(dst, ring);
        self.rhs.add_to(dst, ring);
    }

    fn assign_to<V>(self, dst: &mut Vector<V, El<R>>)
        where V: VectorViewMut<El<R>>
    {
        assert_eq!(self.len(), dst.len());
        self.lhs.assign_to(dst);
        self.rhs.add_to(dst, &self.ring);
    }
}

pub struct VectorNeg<R, V> 
    where R: Ring, V: VecFn<El<R>>
{
    ring: R,
    val: V
}

impl<R, V> VectorNeg<R, V>
    where R: Ring, V: VecFn<El<R>>
{
    pub fn new(val: V, ring: R) -> Self {
        VectorNeg { ring, val }
    }
}

impl<R, V> VecFn<El<R>> for VectorNeg<R, V>
    where R: Ring, V: VecFn<El<R>>
{
    fn len(&self) -> usize {
        self.val.len()
    }

    fn at(&self, i: usize) -> El<R> {
        self.ring.neg(self.val.at(i))
    }
}

pub struct VectorScaled<R, V> 
    where R: Ring, V: VecFn<El<R>>
{
    ring: R,
    coeff: El<R>,
    val: V
}

impl<R, V> VectorScaled<R, V>
    where R: Ring, V: VecFn<El<R>>
{
    pub fn new(val: V, coeff: El<R>, ring: R) -> Self {
        VectorScaled { ring, coeff, val }
    }
}

impl<R, M> VecFn<El<R>> for VectorScaled<R, M>
    where R: Ring, M: VecFn<El<R>>
{
    fn len(&self) -> usize {
        self.val.len()
    }

    fn at(&self, i: usize) -> El<R> {
        self.ring.mul_ref(&self.val.at(i), &self.coeff)
    }
}

pub struct VectorKronecker<R, V1, V2>
    where R: Ring, V1: VectorView<El<R>>, V2: VectorView<El<R>>
{
    ring: R,
    lhs: Vector<V1, El<R>>,
    rhs: Vector<V2, El<R>>
}

impl<R, V1, V2> VectorKronecker<R, V1, V2>
    where R: Ring, V1: VectorView<El<R>>, V2: VectorView<El<R>>
{
    pub fn new(lhs: Vector<V1, El<R>>, rhs: Vector<V2, El<R>>, ring: R) -> Self {
        VectorKronecker { ring, lhs, rhs }
    }
}

impl<R, V1, V2> VecFn<El<R>> for VectorKronecker<R, V1, V2>
    where R: Ring, V1: VectorView<El<R>>, V2: VectorView<El<R>>
{
    fn len(&self) -> usize {
        self.lhs.len() * self.rhs.len()
    }

    fn at(&self, i: usize) -> El<R> {
        let lhs_i = i / self.rhs.len();
        let rhs_i = i % self.rhs.len();
        return self.ring.mul_ref(self.lhs.at(lhs_i), self.rhs.at(rhs_i));
    }
}

pub struct VectorMapped<V, T, U, F> 
    where V: VecFn<T>, F: Fn(T) -> U
{
    base: V,
    from: PhantomData<T>,
    to: PhantomData<U>,
    func: F
}

impl<V, T, U, F> VectorMapped<V, T, U, F> 
    where V: VecFn<T>, F: Fn(T) -> U
{
    pub fn new(val: V, func: F) -> Self {
        VectorMapped {
            base: val,
            from: PhantomData,
            to: PhantomData,
            func: func
        }
    }
}

impl<V, T, U, F> VecFn<U> for VectorMapped<V, T, U, F> 
    where V: VecFn<T>, F: Fn(T) -> U
{
    fn len(&self) -> usize {
        self.base.len()
    }

    fn at(&self, i: usize) -> U {
        (self.func)(self.base.at(i))
    }
}

#[test]
fn test_kronecker() {
    let a = Vector::from_array([1, 2, 3]);
    let b = Vector::from_array([1, 1, 0]);
    let expected = Vector::from_array([1, 1, 0, 2, 2, 0, 3, 3, 0]);
    assert_eq!(expected, VectorKronecker::new(a, b, i64::RING).compute());
}