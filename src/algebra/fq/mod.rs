use super::super::ring::*;
use super::super::bigint::*;
use super::super::wrapper::*;

pub trait FiniteRingIterFn<R: FiniteRing>: Clone {

    fn next(&mut self, ring: &R) -> Option<R::El>;
}

impl<'a, R: FiniteRing, F: FiniteRingIterFn<R>> FiniteRingIterFn<&'a R> for F {
    
    fn next(&mut self, ring: &&'a R) -> Option<R::El> {
        <F as FiniteRingIterFn<R>>::next(self, *ring)
    }
}

pub trait FiniteRing : Ring {

    type IterFn: FiniteRingIterFn<Self>;

    fn size(&self) -> BigInt;
    fn iter_fn(&self) -> Self::IterFn;
}

impl<'a, R> FiniteRing for &'a R
    where R: FiniteRing
{
    type IterFn = R::IterFn;

    fn size(&self) -> BigInt { (**self).size() }
    fn iter_fn(&self) -> Self::IterFn { (**self).iter_fn() }
}

pub struct WrappingRingIterFn<R: FiniteRing> {
    base_fn: R::IterFn
}

impl<R: FiniteRing> Clone for WrappingRingIterFn<R> {

    fn clone(&self) -> Self {
        WrappingRingIterFn {
            base_fn: self.base_fn.clone()
        }
    }
}

impl<R: FiniteRing> FiniteRingIterFn<WrappingRing<R>> for WrappingRingIterFn<R> {

    fn next(&mut self, ring: &WrappingRing<R>) -> Option<RingElWrapper<R>> {
        let el = self.base_fn.next(ring.wrapped_ring());
        el.map(|e| ring.wrapped_ring().bind_by_value(e))
    }
}

impl<R> FiniteRing for WrappingRing<R>
    where R: FiniteRing
{
    type IterFn = WrappingRingIterFn<R>;

    fn size(&self) -> BigInt {
        self.wrapped_ring().size()
    }
    
    fn iter_fn(&self) -> Self::IterFn {
        WrappingRingIterFn {
            base_fn: self.wrapped_ring().iter_fn()
        }
    }
}

#[derive(Clone)]
pub struct FiniteRingElementIter<R>
    where R: FiniteRing
{
    ring: R,
    iter_fn: R::IterFn
}

impl<R> Iterator for FiniteRingElementIter<R>
    where R: FiniteRing
{
    type Item = R::El;

    fn next(&mut self) -> Option<R::El> {
        self.iter_fn.next(&self.ring)
    }
}

pub fn elements<R: FiniteRing>(ring: R) -> FiniteRingElementIter<R> {
    FiniteRingElementIter {
        iter_fn: ring.iter_fn(),
        ring: ring
    }
}

pub mod zn_big;
pub mod zn_small;
pub mod fq_small;