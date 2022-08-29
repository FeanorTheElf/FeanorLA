use super::ring::*;
use super::integer::*;
use super::wrapper::*;

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
    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32;
}

impl<'a, R> FiniteRing for &'a R
    where R: FiniteRing
{
    type IterFn = R::IterFn;

    fn size(&self) -> BigInt { (**self).size() }
    fn iter_fn(&self) -> Self::IterFn { (**self).iter_fn() }

    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32 {
            (**self).random_element(rng)
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

pub fn finite_field_elements<R: FiniteRing>(ring: R) -> FiniteRingElementIter<R> {
    FiniteRingElementIter {
        iter_fn: ring.iter_fn(),
        ring: ring
    }
}

pub mod zn_big;
pub mod zn_small;
pub mod fq_small;

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
        el.map(|e| ring.from(e))
    }
}

impl<R> FiniteRing for WrappingRing<R>
    where R: FiniteRing
{
    type IterFn = WrappingRingIterFn<R>;

    fn random_element<G>(&self, rng: G) -> RingElWrapper<R>
        where G: FnMut() -> u32
    {
        self.from(self.wrapped_ring().random_element(rng))
    }

    fn size(&self) -> BigInt {
        self.wrapped_ring().size()
    }

    fn iter_fn(&self) -> Self::IterFn {
        WrappingRingIterFn {
            base_fn: self.wrapped_ring().iter_fn()
        }
    }
}

impl<R> WrappingRing<R>
    where R: FiniteRing
{
    pub fn elements(&self) -> FiniteRingElementIter<&WrappingRing<R>> {
        finite_field_elements(&self)
    }
}