use super::prelude::*;
use super::wrapper::*;

pub trait FiniteRingIterFn<R: FiniteRing>: Clone {

    fn next(&mut self, ring: &R) -> Option<R::El>;
}

impl<R, F> FiniteRingIterFn<R> for F 
    where R: RingDecorator, R::DecoratedRing: FiniteRing, F: FiniteRingIterFn<R::DecoratedRing>
{
    
    fn next(&mut self, ring: &R) -> Option<El<R>> {
        <F as FiniteRingIterFn<R::DecoratedRing>>::next(self, ring.decorated_ring())
    }
}

pub trait FiniteRing: Ring {

    type IterFn: FiniteRingIterFn<Self>;

    fn size(&self) -> BigInt;
    fn iter_fn(&self) -> Self::IterFn;
    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32;
}

pub trait IntegerQuotientRing: FiniteRing {

    type LiftingRing: IntegerRing;

    fn lifting_ring(&self) -> Self::LiftingRing;
    fn lift(&self, x: &El<Self>, ring: &Self::LiftingRing) -> El<Self::LiftingRing>;
}

impl<'a, R> FiniteRing for R
    where R: RingDecorator, R::DecoratedRing: FiniteRing
{
    type IterFn = <R::DecoratedRing as FiniteRing>::IterFn;

    fn size(&self) -> BigInt { self.decorated_ring().size() }
    fn iter_fn(&self) -> Self::IterFn { self.decorated_ring().iter_fn() }

    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32 
    {
        self.decorated_ring().random_element(rng)
    }
}

impl<R> IntegerQuotientRing for R
    where R: RingDecorator, R::DecoratedRing: IntegerQuotientRing
{
    type LiftingRing = <R::DecoratedRing as IntegerQuotientRing>::LiftingRing;

    fn lifting_ring(&self) -> Self::LiftingRing { self.decorated_ring().lifting_ring() }
    fn lift(&self, x: &El<Self>, ring: &Self::LiftingRing) -> El<Self::LiftingRing> { self.decorated_ring().lift(x, ring) }
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