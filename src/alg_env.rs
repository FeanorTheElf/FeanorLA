use super::alg::*;
use super::algebra::bigint::BigInt;
use super::alg_macros::*;

use std::marker::PhantomData;

#[derive(Debug)]
pub struct RingReferencingEl<'a, R: ?Sized, A: RingAxioms>
    where R: Ring
{
    ring: &'a R,
    axioms: PhantomData<A>,
    el: R::El
}

impl<'a, R, A> Copy for RingReferencingEl<'a, R, A>
where R: Ring, R::El: Copy, A: RingAxioms
{}

impl<'a, R, A> Clone for RingReferencingEl<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn clone(&self) -> Self {
        RingReferencingEl {
            ring: self.ring,
            axioms: PhantomData,
            el: self.el.clone()
        }
    }
}

pub trait BindableElementRing: Ring {

    fn bind<A: RingAxioms>(&self, el: Self::El) -> RingReferencingEl<Self, A>;
}

impl<R: Ring> BindableElementRing for R {
    
    fn bind<A: RingAxioms>(&self, el: Self::El) -> RingReferencingEl<Self, A> {
        assert!(self.is_euclidean() == A::is_euclidean());
        assert!(self.is_field() == A::is_field());
        assert!(self.is_integral() == A::is_integral());
        RingReferencingEl::create(self, el)
    }
}

impl<'a, R, A> RingReferencingEl<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn create(ring: &'a R, val: R::El) -> RingReferencingEl<'a, R, A> {
        return RingReferencingEl {
            ring: ring,
            el: val,
            axioms: PhantomData
        };
    }

    pub fn val(&self) -> &R::El {
        &self.el
    }

    pub fn val_mut(&mut self) -> &mut R::El {
        &mut self.el
    }
}

struct WrappingRing<'a, R, A> 
{
    ring: PhantomData<&'a R>,
    axioms: PhantomData<A>
}

impl<'a, R, A> WrappingRing<'a, R, A> 
where R: 'a + Ring + std::fmt::Debug, A: RingAxioms + std::fmt::Debug 
{
    const RING: WrappingRing<'a, R, A> = WrappingRing {
        ring: PhantomData,
        axioms: PhantomData
    };
}

impl<'a, R, A> Ring for WrappingRing<'a, R, A> 
where R: 'a + Ring + std::fmt::Debug, A: RingAxioms + std::fmt::Debug
{
    type El = RingReferencingEl<'a, R, A>;
    
    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.add_ref(lhs.el, rhs.val()))
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.mul_ref(lhs.val(), rhs.val()))
    }

    fn neg(&self, val: Self::El) -> Self::El  {
        RingReferencingEl::create(val.ring, val.ring.neg(val.el))
    }

    fn zero(&self) -> Self::El { 
        panic!("`WrappingRing` does not provide the constants from the ring; Use `base_ring.bind(base_ring.zero())` instead.") 
    }

    fn one(&self) -> Self::El { 
        panic!("`WrappingRing` does not provide the constants from the ring; Use `base_ring.bind(base_ring.one())` instead.")
    }

    fn unspecified_element(&self) -> Self::El {
        panic!("`WrappingRing` does not provide the constants from the ring; Use `base_ring.bind(base_ring.unspecified_element())` instead.")
    }

    fn from_z(&self, _x: i64) -> Self::El { 
        panic!("`WrappingRing` does not provide the constants from the ring; Use `base_ring.bind(base_ring.unspecified_element())` instead.")
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        lhs.ring.eq(lhs.val(), rhs.val())
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_fst(lhs.val(), rhs.el))
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_snd(lhs.el, rhs.val()))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        RingReferencingEl::create(lhs.ring, lhs.ring.add(lhs.el, rhs.el))
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.mul(lhs.el, rhs.el))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.sub(lhs.el, rhs.el))
    }

    fn pow(&self, basis: Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        RingReferencingEl::create(basis.ring, basis.ring.pow(basis.el, exp))
    }

    fn pow_big(&self, basis: Self::El, exp: BigInt) -> Self::El 
        where Self::El: Clone
    {
        RingReferencingEl::create(basis.ring, basis.ring.pow_big(basis.el, exp))
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        val.ring.is_zero(val.val())
    }

    fn is_one(&self, val: &Self::El) -> bool {
        val.ring.is_one(val.val())
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        val.ring.is_neg_one(val.val())
    }

    fn is_integral(&self) -> bool {
        A::is_integral()
    }

    fn is_euclidean(&self) -> bool {
        A::is_euclidean()
    }

    fn is_field(&self) -> bool {
        A::is_field()
    }
    
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        let (quo, rem) = lhs.ring.euclidean_div_rem(lhs.el, rhs.val());
        (
            RingReferencingEl::create(lhs.ring, quo),
            RingReferencingEl::create(lhs.ring, rem)
        )
    }

    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_rem(lhs.el, rhs.val()))
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_div(lhs.el, rhs.val()))
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingReferencingEl::create(lhs.ring, lhs.ring.div(lhs.el, rhs.val()))
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        el.ring.format(el.val(), f, in_prod)
    }

    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        el.ring.format_in_brackets(el.val(), f)
    }
}

impl_ring_el!{ 
    RingReferencingEl<'a, R, A>; 
    WrappingRing::<'a, R, A>::RING; 
    A; 
    'a, R: 'a + Ring + std::fmt::Debug, A: RingAxioms + std::fmt::Debug
}
impl_euclidean_ring_el!{ 
    RingReferencingEl<'a, R, RingAxiomsEuclideanRing>; 
    WrappingRing::<'a, R, RingAxiomsEuclideanRing>::RING;
    'a, R: 'a + Ring + std::fmt::Debug 
}
impl_field_ring_el!{ 
    RingReferencingEl<'a, R, RingAxiomsField>; 
    WrappingRing::<'a, R, RingAxiomsField>::RING; 
    'a, R: 'a + Ring + std::fmt::Debug 
}

#[test]
fn test_exper() {
    let ring = StaticRing::<i64>::RING;
    let a = ring.bind::<RingAxiomsEuclideanRing>(4);
    let b = ring.bind(10);

    let mut c = a + (b * a) + a + a;
    c = c / b;
    c = c * a;
    let d = c;

    assert_eq!(ring.bind(20), d);
}