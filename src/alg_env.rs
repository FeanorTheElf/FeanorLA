use super::alg::*;
use super::algebra::bigint::BigInt;
use super::alg_macros::*;

use std::marker::PhantomData;

pub struct RingReferencingElBase<'a, R: ?Sized, A: RingAxioms>
    where R: Ring
{
    ring: &'a R,
    axioms: PhantomData<A>,
    el: R::El
}

impl<'a, R: ?Sized, A: RingAxioms> std::fmt::Debug for RingReferencingElBase<'a, R, A> 
    where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "El {} referencing ring of type {}", self.ring.display(&self.el), std::any::type_name::<R>())
    }
}

impl<'a, R, A> Copy for RingReferencingElBase<'a, R, A>
where R: Ring, R::El: Copy, A: RingAxioms
{}

impl<'a, R, A> Clone for RingReferencingElBase<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn clone(&self) -> Self {
        RingReferencingElBase {
            ring: self.ring,
            axioms: PhantomData,
            el: self.el.clone()
        }
    }
}

pub enum RingReferencingEl<'a, R: ?Sized, A: RingAxioms>
    where R: Ring
{
    El(RingReferencingElBase<'a, R, A>),
    Integer(i64)
}

impl<'a, R: ?Sized, A: RingAxioms> std::fmt::Debug for RingReferencingEl<'a, R, A> 
    where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RingReferencingEl::El(el) => write!(f, "{:?}", el),
            RingReferencingEl::Integer(x) => write!(f, "Integer {} as RingReferencingEl", *x)
        }
    }
}

impl<'a, R, A> Copy for RingReferencingEl<'a, R, A>
where R: Ring, R::El: Copy, A: RingAxioms
{}

impl<'a, R, A> Clone for RingReferencingEl<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn clone(&self) -> Self {
        match self {
            RingReferencingEl::El(el) => RingReferencingEl::El(el.clone()),
            RingReferencingEl::Integer(x) => RingReferencingEl::Integer(*x)
        }
    }
}

pub trait BindableElementRing: Ring {

    fn bind<A: RingAxioms>(&self, el: Self::El) -> RingReferencingEl<Self, A>;
}

impl<R: Ring> BindableElementRing for R {
    
    fn bind<A: RingAxioms>(&self, el: Self::El) -> RingReferencingEl<Self, A> {
        debug_assert!(self.is_euclidean() == A::is_euclidean());
        debug_assert!(self.is_field() == A::is_field());
        debug_assert!(self.is_integral() == A::is_integral());
        RingReferencingEl::create(self, el)
    }
}

impl<'a, R, A> RingReferencingEl<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn create(ring: &'a R, val: R::El) -> RingReferencingEl<'a, R, A> {
        return RingReferencingEl::El(RingReferencingElBase {
            ring: ring,
            el: val,
            axioms: PhantomData
        });
    }

    ///
    /// Returns the underlying ring element. For some integer constants,
    /// the underlying ring might not be known, in this case the integer
    /// value is returned using the Error case.
    /// 
    pub fn val(&self) -> Result<&R::El, i64> {
        match self {
            RingReferencingEl::El(el) => Ok(&el.el),
            RingReferencingEl::Integer(x) => Err(*x)
        }
    }

    pub fn unwrap(&self) -> &R::El {
        self.val().unwrap()
    }

    ///
    /// Returns a mutable reference to the underlying ring element. 
    /// For some integer constants, the underlying ring might not be 
    /// known, in this case the integer value is returned using the
    /// Error case.
    /// 
    /// Mutating the returned element will mutate this object element.
    /// 
    pub fn val_mut(&mut self) -> Result<&mut R::El, &mut i64> {
        match self {
            RingReferencingEl::El(el) => Ok(&mut el.el),
            RingReferencingEl::Integer(x) => Err(&mut *x)
        }
    }

    fn val_ring(self, ring: &R) -> R::El {
        match self {
            RingReferencingEl::El(el) => el.el.clone(),
            RingReferencingEl::Integer(x) => ring.from_z(x)
        }
    }

    fn ring(&self) -> Option<&'a R> {
        match self {
            RingReferencingEl::El(el) => Some(el.ring),
            RingReferencingEl::Integer(_) => None
        }
    }

    fn integer(&self) -> Option<i64> {
        match self {
            RingReferencingEl::El(_) => None,
            RingReferencingEl::Integer(x) => Some(*x)
        }
    }
}

impl<'a, R, A> From<i8> for RingReferencingEl<'a, R, A>
where R: Ring, A: RingAxioms
{
    fn from(_: i8) -> Self {
        panic!("`RingReferencingEl` does not provide the constants from the ring; Use `base_ring.bind(base_ring.from_z())` instead.")
    }
}

pub struct WrappingRing<'a, R, A> 
{
    ring: PhantomData<&'a R>,
    axioms: PhantomData<A>
}

impl<'a, R, A> WrappingRing<'a, R, A> 
where R: 'a + Ring, A: RingAxioms 
{
    const RING: WrappingRing<'a, R, A> = WrappingRing {
        ring: PhantomData,
        axioms: PhantomData
    };
}

impl<'a, R, A> std::fmt::Debug for WrappingRing<'a, R, A> {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "WrappingRing for {}",  std::any::type_name::<R>())
    }
}

impl<'a, R, A> Clone for WrappingRing<'a, R, A> {

    fn clone(&self) -> Self {
        WrappingRing {
            ring: PhantomData,
            axioms: PhantomData
        }
    }
}

fn panic_equality_impossible() -> ! {
    panic!("Two RingReferencingEls are being compared, but both are just constant integers and contain no data about the underlying ring. Since we dont know whether the embedding Z -> R is injective, we cannot compare them. This is a fundamental issue with the RingReferencingEl-approach, as on the one hand each element must contain a reference to the underlying ring, but on the other hand, we want to provide integer constants without giving access to such a reference. When creating the constants, try using `base_ring.bind(base_ring.from_z(...))` instead of `WrappingRing.from_z(...)` which should fix the issue.")
}

impl<'a, R, A> Ring for WrappingRing<'a, R, A> 
where R: 'a + Ring, A: RingAxioms
{
    type El = RingReferencingEl<'a, R, A>;
    
    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.add_ref(lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.add_ref(lhs.ring.from_z(*x), &lhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.add_ref(rhs.ring.from_z(x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(x + *y)
        }
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.mul_ref(&lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.mul_ref(&lhs.ring.from_z(*x), &lhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.mul_ref(&rhs.ring.from_z(*x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(*x * *y)
        }
    }

    fn neg(&self, val: Self::El) -> Self::El  {
        match val {
            RingReferencingEl::El(val) => RingReferencingEl::create(val.ring, val.ring.neg(val.el)),
            RingReferencingEl::Integer(x) => RingReferencingEl::Integer(-x)
        }
    }

    fn zero(&self) -> Self::El {
        RingReferencingEl::Integer(0)
    }

    fn one(&self) -> Self::El {
        RingReferencingEl::Integer(1)
    }

    fn from_z(&self, x: i64) -> Self::El { 
        RingReferencingEl::Integer(x)
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => lhs.ring.eq(&lhs.el, &rhs.el),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => lhs.ring.eq(&lhs.ring.from_z(*x), &lhs.el),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => rhs.ring.eq(&rhs.ring.from_z(*x), &rhs.el),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y)) if x == y
                 => true,
            (RingReferencingEl::Integer(_), RingReferencingEl::Integer(_))
                 => panic_equality_impossible()
        }
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_fst(&lhs.el, rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_fst(&lhs.el, lhs.ring.from_z(x))),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.sub_ref_fst(&rhs.ring.from_z(*x), rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(*x - y)
        }
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_snd(lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.sub_ref_snd(lhs.el, &lhs.ring.from_z(*x))),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.sub_ref_snd(rhs.ring.from_z(x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(x - *y)
        }
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        if let Some(ring) = lhs.ring() {
            RingReferencingEl::create(ring, ring.add(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else if let Some(ring) = rhs.ring() {
            RingReferencingEl::create(ring, ring.add(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else {
            RingReferencingEl::Integer(lhs.integer().unwrap() + rhs.integer().unwrap())
        }
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        if let Some(ring) = lhs.ring() {
            RingReferencingEl::create(ring, ring.mul(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else if let Some(ring) = rhs.ring() {
            RingReferencingEl::create(ring, ring.mul(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else {
            RingReferencingEl::Integer(lhs.integer().unwrap() * rhs.integer().unwrap())
        }
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        if let Some(ring) = lhs.ring() {
            RingReferencingEl::create(ring, ring.sub(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else if let Some(ring) = rhs.ring() {
            RingReferencingEl::create(ring, ring.sub(lhs.val_ring(ring), rhs.val_ring(ring)))
        } else {
            RingReferencingEl::Integer(lhs.integer().unwrap() - rhs.integer().unwrap())
        }
    }

    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        match basis {
            RingReferencingEl::El(el) => RingReferencingEl::create(&el.ring, el.ring.pow(&el.el, exp)),
            RingReferencingEl::Integer(x) => RingReferencingEl::Integer(x.pow(exp))
        }
    }

    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El 
        where Self::El: Clone
    {
        match basis {
            RingReferencingEl::El(el) => RingReferencingEl::create(&el.ring, el.ring.pow_big(&el.el, exp)),
            RingReferencingEl::Integer(x) => RingReferencingEl::Integer(i64::RING.pow_big(x, exp))
        }
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        match val {
            RingReferencingEl::El(val) => val.ring.is_zero(&val.el),
            RingReferencingEl::Integer(0) => true,
            RingReferencingEl::Integer(_) => panic_equality_impossible(),
        }
    }

    fn is_one(&self, val: &Self::El) -> bool {
        match val {
            RingReferencingEl::El(val) => val.ring.is_zero(&val.el),
            RingReferencingEl::Integer(1) => true,
            RingReferencingEl::Integer(_) => panic_equality_impossible(),
        }
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        match val {
            RingReferencingEl::El(val) => val.ring.is_zero(&val.el),
            RingReferencingEl::Integer(-1) => true,
            RingReferencingEl::Integer(_) => panic_equality_impossible(),
        }
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
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs)) => {
                let (quo, rem) = lhs.ring.euclidean_div_rem(lhs.el, &rhs.el);
                (
                    RingReferencingEl::create(&lhs.ring, quo),
                    RingReferencingEl::create(&lhs.ring, rem)
                )
            },
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x)) => {
                let (quo, rem) = lhs.ring.euclidean_div_rem(lhs.el, &lhs.ring.from_z(*x));
                (
                    RingReferencingEl::create(&lhs.ring, quo),
                    RingReferencingEl::create(&lhs.ring, rem)
                )
            },
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs)) => {
                let (quo, rem) = rhs.ring.euclidean_div_rem(rhs.ring.from_z(x), &rhs.el);
                (
                    RingReferencingEl::create(&rhs.ring, quo),
                    RingReferencingEl::create(&rhs.ring, rem)
                )
            },
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y)) => (
                RingReferencingEl::Integer(x / *y),
                RingReferencingEl::Integer(x % *y)
            )
        }
    }

    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_rem(lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_rem(lhs.el, &lhs.ring.from_z(*x))),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.euclidean_rem(rhs.ring.from_z(x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(x % *y)
        }
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_div(lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.euclidean_div(lhs.el, &lhs.ring.from_z(*x))),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.euclidean_div(rhs.ring.from_z(x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(x / *y)
        }
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        match (lhs, rhs) {
            (RingReferencingEl::El(lhs), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.div(lhs.el, &rhs.el)),
            (RingReferencingEl::El(lhs), RingReferencingEl::Integer(x))
                 => RingReferencingEl::create(lhs.ring, lhs.ring.div(lhs.el, &lhs.ring.from_z(*x))),
            (RingReferencingEl::Integer(x), RingReferencingEl::El(rhs))
                 => RingReferencingEl::create(rhs.ring, rhs.ring.div(rhs.ring.from_z(x), &rhs.el)),
            (RingReferencingEl::Integer(x), RingReferencingEl::Integer(y))
                 => RingReferencingEl::Integer(x / *y)
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        match el {
            RingReferencingEl::El(el) => el.ring.format(&el.el, f, in_prod),
            RingReferencingEl::Integer(x) => write!(f, "{}", x)
        }
    }

    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match el {
            RingReferencingEl::El(el) => el.ring.format_in_brackets(&el.el, f),
            RingReferencingEl::Integer(x) => write!(f, "({})", x)
        }
    }
}

impl_ring_el!{ 
    RingReferencingEl<'a, R, A>; 
    WrappingRing::<'a, R, A>::RING; WrappingRing<'a, R, A>;
    A; 
    'a, R: 'a + Ring, A: RingAxioms
}
impl_euclidean_ring_el!{ 
    RingReferencingEl<'a, R, RingAxiomsEuclideanRing>; 
    WrappingRing::<'a, R, RingAxiomsEuclideanRing>::RING;
    'a, R: 'a + Ring 
}
impl_field_ring_el!{ 
    RingReferencingEl<'a, R, RingAxiomsField>; 
    WrappingRing::<'a, R, RingAxiomsField>::RING; 
    'a, R: 'a + Ring 
}

impl<'a, R: ?Sized, A: RingAxioms> std::ops::Add<i64> for RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn add(self, rhs: i64) -> Self::Output {
        self + WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, 'b, R: ?Sized, A: RingAxioms> std::ops::Add<i64> for &'b RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn add(self, rhs: i64) -> Self::Output {
        self + WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, R: ?Sized, A: RingAxioms> std::ops::Mul<i64> for RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn mul(self, rhs: i64) -> Self::Output {
        self * WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, 'b, R: ?Sized, A: RingAxioms> std::ops::Mul<i64> for &'b RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn mul(self, rhs: i64) -> Self::Output {
        self * WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, R: ?Sized, A: RingAxioms> std::ops::Sub<i64> for RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn sub(self, rhs: i64) -> Self::Output {
        self - WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, 'b, R: ?Sized, A: RingAxioms> std::ops::Sub<i64> for &'b RingReferencingEl<'a, R, A>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, A>;
    
    fn sub(self, rhs: i64) -> Self::Output {
        self - WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, R: ?Sized> std::ops::Div<i64> for RingReferencingEl<'a, R, RingAxiomsField>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, RingAxiomsField>;
    
    fn div(self, rhs: i64) -> Self::Output {
        self / WrappingRing::RING.from_z(rhs)
    }
}

impl<'a, 'b, R: ?Sized> std::ops::Div<i64> for &'b RingReferencingEl<'a, R, RingAxiomsField>
    where R: Ring
{
    type Output = RingReferencingEl<'a, R, RingAxiomsField>;
    
    fn div(self, rhs: i64) -> Self::Output {
        self.clone() / WrappingRing::RING.from_z(rhs)
    }
}

#[test]
fn test_ring_referencing_el_operations() {
    let ring = StaticRing::<i64>::RING;
    let a = ring.bind::<RingAxiomsEuclideanRing>(4);
    let b = ring.bind(10);

    let mut c = a + (b * a) + a + a;
    c = c / b;
    c = c * a;
    let d = c;

    assert_eq!(ring.bind(20), d);
}

#[test]
#[should_panic]
fn test_throw_on_wrong_axioms() {
    let ring = StaticRing::<i64>::RING;
    let a = ring.bind::<RingAxiomsField>(4);
    let b = ring.bind(10);
    let _ = a / b;
}