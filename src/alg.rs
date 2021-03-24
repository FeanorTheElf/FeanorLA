use std::ops::{ 
    Add, Mul, Sub, Neg, Div, Rem, AddAssign, MulAssign, 
    DivAssign, SubAssign, RemAssign 
};

use super::algebra::bigint::*;

pub trait Zero: Sized + Add<Self, Output = Self> {
    fn zero() -> Self;
}

impl Zero for i8 {
    fn zero() -> Self {
        0
    }
}

impl Zero for i16 {
    fn zero() -> Self {
        0
    }
}

impl Zero for i32 {
    fn zero() -> Self {
        0
    }
}

impl Zero for i64 {
    fn zero() -> Self {
        0
    }
}

impl Zero for i128 {
    fn zero() -> Self {
        0
    }
}

impl Zero for u8 {
    fn zero() -> Self {
        0
    }
}

impl Zero for u16 {
    fn zero() -> Self {
        0
    }
}

impl Zero for u32 {
    fn zero() -> Self {
        0
    }
}

impl Zero for u64 {
    fn zero() -> Self {
        0
    }
}

impl Zero for u128 {
    fn zero() -> Self {
        0
    }
}

impl Zero for f32 {
    fn zero() -> Self {
        0.
    }
}

impl Zero for f64 {
    fn zero() -> Self {
        0.
    }
}

pub trait One: Sized + Mul<Self, Output = Self> {
    fn one() -> Self;
}

impl One for i8 {
    fn one() -> Self {
        1
    }
}

impl One for i16 {
    fn one() -> Self {
        1
    }
}

impl One for i32 {
    fn one() -> Self {
        1
    }
}

impl One for i64 {
    fn one() -> Self {
        1
    }
}

impl One for i128 {
    fn one() -> Self {
        1
    }
}

impl One for u8 {
    fn one() -> Self {
        1
    }
}

impl One for u16 {
    fn one() -> Self {
        1
    }
}

impl One for u32 {
    fn one() -> Self {
        1
    }
}

impl One for u64 {
    fn one() -> Self {
        1
    }
}

impl One for u128 {
    fn one() -> Self {
        1
    }
}

impl One for f32 {
    fn one() -> Self {
        1.
    }
}

impl One for f64 {
    fn one() -> Self {
        1.
    }
}

pub trait Root {
    fn sqrt(self) -> Self;
}

impl Root for f32 {
    fn sqrt(self) -> f32 {
        self.sqrt()
    }
}

impl Root for f64 {
    fn sqrt(self) -> f64 {
        self.sqrt()
    }
}

///
/// Note that the Ring axioms are mutually exclusive, so a ring providing
/// euclidean division should be RingAxiomsEuclideanRing but not 
/// RingAxiomsIntegralRing even though it is integral
/// 
pub trait RingAxioms {
    fn is_integral() -> bool;
    fn is_euclidean() -> bool;
    fn is_field() -> bool;
}

pub struct RingAxiomsRing;
pub struct RingAxiomsIntegralRing;
pub struct RingAxiomsEuclideanRing;
pub struct RingAxiomsField;

impl RingAxioms for RingAxiomsRing {
    fn is_integral() -> bool { false }
    fn is_euclidean() -> bool { false }
    fn is_field() -> bool { false }
}

impl RingAxioms for RingAxiomsIntegralRing {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { false }
    fn is_field() -> bool { false }
}

impl RingAxioms for RingAxiomsEuclideanRing {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { true }
    fn is_field() -> bool { false }
}

impl RingAxioms for RingAxiomsField {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { false }
    fn is_field() -> bool { true }
}

///
/// Use this for rings whose structure is determined completely at compile time
/// 
/// Multiplication must be commutative
/// 
pub trait RingEl: 
    Clone + Sized + Add<Output = Self> + Mul<Output = Self> + 
    AddAssign + PartialEq + Zero + One + Neg<Output = Self> + 
    Sub<Output = Self> + SubAssign + std::fmt::Debug
{
    type Axioms: RingAxioms;
}

///
/// Division must satisfy the invariant x = (x/y) * y + (x%y)
/// 
pub trait EuclideanRingEl: 
    RingEl<Axioms = RingAxiomsEuclideanRing> + Rem<Output = Self> + 
    RemAssign + Div<Output = Self> + DivAssign 
{
    ///
    /// Computes (returned, self) := (self / rhs, self % rhs) and returns returned.
    /// Can be faster than computing both separately
    /// 
    fn div_rem(&mut self, rhs: Self) -> Self;
}

pub trait FieldEl: 
    RingEl<Axioms = RingAxiomsField> + MulAssign + Div<Output = Self> + DivAssign 
{}

impl RingEl for i8 {
    type Axioms = RingAxiomsEuclideanRing;
}

impl RingEl for i16 {
    type Axioms = RingAxiomsEuclideanRing;
}

impl RingEl for i32 {
    type Axioms = RingAxiomsEuclideanRing;
}

impl RingEl for i64 {
    type Axioms = RingAxiomsEuclideanRing;
}

impl RingEl for i128 {
    type Axioms = RingAxiomsEuclideanRing;
}


impl EuclideanRingEl for i8 {
    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }
}

impl EuclideanRingEl for i16 {
    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }
}

impl EuclideanRingEl for i32 {
    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }
}

impl EuclideanRingEl for i64 {
    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }
}

impl EuclideanRingEl for i128 {
    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }
}

impl RingEl for f32 {
    type Axioms = RingAxiomsField;
}

impl RingEl for f64 {
    type Axioms = RingAxiomsField;
}

impl FieldEl for f32 {}
impl FieldEl for f64 {}

pub trait Integer: Clone + Eq + Ord + EuclideanRingEl {}

impl Integer for i8 {}
impl Integer for i16 {}
impl Integer for i32 {}
impl Integer for i64 {}
impl Integer for i128 {}

pub struct RingElDisplay<'a, R: ?Sized> 
    where R: Ring
{
    ring: &'a R,
    el: &'a R::El
}

impl<'a, R> std::fmt::Display for RingElDisplay<'a, R>
where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.ring.format(self.el, f, false)
    }
}

///
/// Trait to represent a commutative ring with one as a collection of operations 
/// on the elements. More abstract functionality (global properties, like ideals) 
/// is not provided, this is mainly an interface that can be used for algorithms 
/// that deal with ring elements.
/// 
pub trait Ring {
    type El: Sized + Clone + std::fmt::Debug;

    //
    // Design rationale: Types that are cheap to copy can be used with the 
    // add/sub()-functions, for types that are expensive to copy, one can use the 
    // add_ref/sub_ref()-functions. However, as each add/sub() / add_ref/sub_ref() 
    // function returns the result by value, for big types, is is usually 
    // most efficient to clone at least one parameter and potentially reuse 
    // the memory.
    //
    // If an operation can improve efficience by consuming both parameters, one should
    // explicitly implement the default-implemented add/sub()-functions.
    //
    // For multiplication, the situation is usually different: Mostly, multiplication
    // is given by some kind of operation on a cartesian product of the
    // components of lhs and rhs, with a following reduce. In this case, one usually
    // cannot profit from getting the parameters by value. If this is not true, then
    // one should implement the default-implemented mul()-function
    //

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;

    ///
    /// Calculates the product of lhs and rhs. Note that multiplication is assumed to
    /// be commutative.
    /// 
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El;
    fn neg(&self, val: Self::El) -> Self::El;
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool;

    ///
    /// Returns a ring element that can be used as drop-in if one needs some 
    /// unspecified element. This is used only in very exceptional cases.
    /// Best do not touch this - chances are, you will never encounter it.
    /// 
    /// # Example
    /// 
    /// A default implementation of add_assign might move out the value from the
    /// mutable reference, then call add_ref and fill the result in again. However,
    /// if the underlying add-call panics, there is no value to fill in again.
    /// Usually, this is irrelevant, as the variable on which add_assign is called
    /// goes out of scope by the panic, if however the panic is caught (might requrie
    /// unsafe code due to UnwindSafe ???), this is not the case and the value can
    /// be accessed later. In this case, it will be filled with invalid().
    /// 
    fn unspecified_element(&self) -> Self::El {
        self.zero()
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        self.add_ref(self.neg(rhs), lhs)
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.neg(self.add_ref(self.neg(lhs), rhs))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.add_ref(lhs, &rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.mul_ref(&lhs, &rhs) }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.add(lhs, self.neg(rhs))
    }

    fn pow(&self, basis: Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        self.pow_big(basis, BigInt::from(exp as i64))
    }

    fn pow_big(&self, basis: Self::El, exp: BigInt) -> Self::El 
        where Self::El: Clone
    {
        assert!(exp >= 0);
        if exp.is_zero() {
            return self.one();
        }

        let mut power = basis;
        let mut result = self.one();
        for i in 0..(exp.log2_floor() + 1) {
            if exp.is_bit_set(i) {
                result = self.mul_ref(&result, &power);
            }
            power = self.mul(power.clone(), power);
        }
        return result;
    }

    fn is_zero(&self, val: &Self::El) -> bool { self.eq(val, &self.zero()) }
    fn is_one(&self, val: &Self::El) -> bool { self.eq(val, &self.one()) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.eq(val, &self.neg(self.one())) }

    ///
    /// Returns whether the ring is integral, so if there for all nonzero a, b it holds
    /// that ab != 0. It it allowed to return false even for integral rings, as determining
    /// whether the ring is integral might be computationally intractable. 
    /// 
    fn is_integral(&self) -> bool;
    ///
    /// Returns whether the ring is euclidean, so whether the euclidean division and remainder
    /// functions are implemented and behave correctly. It it allowed to return false even for 
    /// euclidean rings, as determining whether the ring is euclidean might be computationally 
    /// intractable, or (more likely) the euclidean division exist mathematically but cannot be 
    /// computed efficiently.
    /// 
    fn is_euclidean(&self) -> bool;
    ///
    /// Returns whether the ring is a field, so whether each nonzero element has a unique 
    /// inverse. It it allowed to return false even for fields, as determining whether the 
    /// ring is a field might be computationally intractable, or (more likely) the inverse
    /// always exists but cannot be computed efficiently.
    /// 
    fn is_field(&self) -> bool;
    
    ///
    /// May panic if the ring is not euclidean (meaning `is_euclidean() == true`). 
    /// The first result is the quotient and the second result the remainder. 
    /// The inequality
    ///  lhs = quo * rhs + rem
    /// must always hold.
    /// 
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El);

    ///
    /// May panic if the ring is not euclidean (meaning `is_euclidean() == true`).
    /// 
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).1 
    }

    ///
    /// May panic if the ring is not euclidean (meaning `is_euclidean() == true`).
    /// 
    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.euclidean_div_rem(lhs, rhs).0
    }

    ///
    /// May panic if the ring is not a field (meaning `is_field() == true`). 
    /// If it does not panic, the result must be valid. For a non-field ring, it therefore 
    /// must panic if rhs does not divide lhs, and if it does, it may either compute the 
    /// correct quotient but may also panic nevertheless.
    /// 
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "{:?}", el)
    }

    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(")?;
        self.format(el, f, false)?;
        write!(f, ")")?;
        return Ok(());
    }
}

pub fn display_ring_el<'a, R>(ring: &'a R, el: &'a R::El) -> RingElDisplay<'a, R> 
    where R: Ring
{
    RingElDisplay {
        ring: ring,
        el: el
    }
}

///
/// Provides a global constant that represents the ring whose elements
/// are the given type. The ring operations are derived from the RingEl
/// implementation on the element type.
/// 
/// Use only for not-expensive-to-copy objects.
/// 
pub struct StaticRingImpl<Axioms, T> 
    where Axioms: RingAxioms, T: RingEl<Axioms = Axioms>
{
    element: std::marker::PhantomData<T>
}

impl<Axioms, T> Copy for StaticRingImpl<Axioms, T> 
    where Axioms: RingAxioms, T: RingEl<Axioms = Axioms>
{}

impl<Axioms, T> Clone for StaticRingImpl<Axioms, T> 
    where Axioms: RingAxioms, T: RingEl<Axioms = Axioms>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<A, T> StaticRingImpl<A, T> 
    where A: RingAxioms, T: RingEl<Axioms = A>
{
    pub const RING: StaticRingImpl<A, T> = StaticRingImpl { 
        element: std::marker::PhantomData 
    };
}

impl<'a, R: Ring> Ring for &'a R {
    type El = R::El;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).add_ref(lhs, rhs) }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { (**self).mul_ref(lhs, rhs) }
    fn neg(&self, val: Self::El) -> Self::El { (**self).neg(val) }
    fn zero(&self) -> Self::El { (**self).zero() }
    fn one(&self) -> Self::El { (**self).one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { (**self).eq(lhs, rhs) }
    fn unspecified_element(&self) -> Self::El { (**self).unspecified_element() }
    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El { (**self).sub_ref_fst(lhs, rhs) }
    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).sub_ref_snd(lhs, rhs) }
    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).add(lhs, rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).mul(lhs, rhs) }
    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).sub(lhs, rhs) }
    fn pow(&self, basis: Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        (**self).pow(basis, exp)
    }
    fn pow_big(&self, basis: Self::El, exp: BigInt) -> Self::El 
        where Self::El: Clone
    {
        (**self).pow_big(basis, exp)
    }
    fn is_zero(&self, val: &Self::El) -> bool { (**self).is_zero(val) }
    fn is_one(&self, val: &Self::El) -> bool { (**self).is_one(val) }
    fn is_neg_one(&self, val: &Self::El) -> bool { (**self).is_neg_one(val) }
    fn is_integral(&self) -> bool { (**self).is_integral() }
    fn is_euclidean(&self) -> bool { (**self).is_euclidean() }
    fn is_field(&self) -> bool { (**self).is_field() }
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { (**self).euclidean_div_rem(lhs, rhs) }
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).euclidean_rem(lhs, rhs) }
    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).euclidean_div(lhs, rhs) }
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).div(lhs, rhs) }
    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result { (**self).format(el, f, in_prod) }
    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result { (**self).format_in_brackets(el, f) }
}

impl<T> Ring for StaticRingImpl<T::Axioms, T>
    where T: RingEl
{
    type El = T;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { lhs + rhs.clone() }
    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs + rhs }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { lhs.clone() * rhs.clone() }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs * rhs }
    fn neg(&self, val: Self::El) -> Self::El { -val }
    fn zero(&self) -> Self::El { T::zero() }
    fn one(&self) -> Self::El { T::one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { lhs == rhs }

    default fn is_integral(&self) -> bool { false }
    default fn is_euclidean(&self) -> bool { false }
    default fn is_field(&self) -> bool { false }

    default fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: &Self::El) -> (Self::El, Self::El) { 
        panic!("Not a euclidean domain!");
    }

    default fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El { 
        panic!("Not a field!");
    }

    default fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).0
    }
    
    default fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).1
    }
}

impl<T> Ring for StaticRingImpl<RingAxiomsIntegralRing, T> 
    where T: RingEl<Axioms = RingAxiomsIntegralRing>
{
    fn is_integral(&self) -> bool { true }
}

impl<T> Ring for StaticRingImpl<RingAxiomsEuclideanRing, T> 
    where T: EuclideanRingEl
{
    fn is_integral(&self) -> bool { true }
    fn is_euclidean(&self) -> bool { true }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { 
        let quo = lhs.div_rem(rhs.clone());
        return (quo, lhs);
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        lhs / rhs.clone()
    }
    
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        lhs % rhs.clone()
    }

    default fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El { 
        panic!("Not a field!");
    }
}

impl<T> Ring for StaticRingImpl<RingAxiomsEuclideanRing, T> 
    where T: Integer
{
    fn div(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        let result = lhs.div_rem(rhs.clone());
        if !self.is_zero(&lhs) {
            panic!("The integers are not a field and tried division by {:?} has remainder {:?}!", rhs, lhs);
        }
        return result;
    }
}

impl<T> Ring for StaticRingImpl<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_integral(&self) -> bool { true }
    fn is_field(&self) -> bool { true }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        lhs / rhs.clone()
    }
}


pub type StaticRing<R> = StaticRingImpl<<R as RingEl>::Axioms, R>;


#[test]
fn test_pow() {
    assert_eq!(81 * 81 * 3, StaticRing::<i64>::RING.pow(3, 9));
}

macro_rules! impl_euclidean_ring_el {
    ($t:ty: $ring_constant:expr) => {
        impl std::ops::Add for $t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add(self, rhs)
            }
        }

        impl std::ops::Add<&$t> for $t {
            type Output = $t;

            fn add(self, rhs: &$t) -> Self::Output {
                ($ring_constant).add_ref(self, rhs)
            }
        }

        impl std::ops::Add<$t> for &$t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add_ref(rhs, self)
            }
        }

        impl std::ops::AddAssign for $t {

            fn add_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add(v, rhs));
            }
        }

        impl std::ops::AddAssign<&$t> for $t {

            fn add_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add_ref(v, rhs));
            }
        }

        impl std::ops::Mul for $t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul(self, rhs)
            }
        }

        impl std::ops::Mul<&$t> for $t {
            type Output = $t;

            fn mul(self, rhs: &$t) -> Self::Output {
                ($ring_constant).mul_ref(&self, rhs)
            }
        }

        impl std::ops::Mul<$t> for &$t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul_ref(&rhs, self)
            }
        }

        impl std::ops::MulAssign for $t {

            fn mul_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).mul(v, rhs));
            }
        }

        impl std::ops::MulAssign<&$t> for $t {

            fn mul_assign(&mut self, rhs: &$t) {
                *self = ($ring_constant).mul_ref(self, rhs);
            }
        }

        impl std::ops::Sub for $t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub(self, rhs)
            }
        }

        impl std::ops::Sub<&$t> for $t {
            type Output = $t;

            fn sub(self, rhs: &$t) -> Self::Output {
                ($ring_constant).sub_ref_snd(self, rhs)
            }
        }

        impl std::ops::Sub<$t> for &$t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub_ref_fst(self, rhs)
            }
        }

        impl std::ops::SubAssign for $t {

            fn sub_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub(v, rhs));
            }
        }

        impl std::ops::SubAssign<&$t> for $t {

            fn sub_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub_ref_snd(v, rhs));
            }
        }

        impl std::ops::Div for $t {
            type Output = $t;

            fn div(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_div(self, &rhs)
            }
        }

        impl std::ops::Rem for $t {
            type Output = $t;

            fn rem(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_rem(self, &rhs)
            }
        }

        impl std::ops::RemAssign for $t {
            fn rem_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_rem(v, &rhs));
            }
        }

        impl std::ops::DivAssign for $t {
            fn div_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_div(v, &rhs));
            }
        }

        impl PartialEq for $t {
            fn eq(&self, rhs: &$t) -> bool {
                ($ring_constant).eq(self, rhs)
            }
        }

        impl std::ops::Neg for $t {
            type Output = $t;

            fn neg(self) -> Self::Output {
                ($ring_constant).neg(self)
            }
        }

        impl Zero for $t {
            fn zero() -> $t {
                ($ring_constant).zero()
            }
        }

        impl One for $t {
            fn one() -> $t {
                ($ring_constant).one()
            }
        }

        impl std::fmt::Display for $t {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                ($ring_constant).format(self, f, false)
            }
        }

        impl RingEl for $t {
            type Axioms = RingAxiomsEuclideanRing;
        }

        impl EuclideanRingEl for $t {

            fn div_rem(&mut self, rhs: $t) -> $t {
                let mut result: Option<$t> = None;
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| {
                    let (quo, rem) = ($ring_constant).euclidean_div_rem(v, &rhs);
                    result = Some(quo);
                    return rem;
                });
                return result.unwrap();
            }
        }
    };
}