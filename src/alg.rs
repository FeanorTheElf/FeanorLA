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
    Sub<Output = Self> + SubAssign 
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

///
/// For types that may have rounding errors
/// 
pub trait Float: Clone + PartialEq + PartialOrd + FieldEl {}

impl Float for f32 {}
impl Float for f64 {}

pub trait Ring {
    type El: Sized + Clone;

    //
    // Design rationale: Types that are cheap to copy can be used with the 
    // op()-functions, for types that are expensive to copy, one can use the 
    // op_ref()-functions. However, as each op() / op_ref() function returns 
    // the result by value, for big types, is is usually most efficient to 
    // clone at least one parameter and potentially reuse the memory.
    //
    // If an operation can improve efficience by consuming both parameters, one should
    // explicitly implement the default-implemented op()-functions.
    //

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;
    fn mul_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;
    fn neg(&self, val: Self::El) -> Self::El;
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool;

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        self.add_ref(self.neg(rhs), lhs)
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.neg(self.add_ref(self.neg(lhs), rhs))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.add_ref(lhs, &rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.mul_ref(lhs, &rhs) }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.add(lhs, self.neg(rhs))
    }

    fn pow(&self, basis: Self::El, exp: u64) -> Self::El 
        where Self::El: Clone
    {
        self.pow_big(basis, BigInt::from(exp))
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
                result = self.mul_ref(result, &power);
            }
            power = self.mul(power.clone(), power);
        }
        return result;
    }

    fn is_zero(&self, val: &Self::El) -> bool { self.eq(val, &self.zero()) }
    fn is_one(&self, val: &Self::El) -> bool { self.eq(val, &self.one()) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.eq(val, &self.neg(self.one())) }

    fn is_integral(&self) -> bool;
    fn is_euclidean(&self) -> bool;
    fn is_field(&self) -> bool;
    
    ///
    /// May panic if the ring is not euclidean. The first result is the quotient
    /// and the second result the remainder. The inequality
    ///  lhs = quo * rhs + rem
    /// must always hold.
    /// 
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: Self::El) -> (Self::El, Self::El);

    ///
    /// May panic if the ring is not euclidean
    /// 
    fn euclidean_rem(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).1 
    }

    ///
    /// May panic if the ring is not euclidean
    /// 
    fn euclidean_div(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.euclidean_div_rem(lhs, rhs).0
    }

    ///
    /// May panic if the ring is not a field. If it does not panic, the result
    /// must be valid. For a non-field ring, it therefore must panic if rhs does not
    /// divide lhs, and if it does, it may either compute the correct quotient but may
    /// also panic nevertheless.
    /// 
    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El;
}

pub struct StaticRing<Axioms, T> 
    where Axioms: RingAxioms, T: RingEl<Axioms = Axioms>
{
    element: std::marker::PhantomData<T>
}

impl<A, T> StaticRing<A, T> 
    where A: RingAxioms, T: RingEl<Axioms = A>
{
    pub const RING: StaticRing<A, T> = StaticRing { 
        element: std::marker::PhantomData 
    };
}

impl<T> Ring for StaticRing<T::Axioms, T>
    where T: RingEl
{
    type El = T;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { lhs + rhs.clone() }
    fn mul_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { lhs * rhs.clone() }
    fn neg(&self, val: Self::El) -> Self::El { -val }
    fn zero(&self) -> Self::El { T::zero() }
    fn one(&self) -> Self::El { T::one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { lhs == rhs }

    default fn is_integral(&self) -> bool { false }
    default fn is_euclidean(&self) -> bool { false }
    default fn is_field(&self) -> bool { false }

    default fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: Self::El) -> (Self::El, Self::El) { 
        panic!("Not a euclidean domain!");
    }

    default fn div(&self, _lhs: Self::El, _rhs: Self::El) -> Self::El { 
        panic!("Not a field!");
    }

    default fn euclidean_div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).0
    }
    
    default fn euclidean_rem(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).1
    }
}

impl<T> Ring for StaticRing<RingAxiomsIntegralRing, T> 
    where T: RingEl<Axioms = RingAxiomsIntegralRing>
{
    fn is_integral(&self) -> bool { true }
}

impl<T> Ring for StaticRing<RingAxiomsEuclideanRing, T> 
    where T: EuclideanRingEl
{
    fn is_integral(&self) -> bool { true }
    fn is_euclidean(&self) -> bool { true }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: Self::El) -> (Self::El, Self::El) { 
        let quo = lhs.div_rem(rhs);
        return (quo, lhs);
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        lhs / rhs
    }
    
    fn euclidean_rem(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        lhs % rhs
    }
}

impl<T> Ring for StaticRing<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_integral(&self) -> bool { true }
    fn is_field(&self) -> bool { true }

    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        lhs / rhs
    }
}

#[test]
fn test_pow() {
    assert_eq!(81 * 81 * 3, StaticRing::<RingAxiomsEuclideanRing, i64>::RING.pow(3, 9));
}