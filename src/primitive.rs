use super::ring::*;
use super::bigint::*;
use super::embedding::*;
use super::ring_property::*;
use std::ops::{ 
    Add, Mul, Sub, Neg, Div,
    AddAssign, MulAssign, SubAssign, DivAssign
};

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

///
/// Note that the Ring axioms are mutually exclusive, so a ring providing
/// euclidean division should be RingAxiomsEuclideanRing but not 
/// RingAxiomsIntegralRing even though it is integral
/// 
pub trait RingAxioms: 'static {
    fn is_integral() -> bool;
    fn is_euclidean() -> bool;
    fn is_field() -> bool;
    fn is_divisibility_computable() -> bool;
    fn is_ufd() -> bool;
}

#[derive(Debug)]
pub struct RingAxiomsEuclidean;
#[derive(Debug)]
pub struct RingAxiomsField;
#[derive(Debug)]
pub struct RingAxiomsIntegral;
#[derive(Debug)]
pub struct RingAxiomsGeneral;

impl RingAxioms for RingAxiomsIntegral {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { false }
    fn is_field() -> bool { false }
    fn is_divisibility_computable() -> bool { false }
    fn is_ufd() -> bool { false }
}

impl RingAxioms for RingAxiomsEuclidean {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { true }
    fn is_field() -> bool { false }
    fn is_divisibility_computable() -> bool { true }
    fn is_ufd() -> bool { false }
}

impl RingAxioms for RingAxiomsField {
    fn is_integral() -> bool { true }
    fn is_euclidean() -> bool { true }
    fn is_field() -> bool { true }
    fn is_divisibility_computable() -> bool { true }
    fn is_ufd() -> bool { true }
}

impl RingAxioms for RingAxiomsGeneral {
    fn is_integral() -> bool { false }
    fn is_euclidean() -> bool { false }
    fn is_field() -> bool { false }
    fn is_divisibility_computable() -> bool { false }
    fn is_ufd() -> bool { false }
}

pub trait RingEl: 
    Clone + Sized + Add<Output = Self> + Mul<Output = Self> + 
    PartialEq + Zero + One + Neg<Output = Self> + 
    Sub<Output = Self> + std::fmt::Debug +
    AddAssign + MulAssign + SubAssign
{
    type Axioms: RingAxioms;
    type RingType: Ring<El = Self> + SingletonRing;
    const RING: Self::RingType;

    fn pow(self, exp: u32) -> Self {
        Self::RING.pow(&self, exp)
    }

    ///
    /// Returns the ring associated to this type, i.e. `Self::RING`.
    /// This function is useful in case you do not explicitly want
    /// to write down the type of an element, e.g. because it is the
    /// result of nesting ring constructors.
    /// 
    fn ring(&self) -> Self::RingType {
        Self::RING
    }
}

pub trait EuclideanEl: 
    RingEl<Axioms = RingAxiomsEuclidean>
{
    ///
    /// Computes (returned, self) := (self / rhs, self % rhs) and returns returned.
    /// Can be faster than computing both separately
    /// 
    fn div_rem(&mut self, rhs: Self) -> Self;

    fn euclidean_deg(&self) -> BigInt;
}

pub trait FieldEl: 
    RingEl<Axioms = RingAxiomsField> + Div<Output = Self>  + DivAssign
{}

impl RingEl for i8 {
    type Axioms = RingAxiomsEuclidean;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl RingEl for i16 {
    type Axioms = RingAxiomsEuclidean;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl RingEl for i32 {
    type Axioms = RingAxiomsEuclidean;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl RingEl for i64 {
    type Axioms = RingAxiomsEuclidean;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl RingEl for i128 {
    type Axioms = RingAxiomsEuclidean;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}


impl EuclideanEl for i8 {

    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }

    fn euclidean_deg(&self) -> BigInt {
        BigInt::from(self.abs() as i64)
    }
}

impl EuclideanEl for i16 {

    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }

    fn euclidean_deg(&self) -> BigInt {
        BigInt::from(self.abs() as i64)
    }
}

impl EuclideanEl for i32 {

    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }

    fn euclidean_deg(&self) -> BigInt {
        BigInt::from(self.abs() as i64)
    }
}

impl EuclideanEl for i64 {

    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }

    fn euclidean_deg(&self) -> BigInt {
        BigInt::from(self.abs() as i64)
    }
}

impl EuclideanEl for i128 {

    fn div_rem(&mut self, rhs: Self) -> Self { 
        let result = *self / rhs;
        *self %= rhs;
        return result;
    }

    fn euclidean_deg(&self) -> BigInt {
        let upper_part = BigInt::from((self.abs() >> 64) as i64);
        let lower_part = BigInt::from((self.abs() & ((1 << 64) - 1)) as i64);
        return BigInt::RING.add(BigInt::RING.mul(BigInt::power_of_two(64), upper_part), lower_part);
    }
}

impl RingEl for f32 {
    type Axioms = RingAxiomsField;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;

    fn pow(self, exp: u32) -> Self {
        if exp > i32::MAX as u32 {
            let value = f32::powi(self, (exp >> 1) as i32);
            if exp & 1 == 1 {
                value * value * self
            } else {
                value * value
            }
        } else {
            f32::powi(self, exp as i32)
        }
    }
}

impl RingEl for f64 {
    type Axioms = RingAxiomsField;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;

    fn pow(self, exp: u32) -> Self {
        if exp > i32::MAX as u32 {
            let value = f64::powi(self, (exp >> 1) as i32);
            if exp & 1 == 1 {
                value * value * self
            } else {
                value * value
            }
        } else {
            f64::powi(self, exp as i32)
        }
    }
}

impl FieldEl for f32 {}
impl FieldEl for f64 {}

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

impl<Axioms, T> std::fmt::Debug for StaticRingImpl<Axioms, T> 
    where Axioms: RingAxioms, T: RingEl<Axioms = Axioms>
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "StaticRing")
    }
}

impl<A, T> StaticRingImpl<A, T> 
    where A: RingAxioms, T: RingEl<Axioms = A>
{
    pub const RING: StaticRingImpl<A, T> = StaticRingImpl { 
        element: std::marker::PhantomData 
    };
}

trait DisplayOrDebug {
    fn fmt(&self, out: &mut std::fmt::Formatter) -> std::fmt::Result;
}

impl<T: std::fmt::Debug> DisplayOrDebug for T {
    default fn fmt(&self, out: &mut std::fmt::Formatter) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, out)
    }
}

impl<T: std::fmt::Display + std::fmt::Debug> DisplayOrDebug for T {
    fn fmt(&self, out: &mut std::fmt::Formatter) -> std::fmt::Result {
        <Self as std::fmt::Display>::fmt(self, out)
    }
}

impl<T> PartialEq for StaticRingImpl<T::Axioms, T> 
    where T: RingEl
{
    fn eq(&self, _rhs: &Self) -> bool {
        true
    }
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

    fn is_noetherian(&self) -> bool { true }
    default fn is_integral(&self) -> RingPropValue { RingPropValue::False }
    default fn is_field(&self) -> RingPropValue { RingPropValue::False }
    default fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El { 
        panic!("Not a field!");
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        <Self::El as DisplayOrDebug>::fmt(el, f)
    }
}

impl<T> Ring for StaticRingImpl<RingAxiomsEuclidean, T> 
    where T: RingEl<Axioms = RingAxiomsEuclidean>
{
    fn is_integral(&self) -> RingPropValue { RingPropValue::True }
}

impl<T> EuclideanInfoRing for StaticRingImpl<RingAxiomsEuclidean, T> 
    where T: EuclideanEl
{
    fn is_euclidean(&self) -> RingPropValue { 
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { 
        let quo = lhs.div_rem(rhs.clone());
        return (quo, lhs);
    }

    fn euclidean_deg(&self, el: Self::El) -> BigInt {
        el.euclidean_deg()
    }
}

impl<T> DivisibilityInfoRing for StaticRingImpl<RingAxiomsEuclidean, T> 
    where T: EuclideanEl
{
    fn is_divisibility_computable(&self) -> bool {
        true
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if self.is_zero(rhs) && !self.is_zero(lhs) {
            return None;
        } else if self.is_zero(rhs) && self.is_zero(lhs) {
            return Some(self.one());
        }
        let (quo, rem) = self.euclidean_div_rem(lhs.clone(), rhs);
        if self.is_zero(&rem) {
            Some(quo)
        } else {
            None
        }
    }
}

impl<T> EuclideanInfoRing for StaticRingImpl<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_euclidean(&self) -> RingPropValue { 
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { 
        assert!(!self.is_zero(rhs));
        (lhs / rhs.clone(), self.zero())
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        assert!(!self.is_zero(rhs));
        lhs / rhs.clone()
    }
    
    fn euclidean_rem(&self, _lhs: Self::El, rhs: &Self::El) -> Self::El { 
        assert!(!self.is_zero(rhs));
        self.zero()
    }

    fn euclidean_deg(&self, el: Self::El) -> BigInt {
        if self.is_zero(&el) {
            return BigInt::RING.zero();
        } else {
            return BigInt::RING.one();
        }
    }
}

impl<T> Ring for StaticRingImpl<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_integral(&self) -> RingPropValue { RingPropValue::True }
    fn is_field(&self) -> RingPropValue { RingPropValue::True }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        lhs / rhs.clone()
    }
}

impl<T> DivisibilityInfoRing for StaticRingImpl<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_divisibility_computable(&self) -> bool {
        true
    }

    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.is_zero(lhs) || !self.is_zero(rhs)
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if !self.is_zero(lhs) && self.is_zero(rhs) {
            None
        } else if self.is_zero(lhs) && self.is_zero(rhs) {
            Some(self.one())
        } else {
            Some(lhs.clone() / rhs.clone())
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        !self.is_zero(el)
    }
}

impl<T> SingletonRing for StaticRingImpl<T::Axioms, T>
    where T: RingEl
{
    fn singleton() -> Self {
        Self::RING
    }
}

pub type StaticRing<R> = StaticRingImpl<<R as RingEl>::Axioms, R>;

impl<T: RingEl> CanonicalEmbeddingInfo<StaticRing<T>> for StaticRing<T> {

    fn has_embedding(&self, _from: &StaticRing<T>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<T>, el: T) -> T {
        el
    }
}

impl<T: RingEl> CanonicalIsomorphismInfo<StaticRing<T>> for StaticRing<T> {

    fn has_isomorphism(&self, _from: &StaticRing<T>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<T>, el: T) -> T {
        el
    }
}

impl CanonicalEmbeddingInfo<BigIntRing> for StaticRing<i64> {

    fn has_embedding(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntRing, el: BigInt) -> Self::El {
        el.to_int().expect("Overflow when embedding BigInt into i64")
    }
}

impl CanonicalIsomorphismInfo<BigIntRing> for StaticRing<i64> {

    fn has_isomorphism(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntRing, el: i64) -> BigInt {
        BigInt::from(el)
    }
}