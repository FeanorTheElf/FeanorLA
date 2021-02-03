use std::ops::{ Add, Mul, Sub, Neg, Div, Rem, AddAssign, MulAssign, DivAssign, SubAssign, RemAssign };

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
/// Multiplication must be commutative
/// 
pub trait RingEl: Sized + Add<Output = Self> + Mul<Output = Self> + AddAssign + PartialEq + Zero + One + Neg<Output = Self> + Sub<Output = Self> + SubAssign {}

pub trait IntegralRingEl: RingEl {}

///
/// Division must satisfy the invariant x = (x/y) * y + (x%y)
/// 
pub trait EuclideanRingEl: IntegralRingEl + Rem<Output = Self> + RemAssign + Div<Output = Self> + DivAssign {

    ///
    /// Computes (returned, self) := (self / rhs, self % rhs) and returns returned.
    /// Can be faster than computing both separately
    /// 
    fn div_rem(&mut self, rhs: Self) -> Self;

}

pub trait FieldEl: IntegralRingEl + MulAssign + Div<Output = Self> + DivAssign {}

impl RingEl for i8 {}
impl RingEl for i16 {}
impl RingEl for i32 {}
impl RingEl for i64 {}
impl RingEl for i128 {}
impl IntegralRingEl for i8 {}
impl IntegralRingEl for i16 {}
impl IntegralRingEl for i32 {}
impl IntegralRingEl for i64 {}
impl IntegralRingEl for i128 {}

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

impl RingEl for f32 {}
impl RingEl for f64 {}
impl IntegralRingEl for f32 {}
impl IntegralRingEl for f64 {}
impl FieldEl for f32 {}
impl FieldEl for f64 {}

///
/// Division means truncating division
/// 
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
    type El: Sized;

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El;
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El;
    fn neg(&self, val: Self::El) -> Self::El;
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool;

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.add(lhs, self.neg(rhs))
    }

    fn pow(&self, basis: Self::El, exp: u64) -> Self::El 
        where Self::El: Clone
    {
        let mut power = basis;
        let mut result = self.one();
        for i in 0..(64 - exp.leading_zeros()) {
            if ((exp >> i) & 1) != 0 {
                result = self.mul(result, power.clone());
            }
            power = self.mul(power.clone(), power);
        }
        return result;
    }
}

pub trait IntegralRing: Ring {}

pub trait EuclideanRing: IntegralRing {

    fn div_rem(&self, lhs: Self::El, rhs: Self::El) -> (Self::El, Self::El);

    fn rem(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.div_rem(lhs, rhs).1 }
    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.div_rem(lhs, rhs).0 }
}

pub trait Field: IntegralRing {
    
    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El;
}

pub struct StaticRing<T: RingEl> {
    element: std::marker::PhantomData<T>
}

impl<T: RingEl> StaticRing<T> {
    
    pub const RING: StaticRing<T> = StaticRing { element: std::marker::PhantomData };

}

impl<T: RingEl> Ring for StaticRing<T> {
    type El = T;

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs + rhs }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs * rhs }
    fn neg(&self, val: Self::El) -> Self::El { -val }
    fn zero(&self) -> Self::El { T::zero() }
    fn one(&self) -> Self::El { T::one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { lhs == rhs }
}

impl<T: IntegralRingEl> IntegralRing for StaticRing<T> {}

impl<T: EuclideanRingEl> EuclideanRing for StaticRing<T> {

    fn div_rem(&self, mut lhs: Self::El, rhs: Self::El) -> (Self::El, Self::El) { 
        let quo = lhs.div_rem(rhs); 
        (quo, lhs)
    }

    fn rem(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs % rhs }
    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs / rhs }
}

impl<T: FieldEl> Field for StaticRing<T> {

    fn div(&self, lhs: Self::El, rhs: Self::El) -> Self::El { lhs / rhs }
}

#[test]
fn test_pow() {
    assert_eq!(81 * 81 * 3, StaticRing::<i64>::RING.pow(3, 9));
}