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
pub trait Ring: Sized + Add<Output = Self> + Mul<Output = Self> + AddAssign + PartialEq + Zero + One + Neg<Output = Self> + Sub<Output = Self> + SubAssign {}

pub trait IntegralRing: Ring {}

pub trait EuclideanRing : IntegralRing + Rem<Output = Self> + RemAssign {}

pub trait Field: IntegralRing + MulAssign + Div<Output = Self> + DivAssign {}

impl Ring for i8 {}
impl Ring for i16 {}
impl Ring for i32 {}
impl Ring for i64 {}
impl Ring for i128 {}
impl IntegralRing for i8 {}
impl IntegralRing for i16 {}
impl IntegralRing for i32 {}
impl IntegralRing for i64 {}
impl IntegralRing for i128 {}
impl EuclideanRing for i8 {}
impl EuclideanRing for i16 {}
impl EuclideanRing for i32 {}
impl EuclideanRing for i64 {}
impl EuclideanRing for i128 {}

impl Ring for f32 {}
impl Ring for f64 {}
impl IntegralRing for f32 {}
impl IntegralRing for f64 {}
impl Field for f32 {}
impl Field for f64 {}

///
/// Division means truncating division
/// 
pub trait Integer: Clone + Eq + Ord + EuclideanRing + Div<Output = Self> + DivAssign {}

impl Integer for i8 {}
impl Integer for i16 {}
impl Integer for i32 {}
impl Integer for i64 {}
impl Integer for i128 {}

///
/// For types that may have rounding errors
/// 
pub trait Float: Clone + PartialEq + PartialOrd + Field {}

impl Float for f32 {}
impl Float for f64 {}