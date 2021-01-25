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

pub trait Ring: Sized + Add<Output = Self> + Mul<Output = Self> + AddAssign + PartialEq + Zero + One + Neg<Output = Self> + Sub<Output = Self> + SubAssign {}

pub trait EuclideanRing : Ring + Rem<Output = Self> + RemAssign {}

pub trait Field: Ring + MulAssign + Div<Output = Self> + DivAssign {}

impl Ring for i8 {}
impl Ring for i16 {}
impl Ring for i32 {}
impl Ring for i64 {}
impl Ring for i128 {}
impl EuclideanRing for i8 {}
impl EuclideanRing for i16 {}
impl EuclideanRing for i32 {}
impl EuclideanRing for i64 {}
impl EuclideanRing for i128 {}

impl Ring for f32 {}
impl Ring for f64 {}
impl Field for f32 {}
impl Field for f64 {}

pub trait Integer: Clone + Eq + Ord + EuclideanRing + /* truncating division */ Div<Output = Self> + DivAssign {}

impl Integer for i8 {}
impl Integer for i16 {}
impl Integer for i32 {}
impl Integer for i64 {}
impl Integer for i128 {}

pub trait Float: Clone + PartialEq + PartialOrd + Field {

    ///
    /// This function should yield a non-negative floating point value that gives an estimate
    /// for the size of the float. Namely, it should hold that for any float type F have
    /// 0_F.stability_abs() == 0 and if |a| >> |b| (much bigger), then a.stability_abs() > b.stability_abs()
    /// 
    /// This can be used by algorithms to optimize numerical stability
    /// 
    fn stability_abs(&self) -> f32;

}

impl Float for f32 {
    
    fn stability_abs(&self) -> f32 {
        self.abs()
    }
}

impl Float for f64 {
    
    fn stability_abs(&self) -> f32 {
        (*self as f32).abs()
    }
}