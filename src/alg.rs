use std::ops::{ Add, Mul, Sub, Neg, Div, Rem, AddAssign, MulAssign, DivAssign, SubAssign, RemAssign };

pub trait Zero {
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

pub trait One {
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

pub trait Semiring: Sized + Add<Output = Self> + Mul<Output = Self> + AddAssign + MulAssign + PartialEq + Zero + One {}

pub trait Ring: Semiring + Neg<Output = Self> + Sub<Output = Self> + SubAssign {}

pub trait CommutativeRing: Ring {}

pub trait Field: CommutativeRing + Div<Output = Self> + DivAssign {}

impl Semiring for u8 {}
impl Semiring for u16 {}
impl Semiring for u32 {}
impl Semiring for u64 {}
impl Semiring for u128 {}

impl Semiring for i8 {}
impl Semiring for i16 {}
impl Semiring for i32 {}
impl Semiring for i64 {}
impl Semiring for i128 {}
impl Ring for i8 {}
impl Ring for i16 {}
impl Ring for i32 {}
impl Ring for i64 {}
impl Ring for i128 {}
impl CommutativeRing for i8 {}
impl CommutativeRing for i16 {}
impl CommutativeRing for i32 {}
impl CommutativeRing for i64 {}
impl CommutativeRing for i128 {}

impl Semiring for f32 {}
impl Semiring for f64 {}
impl Ring for f32 {}
impl Ring for f64 {}
impl CommutativeRing for f32 {}
impl CommutativeRing for f64 {}
impl Field for f32 {}
impl Field for f64 {}

pub trait Integer: Add + Sub + Mul + Neg + Rem + AddAssign + MulAssign + SubAssign + RemAssign {}

impl Integer for i8 {}
impl Integer for i16 {}
impl Integer for i32 {}
impl Integer for i64 {}
impl Integer for i128 {}