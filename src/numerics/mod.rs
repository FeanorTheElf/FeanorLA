use super::alg::FieldEl;

///
/// For types that may have rounding errors
/// 
pub trait Float: Clone + PartialEq + PartialOrd + FieldEl {}

impl Float for f32 {}
impl Float for f64 {}


#[cfg(test)]
#[macro_use]
pub mod macros;
pub mod simplex;
pub mod qr;