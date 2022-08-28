pub mod primitive_rational;

pub use primitive_rational::r64;
use super::embedding::*;
use super::integer::*;
use super::primitive::*;
use super::fraction_field::*;

pub trait RationalField: FractionField<BaseRing: IntegerRing> + CanonicalIsomorphismInfo<StaticRing<r64>> {}

impl<'a, R: RationalField> RationalField for &'a R {}