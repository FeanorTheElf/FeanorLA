use super::super::ring::*;
use super::super::bigint::*;

pub trait FiniteRing : DivisibilityInfoRing {

    fn characteristic(&self) -> BigInt;
    fn size(&self) -> BigInt;
}

impl<'a, R> FiniteRing for &'a R
    where R: FiniteRing
{
    fn characteristic(&self) -> BigInt { (**self).characteristic() }
    fn size(&self) -> BigInt { (**self).size() }
}

pub mod zn_big;
pub mod zn_small;
pub mod fq_small;