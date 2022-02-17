use super::super::ring::*;
use super::super::bigint::*;

pub trait FiniteRing : DivisibilityInfoRing {

    fn characteristic(&self) -> BigInt;
    fn size(&self) -> BigInt;
}

pub mod zn_big;
pub mod zn_small;
pub mod fq_small;