pub mod primitive_rational;

pub use primitive_rational::r64;
use super::ring::*;
use super::embedding::*;
use super::integer::*;
use super::primitive::*;

pub trait RationalField: Ring + CanonicalIsomorphismInfo<StaticRing<r64>> + RingExtension<BaseRing: IntegerRing> {

    fn num(&self, el: &El<Self>) -> El<Self::BaseRing>;
    fn den(&self, el: &El<Self>) -> El<Self::BaseRing>;
}

impl RationalField for StaticRing<r64> {
    
    fn num(&self, el: &El<Self>) -> El<Self::BaseRing> {
        r64::num(el)
    }

    fn den(&self, el: &El<Self>) -> El<Self::BaseRing> {
        r64::num(el)
    }
}

impl<'a, R: RationalField> RationalField for &'a R {

    fn num(&self, el: &El<Self>) -> El<Self::BaseRing> { (**self).num(el) }
    fn den(&self, el: &El<Self>) -> El<Self::BaseRing> { (**self).den(el) }
}