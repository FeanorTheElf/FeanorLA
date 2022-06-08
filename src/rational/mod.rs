pub mod primitive_rational;

pub use primitive_rational::r64;
use super::ring::*;
use super::embedding::*;
use super::primitive::*;
use super::wrapper::*;
use super::integer::*;

pub trait RationalField: Ring + CanonicalIsomorphismInfo<StaticRing<r64>> + CanonicalEmbeddingInfo<Self::UnderlyingIntegers> {

    type UnderlyingIntegers: IntegerRing;

    fn num(&self, el: &Self::El) -> El<Self::UnderlyingIntegers>;
    fn den(&self, el: &Self::El) -> El<Self::UnderlyingIntegers>;
    fn underlying_integers(&self) -> Self::UnderlyingIntegers;
}

impl RationalField for StaticRing<r64> {
    
    type UnderlyingIntegers = StaticRing<i64>;

    fn num(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        r64::num(el)
    }

    fn den(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        r64::num(el)
    }

    fn underlying_integers(&self) -> Self::UnderlyingIntegers {
        i64::RING
    }
}

impl<R: RationalField> RationalField for WrappingRing<R> {
    
    type UnderlyingIntegers = WrappingRing<R::UnderlyingIntegers>;

    fn num(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        self.wrapped_ring().underlying_integers().bind_by_value(self.wrapped_ring().num(el.val()))
    }

    fn den(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        self.wrapped_ring().underlying_integers().bind_by_value(self.wrapped_ring().den(el.val()))
    }

    fn underlying_integers(&self) -> Self::UnderlyingIntegers {
        self.wrapped_ring().underlying_integers().bind_ring_by_value()
    }
}