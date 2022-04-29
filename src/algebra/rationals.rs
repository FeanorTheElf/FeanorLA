use super::super::ring::*;
use super::super::embedding::*;
use super::super::primitive::*;
use super::super::wrapper::*;
use super::rat::*;
use super::integer::*;

pub trait RationalField: Ring + CanonicalIsomorphismInfo<StaticRing<r64>> + CanonicalEmbeddingInfo<Self::UnderlyingIntegers> {

    type UnderlyingIntegers: IntegerRing;

    fn num(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El;
    fn den(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El;
    fn underlying_integers(&self) -> Self::UnderlyingIntegers;
}

impl RationalField for StaticRing<r64> {
    
    type UnderlyingIntegers = StaticRing<i64>;

    fn num(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El {
        r64::num(el)
    }

    fn den(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El {
        r64::num(el)
    }

    fn underlying_integers(&self) -> Self::UnderlyingIntegers {
        i64::RING
    }
}

impl<R: RationalField> RationalField for WrappingRing<R> {
    
    type UnderlyingIntegers = WrappingRing<R::UnderlyingIntegers>;

    fn num(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El {
        self.wrapped_ring().underlying_integers().bind_by_value(self.wrapped_ring().num(el.val()))
    }

    fn den(&self, el: &Self::El) -> <Self::UnderlyingIntegers as Ring>::El {
        self.wrapped_ring().underlying_integers().bind_by_value(self.wrapped_ring().den(el.val()))
    }

    fn underlying_integers(&self) -> Self::UnderlyingIntegers {
        self.wrapped_ring().underlying_integers().bind_ring_by_value()
    }
}