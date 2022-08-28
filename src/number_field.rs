use super::ring::*;
use super::ring_property::*;
use super::embedding::*;
use super::ring_extension::simple_extension::*;
use super::integer::*;
use super::rational::*;
use super::fraction_field::fraction_field_impl::*;
use super::poly::uni_var::*;

pub trait NumberField: Ring + CanonicalIsomorphismInfo<SimpleRingExtension<FractionFieldImpl<BigIntRing>>> + RingExtension<BaseRing: RationalField> {

    fn is_number_field(&self) -> RingPropValue;
    fn degree(&self) -> usize;
    fn primitive_element(&self) -> El<Self>;
    fn mipo(&self, el: &El<Self>) -> El<PolyRingImpl<<Self as RingExtension>::BaseRing>>;
}

impl<'a, K: NumberField> CanonicalEmbeddingInfo<SimpleRingExtension<FractionFieldImpl<BigIntRing>>> for &'a K {

    fn has_embedding(&self, from: &SimpleRingExtension<FractionFieldImpl<BigIntRing>>) -> RingPropValue {
        K::has_embedding(*self, from)
    }

    fn embed(&self, from: &SimpleRingExtension<FractionFieldImpl<BigIntRing>>, el: El<SimpleRingExtension<FractionFieldImpl<BigIntRing>>>) -> El<Self> {
        K::embed(*self, from, el)
    }
}

impl<'a, K: NumberField> CanonicalIsomorphismInfo<SimpleRingExtension<FractionFieldImpl<BigIntRing>>> for &'a K {

    fn has_isomorphism(&self, from: &SimpleRingExtension<FractionFieldImpl<BigIntRing>>) -> RingPropValue {
        K::has_isomorphism(*self, from)
    }

    fn preimage(&self, from: &SimpleRingExtension<FractionFieldImpl<BigIntRing>>, el: El<Self>) -> El<SimpleRingExtension<FractionFieldImpl<BigIntRing>>> {
        K::preimage(*self, from, el)
    }
}

impl<'a, K: NumberField> NumberField for &'a K {

    fn is_number_field(&self) -> RingPropValue { (**self).is_number_field() }
    fn degree(&self) -> usize { (**self).degree() }
    fn primitive_element(&self) -> El<Self> { (**self).primitive_element() }
    fn mipo(&self, el: &El<Self>) -> El<PolyRingImpl<<Self as RingExtension>::BaseRing>> { (**self).mipo(el) }
}

