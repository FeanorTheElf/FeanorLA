use super::ring::*;
use super::ring_property::*;
use super::embedding::*;
use super::finite_extension::finite_extension_impl::*;
use super::integer::bigint_soo::*;
use super::rational::*;
use super::fraction_field::fraction_field_impl::*;
use super::poly::uni_var::*;

pub trait NumberField: Ring + CanonicalIsomorphismInfo<FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>> + RingExtension<BaseRing: RationalField> {

    fn is_number_field(&self) -> RingPropValue;
    fn degree(&self) -> usize;
    fn primitive_element(&self) -> El<Self>;
    fn mipo(&self, el: &El<Self>) -> El<PolyRingImpl<<Self as RingExtension>::BaseRing>>;
}

impl<'a, K: NumberField> CanonicalEmbeddingInfo<FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>> for &'a K {

    fn has_embedding(&self, from: &FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>) -> RingPropValue {
        K::has_embedding(*self, from)
    }

    fn embed(&self, from: &FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>, el: El<FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>>) -> El<Self> {
        K::embed(*self, from, el)
    }
}

impl<'a, K: NumberField> CanonicalIsomorphismInfo<FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>> for &'a K {

    fn has_isomorphism(&self, from: &FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>) -> RingPropValue {
        K::has_isomorphism(*self, from)
    }

    fn preimage(&self, from: &FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>, el: El<Self>) -> El<FiniteExtensionImpl<FractionFieldImpl<BigIntSOORing>>> {
        K::preimage(*self, from, el)
    }
}

impl<'a, K: NumberField> NumberField for &'a K {

    fn is_number_field(&self) -> RingPropValue { (**self).is_number_field() }
    fn degree(&self) -> usize { (**self).degree() }
    fn primitive_element(&self) -> El<Self> { (**self).primitive_element() }
    fn mipo(&self, el: &El<Self>) -> El<PolyRingImpl<<Self as RingExtension>::BaseRing>> { (**self).mipo(el) }
}

