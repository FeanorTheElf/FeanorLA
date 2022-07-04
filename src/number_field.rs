use super::ring::*;
use super::ring_property::*;
use super::embedding::*;
use super::ring_extension::simple_extension::*;
use super::integer::*;
use super::rational::*;
use super::fraction_field::*;
use super::poly::uni_var::*;

pub trait NumberField: Ring + CanonicalIsomorphismInfo<SimpleRingExtension<FieldOfFractions<BigIntRing>>> + RingExtension<BaseRing: RationalField> {

    fn is_number_field(&self) -> RingPropValue;
    fn degree(&self) -> usize;
    fn primitive_element(&self) -> El<Self>;
    fn mipo(&self, el: &El<Self>) -> El<PolyRing<<Self as RingExtension>::BaseRing>>;
}

impl<'a, K: NumberField> CanonicalEmbeddingInfo<SimpleRingExtension<FieldOfFractions<BigIntRing>>> for &'a K {

    fn has_embedding(&self, from: &SimpleRingExtension<FieldOfFractions<BigIntRing>>) -> RingPropValue {
        K::has_embedding(*self, from)
    }

    fn embed(&self, from: &SimpleRingExtension<FieldOfFractions<BigIntRing>>, el: El<SimpleRingExtension<FieldOfFractions<BigIntRing>>>) -> El<Self> {
        K::embed(*self, from, el)
    }
}

impl<'a, K: NumberField> CanonicalIsomorphismInfo<SimpleRingExtension<FieldOfFractions<BigIntRing>>> for &'a K {

    fn has_isomorphism(&self, from: &SimpleRingExtension<FieldOfFractions<BigIntRing>>) -> RingPropValue {
        K::has_isomorphism(*self, from)
    }

    fn preimage(&self, from: &SimpleRingExtension<FieldOfFractions<BigIntRing>>, el: El<Self>) -> El<SimpleRingExtension<FieldOfFractions<BigIntRing>>> {
        K::preimage(*self, from, el)
    }
}

impl<'a, K: NumberField> NumberField for &'a K {

    fn is_number_field(&self) -> RingPropValue { (**self).is_number_field() }
    fn degree(&self) -> usize { (**self).degree() }
    fn primitive_element(&self) -> El<Self> { (**self).primitive_element() }
    fn mipo(&self, el: &El<Self>) -> El<PolyRing<<Self as RingExtension>::BaseRing>> { (**self).mipo(el) }
}

