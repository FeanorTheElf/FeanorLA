use super::super::prelude::*;
use super::*;

pub trait GaloisExtensionInfo: RingExtension {

    type GaloisGroupEl;

    fn is_field_extension(&self) -> RingPropValue {
        self.is_field() & self.base_ring().is_field()
    }

    fn base_field(&self) -> &Self::BaseRing {
        assert!(self.is_field_extension().can_use());
        self.base_ring()
    }

    fn is_galois_extension(&self) -> RingPropValue;

    fn apply(&self, x: El<Self>, automorphism: Self::GaloisGroupEl) -> El<Self>;
}