use super::*;
use super::super::fq::*;
use super::super::finite_field_sqrt;

impl<K: FiniteRing + HashableElRing> EllipticCurve<WrappingRing<K>> {

    pub fn random_affine_point<F>(&self, rng: F) -> EllipticCurvePoint<WrappingRing<K>>
        where F: FnMut() -> u32
    {
        self.random_affine_point_extension(self.base_field(), rng)
    }

    pub fn random_affine_point_extension<L, F>(&self, ext_field: &WrappingRing<L>, mut rng: F) -> EllipticCurvePoint<WrappingRing<L>>
        where F: FnMut() -> u32, L: FiniteRing + CanonicalEmbeddingInfo<K> + HashableElRing
    {
        assert!(ext_field.is_field().can_use());
        loop {
            let x = ext_field.random_element(&mut rng);
            if let Some(y) = finite_field_sqrt::sqrt(x.pow(3) + &x * ext_field.embed(self.base_field(), self.A.clone()) + ext_field.embed(self.base_field(), self.B.clone()), &ext_field) {
                return EllipticCurvePoint::Affine(x, y);
            }
        }
    }
}