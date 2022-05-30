use super::*;
use super::super::fq::*;
use super::super::finite_field_sqrt;

impl<K: FiniteRing + HashableElRing> EllipticCurve<WrappingRing<K>> {

    pub fn random_affine_point<F>(&self, mut rng: F) -> EllipticCurvePoint<WrappingRing<K>>
        where F: FnMut() -> u32
    {
        loop {
            let x = self.base_field().random_element(&mut rng);
            if let Some(y) = finite_field_sqrt::sqrt(x.pow(3) + &x * &self.A + &self.B, self.base_field()) {
                return EllipticCurvePoint::Affine(x, y);
            }
        }
    }
}