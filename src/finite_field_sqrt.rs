use super::prelude::*;
use super::fq::*;
use super::ring_extension::simple_extension::*;

use std::hash::Hasher;
use std::collections::hash_map::DefaultHasher;

pub fn sqrt<F: FiniteRing + HashableElRing>(x: El<F>, field: F) -> Option<El<F>> {
    assert!(field.is_field().can_use());
    assert!(field.characteristic() != 2);

    let mut hasher = DefaultHasher::new();
    field.hash(&mut hasher, &x);
    let mut rng = oorandom::Rand32::new(hasher.finish());

    let n = BigInt::RING.quotient(&(field.size() - 1), &BigInt::from(2)).unwrap();
    if !field.is_one(&field.pow_big(&x, &n)) {
        return None;
    }
    let ring: SimpleRingExtension<_, _, VectorArray<_, 2>> = SimpleRingExtension::new(&field, Vector::from_array([x, field.zero()]));
    loop {
        let g = ring.random_element(|| rng.rand_u32());
        let h = ring.pow_big(&g, &n);
        if !ring.is_one(&h) && !ring.is_neg_one(&h) {
            let factor = ring.sub(h, ring.one());
            return Some(field.neg(field.div(factor.at(0).clone(), factor.at(1))));
        }
    }
}

#[cfg(test)]
use super::fq::fq_small::*;
#[cfg(test)]
use super::fq::zn_small::*;

#[test]
fn test_zn_sqrt() {
    let field = ZnEl::<73>::RING;
    let a = field.from_z(42);
    let root_a2 = sqrt(field.pow(&a, 2), &field).unwrap();
    assert!(field.eq(&a, &root_a2) || field.eq(&field.neg(a), &root_a2));
}

#[test]
fn test_fq_sqrt() {
    let field = F49;
    let a = field.from_z(11);
    let root_a2 = sqrt(field.pow(&a, 2), &field).unwrap();
    assert!(field.eq(&a, &root_a2) || field.eq(&field.neg(a), &root_a2));
}