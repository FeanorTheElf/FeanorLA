use super::prelude::*;
use super::fq::*;
use super::integer::*;
use super::la::vec::*;
use super::finite_extension::finite_extension_impl::*;

use std::hash::Hasher;
use std::collections::hash_map::DefaultHasher;

///
/// Computes a square root in a finite field of a quadratic residue.
/// 
/// # Current approach
/// 
/// To find whether the value is a quadratic residue, compute `x^((q - 1)/2)`.
/// Then x is a quadratic residue if and only if this is not -1, as the multiplicative
/// group of a finite field is cyclic.
/// 
/// To find the root, consider the ring
/// ```text
/// Fq[T]/(T^2 - x) ~ Fq x Fq
/// ```
/// Now we take random elements y in that ring and hope that the the components in the
/// CRT representation on the right-hand side have order dividing resp. not dividing
/// `(q - 1)/2`. In this case, we find that `gcd((T^2 - x), y^((q - 1)/2))` gives the
/// root of x.
/// 
pub fn sqrt<F: FiniteRing + HashableElRing>(x: El<F>, field: F) -> Option<El<F>> {
    assert!(field.is_field().can_use());
    assert!(field.characteristic() != 2);
    
    if field.is_zero(&x) {
        return Some(field.zero());
    }

    let mut hasher = DefaultHasher::new();
    field.hash(&mut hasher, &x);
    let mut rng = oorandom::Rand32::new(hasher.finish());

    let n = BigInt::RING.quotient(&BigInt::RING.sub(field.size(), BigInt::RING.one()), &BigInt::from(2)).unwrap();
    if field.is_neg_one(&field.pow_big(&x, &n)) {
        return None;
    }
    let ring: FiniteExtensionImpl<_, _, VectorArray<_, 2>> = FiniteExtensionImpl::new(&field, Vector::from_array([x, field.zero()]), "α");
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
    assert!(field.is_eq(&a, &root_a2) || field.is_eq(&field.neg(a), &root_a2));
    assert!(field.is_zero(&sqrt(field.zero(), &field).unwrap()));
}

#[test]
fn test_fq_sqrt() {
    let field = F49;
    let a = field.from_z(11);
    let root_a2 = sqrt(field.pow(&a, 2), &field).unwrap();
    assert!(field.is_eq(&a, &root_a2) || field.is_eq(&field.neg(a), &root_a2));
}