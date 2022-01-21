#![allow(non_snake_case)]

use super::super::alg::*;
use super::super::la::mat::*;
use super::ring_ext::*;
use super::poly::*;

pub struct EllipticCurve<K: Ring> {
    base_field: K,
    A: K::El,
    B: K::El
}

pub type CoordRing<K> = SimpleRingExtension<PolyRing<K>, VectorOwned<<PolyRing<K> as Ring>::El>>;

impl<K: Ring> EllipticCurve<K> {

    pub fn new(base_field: K, A: K::El, B: K::El) -> Self {
        debug_assert!(base_field.is_field());
        EllipticCurve {
            base_field, A, B
        }
    }

    ///
    /// Returns the coordinate ring together with the `x` and `y` elements of
    /// this elliptic curve.
    /// 
    pub fn coordinate_ring(&self) -> (CoordRing<K>, <CoordRing<K> as Ring>::El, <CoordRing<K> as Ring>::El) {
        let poly_ring = PolyRing::adjoint(self.base_field.clone(), "X");
        let mipo = VectorOwned::from_array([poly_ring.zero(), 
            poly_ring.add(
                poly_ring.pow(&poly_ring.unknown(), 3),
                poly_ring.add(
                    poly_ring.mul(poly_ring.from(self.A.clone()), poly_ring.unknown()),
                    poly_ring.from(self.B.clone())
                )
            )
        ]);
        let poly_ring_x = poly_ring.unknown();
        let result = SimpleRingExtension::new(poly_ring, Vector::new(mipo));
        let X = result.from(poly_ring_x);
        let Y = result.generator();
        return (result, X, Y);
    }

    pub fn j_invariant(&self) -> K::El {
        let A_cubed = self.base_field.pow(&self.A, 3);
        return self.base_field.div(
            self.base_field.mul(self.base_field.from_z(1728 * 4), A_cubed.clone()), 
            &self.base_field.add(
                self.base_field.mul(self.base_field.from_z(27), self.base_field.pow(&self.B, 2)),
                self.base_field.mul(self.base_field.from_z(4), A_cubed)
            )
        );
    }

    pub fn is_isomorphic(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_field.eq(&self.j_invariant(), &rhs.j_invariant())
    }
}