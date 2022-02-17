#![allow(non_snake_case)]
use super::super::ring::*;
use super::super::embedding::*;
use super::super::la::mat::*;
use super::ring_ext::*;
use super::poly::*;
use super::fractions::*;

pub struct EllipticCurve<K: Ring> {
    base_ring: K,
    A: K::El,
    B: K::El
}

pub type CoordRing<K> = SimpleRingExtension<PolyRing<K>, VectorOwned<<PolyRing<K> as Ring>::El>>;
pub type FunctionField<K> = FieldOfFractions<CoordRing<K>>;

pub enum EllipticCurvePoint<K> 
    where K: Ring
{
    Affine(K::El, K::El), Infinity
}

impl<K: Ring> EllipticCurve<K> {

    pub fn new(base_ring: K, A: K::El, B: K::El) -> Self {
        EllipticCurve {
            base_ring, A, B
        }
    }

    ///
    /// Returns the coordinate ring together with the `x` and `y` elements of
    /// this elliptic curve.
    /// 
    pub fn coordinate_ring(&self) -> (CoordRing<K>, <CoordRing<K> as Ring>::El, <CoordRing<K> as Ring>::El) {
        let poly_ring = PolyRing::adjoint(self.base_ring.clone(), "X");
        let mipo = Vector::from_array([poly_ring.zero(), 
            poly_ring.add(
                poly_ring.pow(&poly_ring.unknown(), 3),
                poly_ring.add(
                    poly_ring.mul(poly_ring.from(self.A.clone()), poly_ring.unknown()),
                    poly_ring.from(self.B.clone())
                )
            )
        ]);
        let poly_ring_x = poly_ring.unknown();
        let result = SimpleRingExtension::new(poly_ring, mipo);
        let X = result.from(poly_ring_x);
        let Y = result.generator();
        return (result, X, Y);
    }

    pub fn function_field(&self) -> (FunctionField<K>, <FunctionField<K> as Ring>::El, <FunctionField<K> as Ring>::El) {
        assert!(self.base_ring.is_field().can_use());
        let (ring, x, y) = self.coordinate_ring();
        let field = FieldOfFractions::new(ring.clone());
        let (x, y) = {
            let incl = field.embedding(&ring);
            (incl(x), incl(y))
        };
        return (field, x, y);
    }

    pub fn j_invariant(&self) -> K::El {
        assert!(self.base_ring.is_field().can_use());
        let A_cubed = self.base_ring.pow(&self.A, 3);
        return self.base_ring.div(
            self.base_ring.mul(self.base_ring.from_z(1728 * 4), A_cubed.clone()), 
            &self.base_ring.add(
                self.base_ring.mul(self.base_ring.from_z(27), self.base_ring.pow(&self.B, 2)),
                self.base_ring.mul(self.base_ring.from_z(4), A_cubed)
            )
        );
    }

    pub fn is_isomorphic(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_ring.eq(&self.j_invariant(), &rhs.j_invariant())
    }

    pub fn point_add(&self, a: EllipticCurvePoint<K>, b: EllipticCurvePoint<K>) -> EllipticCurvePoint<K> {
        assert!(self.base_ring.is_field().can_use());
        match (a, b) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => EllipticCurvePoint::Infinity,
            (EllipticCurvePoint::Affine(x, y), EllipticCurvePoint::Infinity) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Affine(x, y)) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if self.base_ring.eq(&x1, &x2) && self.base_ring.eq(&y1, &y2) => {
                let lambda = self.base_ring.div(
                    self.base_ring.add_ref(
                        self.base_ring.mul(self.base_ring.from_z(3), self.base_ring.mul_ref(&x1, &x1)),
                        &self.A
                    ),
                    &self.base_ring.mul(self.base_ring.from_z(2), y1)
                );
                let x = self.base_ring.add(
                    self.base_ring.neg(self.base_ring.add_ref(x1, &x2)), 
                    self.base_ring.mul_ref(&lambda, &lambda)
                );
                EllipticCurvePoint::Affine(
                    x.clone(),
                    self.base_ring.add(
                        self.base_ring.neg(y2), 
                        self.base_ring.mul(lambda, self.base_ring.sub(x2, x))
                    )
                )
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if !self.base_ring.eq(&x1, &x2) => {
                let lambda = self.base_ring.div(
                    self.base_ring.sub_ref_fst(&y1, y2),
                    &self.base_ring.sub_ref_fst(&x1, x2.clone())
                );
                let x = self.base_ring.add(
                    self.base_ring.neg(self.base_ring.add_ref(x2, &x1)), 
                    self.base_ring.mul_ref(&lambda, &lambda)
                );
                EllipticCurvePoint::Affine(
                    x.clone(),
                    self.base_ring.add(
                        self.base_ring.neg(y1), 
                        self.base_ring.mul(lambda, self.base_ring.sub(x1, x))
                    )
                )
            },
            (EllipticCurvePoint::Affine(_, _), EllipticCurvePoint::Affine(_, _)) => EllipticCurvePoint::Infinity
        }
    }

    pub fn extend_ring<'a, S, F>(self, new_ring: S, incl: F) -> (EllipticCurve<S>, impl Fn(EllipticCurvePoint<K>) -> EllipticCurvePoint<S>)
        where S: Ring, F: Fn(K::El) -> S::El
    {
        (EllipticCurve {
            base_ring: new_ring,
            A: incl(self.A),
            B: incl(self.B)
        },
        move |point| match point {
            EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(incl(x), incl(y)),
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity
        })
    } 
}
