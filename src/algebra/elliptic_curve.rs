#![allow(non_snake_case)]
use super::super::ring::*;
use super::super::bigint::*;
use super::super::embedding::*;
use super::super::wrapper::*;
use super::super::la::mat::*;
use super::ring_ext::*;
use super::poly::*;
use super::fractions::*;
use super::eea::*;

#[derive(Debug)]
pub struct EllipticCurve<K: Ring> {
    base_ring: K,
    A: K::El,
    B: K::El
}

impl<K: Ring> PartialEq for EllipticCurve<K> {

    fn eq(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_ring.eq(&self.A, &rhs.A) && self.base_ring.eq(&self.B, &rhs.B)
    }
}

impl<K: Ring> std::fmt::Display for EllipticCurve<K> {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Elliptic curve defined by Y^2 = X^3 + {} X + {}", self.base_ring.display(&self.A), self.base_ring.display(&self.B))
    }
}

pub type CoordRing<K> = SimpleRingExtension<PolyRing<K>, VectorOwned<<PolyRing<K> as Ring>::El>>;
pub type FunctionField<K> = FieldOfFractions<CoordRing<K>>;

#[derive(Clone)]
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

    pub fn is_on_curve(&self, point: &EllipticCurvePoint<K>) -> bool {
        match point {
            EllipticCurvePoint::Infinity => true,
            EllipticCurvePoint::Affine(x, y) => self.base_ring.eq(
                &self.base_ring.pow(y, 2),
                &self.base_ring.add(
                    self.base_ring.pow(x, 3), 
                    self.base_ring.add(
                        self.base_ring.mul_ref(&self.A, x),
                        self.B.clone()
                    ))
            )
        }
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

    pub fn points_eq(&self, P: &EllipticCurvePoint<K>, Q: &EllipticCurvePoint<K>) -> bool {
        match (P, Q) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => true,
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) => self.base_ring.eq(x1, x2) && self.base_ring.eq(y1, y2),
            _ => false
        }
    }
}

type QType = FieldOfFractions<BigIntRing>;

impl EllipticCurve<QType> 
{
    ///
    /// Returns an isomorphic elliptic curve E': y^2 = x^3 + Ax + B with A, B in Z
    /// with minimal discriminant and isomorphisms f: E -> E' and f^-1: E' -> E.
    /// 
    pub fn isomorphic_curve_over_z(&self) -> (Self, impl Fn(EllipticCurvePoint<QType>) -> EllipticCurvePoint<QType>, impl Fn(EllipticCurvePoint<QType>) -> EllipticCurvePoint<QType>) {
        let (A_num, A_den) = self.base_ring.reduce(self.A.clone());
        let (B_num, B_den) = self.base_ring.reduce(self.B.clone());

        let Z = BigInt::RING;

        let u_num_A = Z.factor(A_den.clone()).into_iter()
            .map(|(factor, power)| factor.pow((power as u32 - 1) / 4 + 1))
            .product::<RingElWrapper<&BigIntRing>>().into_val();
        let u_num_B = Z.factor(B_den.clone()).into_iter()
            .map(|(factor, power)| factor.pow((power as u32 - 1) / 6 + 1))
            .product::<RingElWrapper<&BigIntRing>>().into_val();
        
        let u_den_A = Z.factor(A_num.clone()).into_iter()
            .map(|(factor, power)| factor.pow(power as u32 / 4))
            .product::<RingElWrapper<&BigIntRing>>().into_val();
        let u_den_B = Z.factor(B_num.clone()).into_iter()
            .map(|(factor, power)| factor.pow(power as u32 / 6))
            .product::<RingElWrapper<&BigIntRing>>().into_val();

        let u = self.base_ring.div(
            self.base_ring.from(lcm(&BigInt::RING, u_num_A, u_num_B)),
            &self.base_ring.from(gcd(&BigInt::RING, u_den_A, u_den_B))
        );
        let u_inv = self.base_ring.div(self.base_ring.one(), &u);
        let base_ring = self.base_ring.clone();
        let base_ring_copy = base_ring.clone();
        return (EllipticCurve {
            base_ring: self.base_ring.clone(),
            A: self.base_ring.mul((A_num, A_den), self.base_ring.pow(&u, 4)),
            B: self.base_ring.mul((B_num, B_den), self.base_ring.pow(&u, 6))
        },
            move |P| match P {
                EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(
                    base_ring.mul(base_ring.pow(&u, 2), x),
                    base_ring.mul(base_ring.pow(&u, 3), y),
                ),
                EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity
            },
            move |P| match P {
                EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(
                    base_ring_copy.mul(base_ring_copy.pow(&u_inv, 2), x),
                    base_ring_copy.mul(base_ring_copy.pow(&u_inv, 3), y),
                ),
                EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity
            },
        )
    }

    pub fn torsion_group(&self) {
    }
}

#[cfg(test)]
use super::super::primitive::*;

#[test]
fn test_isomorphic_curve_over_z() {
    let Q = QType::new(BigInt::RING);
    let i = Q.clone().z_embedding(i64::RING);
    let A = Q.div(Q.pow(&i(3), 5), &Q.pow(&i(2), 4));
    let B = Q.div(Q.pow(&i(3), 6), &Q.pow(&i(2), 6));
    let curve = EllipticCurve::new(Q.clone(), A, B);
    let expected_curve = EllipticCurve::new(Q.clone(), i(3), i(1));
    let (actual_curve, f, finv) = curve.isomorphic_curve_over_z();
    assert_eq!(expected_curve, actual_curve);
    let P = EllipticCurvePoint::Affine(i(0), i(1));
    assert!(actual_curve.is_on_curve(&P));
    assert!(curve.is_on_curve(&finv(P.clone())));
    assert!(actual_curve.points_eq(&P, &f(finv(P.clone()))));
}