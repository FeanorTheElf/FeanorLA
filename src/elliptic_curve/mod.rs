#![allow(non_snake_case)]

pub mod rational_torsion_group;

use super::ring::*;
use super::bigint::*;
use super::embedding::*;
use super::integer::*;
use super::wrapper::*;
use super::la::mat::*;
use super::fq::*;
use super::combinatorics::iters::*;
use super::ring_extension::*;
use super::poly::*;
use super::fraction_field::*;

#[derive(Debug, Clone)]
pub struct EllipticCurve<K: Ring> {
    base_field: K,
    A: K::El,
    B: K::El
}

impl<K> PartialEq for EllipticCurve<K> 
    where K: CanonicalIsomorphismInfo<K>
{

    fn eq(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_field().eq(&self.A, &rhs.A) && self.base_field().eq(&self.B, &rhs.B)
    }
}

impl<K> Eq for EllipticCurve<K>
    where K: CanonicalIsomorphismInfo<K>
{}

impl<K> std::fmt::Display for EllipticCurve<K>
    where K: CanonicalIsomorphismInfo<K>
{

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Elliptic curve defined by Y^2 = X^3 + {} X + {}", self.base_field().display(&self.A), self.base_field().display(&self.B))
    }
}

pub type CoordRing<K: Ring> = SimpleRingExtension<PolyRing<K>, VectorArray<El<PolyRing<K>>, 2>, VectorArray<El<PolyRing<K>>, 2>>;
pub type FunctionField<K: Ring> = FieldOfFractions<CoordRing<K>>;

#[derive(Clone)]
pub enum EllipticCurvePoint<K> 
    where K: Ring
{
    Affine(K::El, K::El), Infinity
}

impl<K> PartialEq for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{
    fn eq(&self, other: &EllipticCurvePoint<WrappingRing<K>>) -> bool {
        match (self, other) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => true,
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) => x1 == x2 && y1 == y2,
            _ => false
        }
    }
}

impl<K> std::fmt::Debug for EllipticCurvePoint<K>
    where K: Ring, K::El: std::fmt::Debug
{
    default fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EllipticCurvePoint::Infinity => write!(f, "o̲"),
            EllipticCurvePoint::Affine(x, y) => write!(f, "({:?}, {:?})", x, y)
        }
    }
}

impl<K> std::fmt::Debug for EllipticCurvePoint<K>
    where K: Ring, K::El: std::fmt::Display + std::fmt::Debug
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EllipticCurvePoint::Infinity => write!(f, "o̲"),
            EllipticCurvePoint::Affine(x, y) => write!(f, "({}, {})", x, y)
        }
    }
}

impl<K> std::fmt::Display for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring, K::El: std::fmt::Display + std::fmt::Debug
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, f)
    }
}

impl<K> Eq for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{}

impl<K> EllipticCurve<K> 
    where K: Ring + CanonicalIsomorphismInfo<K>
{
    pub fn new(base_field: K, A: K::El, B: K::El) -> Self {
        assert!(base_field.is_field().can_use());
        let result = EllipticCurve {
            base_field, A, B
        };
        assert!(!result.base_field().is_zero(&result.discriminant()));
        return result;
    }

    pub fn lift_homomorphism<L, F>(&self, _other: &EllipticCurve<L>, f: F) -> impl Fn(EllipticCurvePoint<K>) -> EllipticCurvePoint<L>
        where L: Ring, F: Fn(K::El) -> L::El
    {
        move |x| match x {
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
            EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(f(x), f(y))
        }
    }

    ///
    /// Returns the coordinate ring together with the `x` and `y` elements of
    /// this elliptic curve.
    /// 
    pub fn coordinate_ring(&self) -> (CoordRing<K>, El<CoordRing<K>>, El<CoordRing<K>>) {
        let poly_ring = PolyRing::adjoint(self.base_field().clone(), "X");
        let mipo = Vector::new(VectorArray::new([poly_ring.zero(), 
            poly_ring.add(
                poly_ring.pow(&poly_ring.unknown(), 3),
                poly_ring.add(
                    poly_ring.mul(poly_ring.from(self.A.clone()), poly_ring.unknown()),
                    poly_ring.from(self.B.clone())
                )
            )
        ]));
        let poly_ring_x = poly_ring.unknown();
        let result = SimpleRingExtension::new(poly_ring, mipo);
        let X = result.from(poly_ring_x);
        let Y = result.generator();
        return (result, X, Y);
    }

    pub fn function_field(&self) -> (FunctionField<K>, El<FunctionField<K>>, El<FunctionField<K>>) {
        assert!(self.base_field().is_field().can_use());
        let (ring, x, y) = self.coordinate_ring();
        let field = FieldOfFractions::new(ring.clone());
        let (x, y) = {
            let incl = embedding(&ring, &field);
            (incl(x), incl(y))
        };
        return (field, x, y);
    }

    pub fn discriminant(&self) -> K::El {
        self.base_field().add(
            self.base_field().mul(self.base_field().from_z(27), self.base_field().pow(&self.B, 2)),
            self.base_field().mul(self.base_field().from_z(4), self.base_field().pow(&self.A, 3))
        )
    }

    pub fn j_invariant(&self) -> K::El {
        assert!(self.base_field().is_field().can_use());
        let A_cubed = self.base_field().pow(&self.A, 3);
        return self.base_field().div(
            self.base_field().mul(self.base_field().from_z(1728 * 4), A_cubed.clone()), 
            &self.discriminant()
        );
    }

    pub fn is_on_curve(&self, point: &EllipticCurvePoint<K>) -> bool {
        match point {
            EllipticCurvePoint::Infinity => true,
            EllipticCurvePoint::Affine(x, y) => self.base_field().eq(
                &self.base_field().pow(y, 2),
                &self.base_field().add(
                    self.base_field().pow(x, 3), 
                    self.base_field().add(
                        self.base_field().mul_ref(&self.A, x),
                        self.B.clone()
                    ))
            )
        }
    }

    ///
    /// Checks whether this curve is isomorphic to rhs over an algebraic extension of the
    /// base field. Note that this does not imply that they are isomorphic over the base field.
    /// 
    pub fn is_isomorphic(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_field().eq(&self.j_invariant(), &rhs.j_invariant())
    }

    pub fn point_add<F>(&self, a: EllipticCurvePoint<F>, b: EllipticCurvePoint<F>, field: &F) -> EllipticCurvePoint<F>
        where F: Ring + CanonicalEmbeddingInfo<K>
    {
        assert!(field.is_field().can_use());
        assert!(field.has_embedding(self.base_field()).can_use());
        match (a, b) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => EllipticCurvePoint::Infinity,
            (EllipticCurvePoint::Affine(x, y), EllipticCurvePoint::Infinity) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Affine(x, y)) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if field.eq(&x1, &x2) && field.eq(&y1, &field.neg(y2.clone())) => {
                EllipticCurvePoint::Infinity
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if field.eq(&x1, &x2) => {
                let lambda = field.div(
                    field.add(
                        field.mul(field.from_z(3), field.mul_ref(&x1, &x1)),
                        field.embed(self.base_field(), self.A.clone())
                    ),
                    &field.mul(field.from_z(2), y1)
                );
                let x = field.add(
                    field.neg(field.add_ref(x1, &x2)), 
                    field.mul_ref(&lambda, &lambda)
                );
                EllipticCurvePoint::Affine(
                    x.clone(),
                    field.add(
                        field.neg(y2), 
                        field.mul(lambda, field.sub(x2, x))
                    )
                )
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if !field.eq(&x1, &x2) => {
                let lambda = field.div(
                    field.sub_ref_fst(&y1, y2),
                    &field.sub_ref_fst(&x1, x2.clone())
                );
                let x = field.add(
                    field.neg(field.add_ref(x2, &x1)), 
                    field.mul_ref(&lambda, &lambda)
                );
                EllipticCurvePoint::Affine(
                    x.clone(),
                    field.add(
                        field.neg(y1), 
                        field.mul(lambda, field.sub(x1, x))
                    )
                )
            },
            (EllipticCurvePoint::Affine(_, _), EllipticCurvePoint::Affine(_, _)) => EllipticCurvePoint::Infinity
        }
    }

    pub fn mul_point<F>(&self, point: &EllipticCurvePoint<F>, n: &BigInt, field: &F) -> EllipticCurvePoint<F> 
        where F: Ring + CanonicalEmbeddingInfo<K>
    {
        assert!(field.is_field().can_use());
        assert!(field.has_embedding(self.base_field()).can_use());
        assert!(*n >= 0);
        if n.is_zero() {
            return EllipticCurvePoint::Infinity;
        }

        let mut result = EllipticCurvePoint::Infinity;
        for i in (0..(BigInt::RING.abs_log2_floor(n) + 1)).rev() {
            if BigInt::RING.abs_is_bit_set(n, i) {
                result = self.point_add(self.point_add(result.clone(), point.clone(), field), result, field);
            } else {
                result = self.point_add(result.clone(), result, field);
            }
        }
        return result;
    }

    pub fn inv_point(&self, point: EllipticCurvePoint<K>) -> EllipticCurvePoint<K> {
        match point {
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
            EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(x, self.base_field().neg(y))
        }
    }

    pub fn extend_ring<'a, S, F>(self, new_ring: S, incl: F) -> (EllipticCurve<S>, impl Fn(EllipticCurvePoint<K>) -> EllipticCurvePoint<S>)
        where S: Ring, F: Fn(K::El) -> S::El
    {
        (EllipticCurve {
            base_field: new_ring,
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
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) => self.base_field().eq(x1, x2) && self.base_field().eq(y1, y2),
            _ => false
        }
    }

    pub fn base_field(&self) -> &K {
        &self.base_field
    }

    pub fn a4(&self) -> &K::El {
        &self.A
    }

    pub fn a6(&self) -> &K::El {
        &self.B
    }
}

impl<K> EllipticCurve<K>
    where K: FiniteRing + CanonicalIsomorphismInfo<K>
{

    pub fn points<'a>(&'a self) -> impl 'a + Iterator<Item = EllipticCurvePoint<K>> {
        cartesian_product(elements(self.base_field()), elements(self.base_field()))
            .map(|(x, y)| EllipticCurvePoint::Affine(x, y))
            .filter(move |p| self.is_on_curve(p))
            .chain(std::iter::once(EllipticCurvePoint::Infinity))
    }
}

impl<QType> std::hash::Hash for EllipticCurvePoint<QType> 
    where QType: HashableElRing + SingletonRing
{

    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            EllipticCurvePoint::Infinity => {
                true.hash(state);
            },
            EllipticCurvePoint::Affine(x, y) => {
                let Q = QType::singleton();
                false.hash(state);
                Q.hash(state, x);
                Q.hash(state, y);
            }
        }
    }
}

#[test]
fn test_elliptic_curve_point_eq() {
    let Q = FieldOfFractions::<BigIntRing>::singleton().bind_ring_by_value();
    let i = z_hom(&Q);
    let E = EllipticCurve::new(Q.clone(), i(0), i(1));
    let P = EllipticCurvePoint::Affine(i(0), i(1));
    let Q = EllipticCurvePoint::Affine(i(0), i(-1));
    assert!(E.is_on_curve(&P));
    assert!(E.is_on_curve(&Q));
    assert_eq!(P, P);
    assert_ne!(P, Q);
}
