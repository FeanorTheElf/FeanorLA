#![allow(non_snake_case)]

pub mod finite_field_curve;
pub mod rational_torsion_group;
pub mod point_count;
// pub mod division_polynomials;

use super::ring::*;
use super::integer::*;
use super::embedding::*;
use super::wrapper::*;
use super::la::mat::*;
use super::fq::*;
use super::combinatorics::iters::*;
use super::ring_extension::simple_extension::*;
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
        self.base_field.eq(&self.A, &rhs.A) && self.base_field.eq(&self.B, &rhs.B)
    }
}

impl<K> Eq for EllipticCurve<K>
    where K: CanonicalIsomorphismInfo<K>
{}

impl<K> std::fmt::Display for EllipticCurve<K>
    where K: CanonicalIsomorphismInfo<K>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Elliptic curve defined by Y^2 = X^3 + {} X + {}", self.base_field.display(&self.A), self.base_field.display(&self.B))
    }
}

pub type CoordRing<K: Ring> = WrappingRing<SimpleRingExtension<PolyRing<K>, VectorArray<El<PolyRing<K>>, 2>, VectorArray<El<PolyRing<K>>, 2>>>;
pub type FunctionField<K: Ring> = WrappingRing<FieldOfFractions<SimpleRingExtension<PolyRing<K>, VectorArray<El<PolyRing<K>>, 2>, VectorArray<El<PolyRing<K>>, 2>>>>;

#[derive(Clone)]
pub enum EllipticCurvePoint<K> 
    where K: Ring
{
    Affine(K::El, K::El), Infinity
}

impl<K> EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{
    pub fn x(&self) -> Option<&El<WrappingRing<K>>> {
        match self {
            EllipticCurvePoint::Affine(x, y) => Some(x),
            EllipticCurvePoint::Infinity => None
        }
    }

    pub fn y(&self) -> Option<&El<WrappingRing<K>>> {
        match self {
            EllipticCurvePoint::Affine(x, y) => Some(y),
            EllipticCurvePoint::Infinity => None
        }
    }
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

impl<K> EllipticCurve<WrappingRing<K>> 
    where K: Ring
{
    pub fn new(base_field: WrappingRing<K>, A: El<WrappingRing<K>>, B: El<WrappingRing<K>>) -> Self {
        assert!(base_field.is_field().can_use());
        let result = EllipticCurve {
            base_field, A, B
        };
        assert!(!result.base_field().is_zero(&result.discriminant()));
        return result;
    }

    pub fn lift_homomorphism<L, F>(&self, _other: &EllipticCurve<WrappingRing<L>>, f: F) -> impl Fn(EllipticCurvePoint<WrappingRing<K>>) -> EllipticCurvePoint<WrappingRing<L>>
        where L: Ring, F: Fn(El<WrappingRing<K>>) -> El<WrappingRing<L>>
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
        let poly_ring = PolyRing::adjoint(self.base_field().wrapped_ring().clone(), "X");
        let mipo = Vector::new(VectorArray::new([poly_ring.zero(), 
            poly_ring.add(
                poly_ring.pow(&poly_ring.unknown(), 3),
                poly_ring.add(
                    poly_ring.mul(poly_ring.from(self.A.clone().into_val()), poly_ring.unknown()),
                    poly_ring.from(self.B.clone().into_val())
                )
            )
        ]));
        let poly_ring_x = poly_ring.unknown();
        let result = SimpleRingExtension::new(poly_ring, mipo).bind_ring_by_value();
        let X = result.from(result.wrapped_ring().from(poly_ring_x));
        let Y = result.from(result.wrapped_ring().generator());
        return (result, X, Y);
    }

    pub fn function_field(&self) -> (FunctionField<K>, El<FunctionField<K>>, El<FunctionField<K>>) {
        assert!(self.base_field().is_field().can_use());
        let (ring, x, y) = self.coordinate_ring();
        let field = FieldOfFractions::new(ring.wrapped_ring().clone()).bind_ring_by_value();
        let (x, y) = {
            let incl = embedding(&ring, &field);
            (incl(x), incl(y))
        };
        return (field, x, y);
    }

    pub fn discriminant(&self) -> El<WrappingRing<K>> {
        self.B.pow(2) * 27 + self.A.pow(3) * 4
    }

    pub fn j_invariant(&self) -> El<WrappingRing<K>> {
        assert!(self.base_field().is_field().can_use());
        self.A.pow(3) * 1728 * 4 / self.discriminant()
    }

    pub fn is_on_curve(&self, point: &EllipticCurvePoint<WrappingRing<K>>) -> bool {
        match point {
            EllipticCurvePoint::Infinity => true,
            EllipticCurvePoint::Affine(x, y) => y.pow(2) == x.pow(3) + x * &self.A + &self.B
        }
    }

    ///
    /// Checks whether this curve is isomorphic to rhs over an algebraic extension of the
    /// base field. Note that this does not imply that they are isomorphic over the base field.
    /// 
    pub fn is_isomorphic(&self, rhs: &EllipticCurve<WrappingRing<K>>) -> bool {
        self.j_invariant() == rhs.j_invariant()
    }

    pub fn point_add<F>(&self, a: EllipticCurvePoint<WrappingRing<F>>, b: EllipticCurvePoint<WrappingRing<F>>, field: &WrappingRing<F>) -> EllipticCurvePoint<WrappingRing<F>>
        where F: Ring + CanonicalEmbeddingInfo<K>
    {
        assert!(field.is_field().can_use());
        assert!(field.has_embedding(self.base_field()).can_use());
        match (a, b) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => EllipticCurvePoint::Infinity,
            (EllipticCurvePoint::Affine(x, y), EllipticCurvePoint::Infinity) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Affine(x, y)) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if x1 == x2 && y1 == -y2.clone() => {
                EllipticCurvePoint::Infinity
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if x1 == x2 => {
                let lambda = (x1.pow(2) * 3 + field.embed(self.base_field(), self.A.clone())) / 2 / y1;
                let x = lambda.pow(2) - x1 - &x2;
                EllipticCurvePoint::Affine(
                    x.clone(),
                    (x2 - x) * lambda - y2
                )
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if !field.eq(&x1, &x2) => {
                let lambda = (y1 - &y2) / (&x1 - &x2);
                let x = lambda.pow(2) - x1 - &x2;
                EllipticCurvePoint::Affine(
                    x.clone(),
                    (x2 - x) * lambda - y2
                )
            },
            (EllipticCurvePoint::Affine(_, _), EllipticCurvePoint::Affine(_, _)) => EllipticCurvePoint::Infinity
        }
    }

    pub fn mul_point<F>(&self, point: &EllipticCurvePoint<WrappingRing<F>>, n: &BigInt, field: &WrappingRing<F>) -> EllipticCurvePoint<WrappingRing<F>> 
        where F: Ring + CanonicalEmbeddingInfo<K>
    {
        assert!(field.is_field().can_use());
        assert!(field.has_embedding(self.base_field()).can_use());
        assert!(*n >= 0);
        if n.is_zero() {
            return EllipticCurvePoint::Infinity;
        }

        let mut result = EllipticCurvePoint::Infinity;
        for i in (0..(n.abs_log2_floor() + 1)).rev() {
            if n.is_bit_set(i) {
                result = self.point_add(self.point_add(result.clone(), point.clone(), field), result, field);
            } else {
                result = self.point_add(result.clone(), result, field);
            }
        }
        return result;
    }

    pub fn inv_point(&self, point: EllipticCurvePoint<WrappingRing<K>>) -> EllipticCurvePoint<WrappingRing<K>> {
        match point {
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
            EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(x, -y)
        }
    }

    pub fn extend_ring<'a, S, F>(self, new_ring: WrappingRing<S>, incl: F) -> (EllipticCurve<WrappingRing<S>>, impl Fn(EllipticCurvePoint<WrappingRing<K>>) -> EllipticCurvePoint<WrappingRing<S>>)
        where S: Ring, F: Fn(El<WrappingRing<K>>) -> El<WrappingRing<S>>
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

    pub fn points_eq(&self, P: &EllipticCurvePoint<WrappingRing<K>>, Q: &EllipticCurvePoint<WrappingRing<K>>) -> bool {
        match (P, Q) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => true,
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) => x1 == x2 && y1 == y2,
            _ => false
        }
    }

    pub fn base_field(&self) -> &WrappingRing<K> {
        &self.base_field
    }

    pub fn a4(&self) -> &El<WrappingRing<K>> {
        &self.A
    }

    pub fn a6(&self) -> &El<WrappingRing<K>> {
        &self.B
    }
}

impl<K> EllipticCurve<WrappingRing<K>>
    where K: FiniteRing + CanonicalIsomorphismInfo<K>
{

    pub fn points<'a>(&'a self) -> impl 'a + Iterator<Item = EllipticCurvePoint<WrappingRing<K>>> {
        cartesian_product(elements(self.base_field()), elements(self.base_field()))
            .map(|(x, y)| EllipticCurvePoint::Affine(x, y))
            .filter(move |p| self.is_on_curve(p))
            .chain(std::iter::once(EllipticCurvePoint::Infinity))
    }
}

impl<K> std::hash::Hash for EllipticCurvePoint<WrappingRing<K>> 
    where K: HashableElRing
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            EllipticCurvePoint::Infinity => {
                true.hash(state);
            },
            EllipticCurvePoint::Affine(x, y) => {
                false.hash(state);
                x.hash(state);
                y.hash(state);
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
