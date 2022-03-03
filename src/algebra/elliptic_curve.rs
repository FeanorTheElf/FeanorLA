#![allow(non_snake_case)]
use super::super::ring::*;
use super::super::bigint::*;
use super::super::embedding::*;
use super::super::wrapper::*;
use super::super::la::mat::*;
use super::super::combinatorics::iters::*;
use super::ring_ext::*;
use super::poly::*;
use super::fractions::*;
use super::eea::*;

use std::collections::HashSet;

#[derive(Debug, Clone)]
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

impl<K: Ring> Eq for EllipticCurve<K> {}

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

impl<K> std::fmt::Debug for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EllipticCurvePoint::Infinity => write!(f, "o̲"),
            EllipticCurvePoint::Affine(x, y) => write!(f, "({}, {})", x, y)
        }
    }
}

impl<K> std::fmt::Display for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        <Self as std::fmt::Debug>::fmt(self, f)
    }
}

impl<K> Eq for EllipticCurvePoint<WrappingRing<K>>
    where K: Ring
{}

impl<K: Ring> EllipticCurve<K> {

    pub fn new(base_ring: K, A: K::El, B: K::El) -> Self {
        let result = EllipticCurve {
            base_ring, A, B
        };
        assert!(!result.base_ring.is_zero(&result.discriminant()));
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
            let incl = embedding(&ring, &field);
            (incl(x), incl(y))
        };
        return (field, x, y);
    }

    pub fn discriminant(&self) -> K::El {
        self.base_ring.add(
            self.base_ring.mul(self.base_ring.from_z(27), self.base_ring.pow(&self.B, 2)),
            self.base_ring.mul(self.base_ring.from_z(4), self.base_ring.pow(&self.A, 3))
        )
    }

    pub fn j_invariant(&self) -> K::El {
        assert!(self.base_ring.is_field().can_use());
        let A_cubed = self.base_ring.pow(&self.A, 3);
        return self.base_ring.div(
            self.base_ring.mul(self.base_ring.from_z(1728 * 4), A_cubed.clone()), 
            &self.discriminant()
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

    ///
    /// Checks whether this curve is isomorphic to rhs over an algebraic extension of the
    /// base field. Note that this does not imply that they are isomorphic over the base field.
    /// 
    pub fn is_isomorphic(&self, rhs: &EllipticCurve<K>) -> bool {
        self.base_ring.eq(&self.j_invariant(), &rhs.j_invariant())
    }

    pub fn point_add(&self, a: EllipticCurvePoint<K>, b: EllipticCurvePoint<K>) -> EllipticCurvePoint<K> {
        assert!(self.base_ring.is_field().can_use());
        match (a, b) {
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Infinity) => EllipticCurvePoint::Infinity,
            (EllipticCurvePoint::Affine(x, y), EllipticCurvePoint::Infinity) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Infinity, EllipticCurvePoint::Affine(x, y)) => EllipticCurvePoint::Affine(x, y),
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if self.base_ring.eq(&x1, &x2) && self.base_ring.eq(&y1, &self.base_ring.neg(y2.clone())) => {
                EllipticCurvePoint::Infinity
            },
            (EllipticCurvePoint::Affine(x1, y1), EllipticCurvePoint::Affine(x2, y2)) if self.base_ring.eq(&x1, &x2) => {
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

    pub fn mul_point(&self, point: &EllipticCurvePoint<K>, n: &BigInt) -> EllipticCurvePoint<K> {
        assert!(*n >= 0);
        if n.is_zero() {
            return EllipticCurvePoint::Infinity;
        }

        let mut result = EllipticCurvePoint::Infinity;
        for i in (0..(n.abs_log2_floor() + 1)).rev() {
            if n.is_bit_set(i) {
                result = self.point_add(self.point_add(result.clone(), point.clone()), result);
            } else {
                result = self.point_add(result.clone(), result);
            }
        }
        return result;
    }

    pub fn inv_point(&self, point: EllipticCurvePoint<K>) -> EllipticCurvePoint<K> {
        match point {
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
            EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(x, self.base_ring.neg(y))
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

type QType = WrappingRing<FieldOfFractions<BigIntRing>>;
type QEl = <QType as Ring>::El;

impl std::hash::Hash for EllipticCurvePoint<QType> {

    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            EllipticCurvePoint::Infinity => {
                true.hash(state);
            },
            EllipticCurvePoint::Affine(x, y) => {
                let Q = QType::singleton();
                false.hash(state);
                let (x_num, x_den) = Q.wrapped_ring().reduce(x.val().clone());
                let (y_num, y_den) = Q.wrapped_ring().reduce(y.val().clone());
                x_num.hash(state);
                x_den.hash(state);
                y_num.hash(state);
                y_den.hash(state);
            }
        }
    }
}

struct IntegralCubic<'a> {
    p: &'a BigInt,
    q: &'a BigInt
}

impl<'a> IntegralCubic<'a> {

    fn new(p: &'a BigInt, q: &'a BigInt) -> Self {
        IntegralCubic {
            p, q
        }
    }

    fn eval(&self, x: &BigInt) -> BigInt {
        BigInt::RING.add_ref(BigInt::RING.add(BigInt::RING.pow(x, 3), BigInt::RING.mul_ref(x, &self.p)), &self.q)
    }

    ///
    /// Calls the given closure on all integral roots of the cubic.
    /// The closure is only called once, even if the root is a multiple roots.
    /// 
    /// The algorithm relies on bisection, and works reliably with very large numbers.
    /// 
    fn calc_integral_roots<F>(&self, mut f: F)
        where F: FnMut(BigInt)
    {
        let mut process_potential_root = |x: BigInt| if self.eval(&x) == 0 { f(x) };
        let Z = BigInt::RING;
        let diff_disc = Z.mul(Z.from_z(-12), self.p.clone());
        if diff_disc <= 0 {
            // zero or one maximum/minimum, so there can be at most one root
            process_potential_root(BigInt::find_zero_floor(|x| self.eval(x), BigInt::ZERO));
        } else {
            let root_size_bound = if BigInt::abs_compare(&self.p, &self.q) == std::cmp::Ordering::Less {
                self.q.clone().abs()
            } else {
                self.p.clone().abs()
            };
            let extremum_floor = BigInt::root_floor(diff_disc, 2).floor_div_small(6);
            // on the intervals [a0, a1], [b0, b1], [c0, c1], the function is monotonous, 
            // hence has at most one root
            let a0 = -root_size_bound.clone();
            let a1 = -extremum_floor.clone();
            let b0 = a1.clone() - 1;
            let b1 = extremum_floor;
            let c0 = b1.clone() + 1;
            let c1 = root_size_bound;
            let a0_val = self.eval(&a0);
            let a1_val = self.eval(&a1);
            let b0_val = self.eval(&b0);
            let b1_val = self.eval(&b1);
            let c0_val = self.eval(&c0);
            let c1_val = self.eval(&c1);
            if a0_val <= 0 && a1_val >= 0 {
                process_potential_root(BigInt::bisect(|x| self.eval(&x), a0, a1));
            }
            if b0_val >= 0 && b1_val <= 0 {
                process_potential_root(BigInt::bisect(|x| -self.eval(&x), b0, b1));
            }
            if c0_val <= 0 && c1_val >= 0 {
                process_potential_root(BigInt::bisect(|x| self.eval(&x), c0, c1));
            }
        }
    }

    fn find_integral_roots(&self) -> impl Iterator<Item = BigInt> {
        let mut roots = [None, None, None];
        let mut i = 0;
        self.calc_integral_roots(|root| {
            roots[i] = Some(root);
            i += 1;
        });
        return std::array::IntoIter::new(roots).filter_map(|x| x);
    }
}

impl EllipticCurve<QType> 
{

    ///
    /// Returns an isomorphic elliptic curve E': y^2 = x^3 + Ax + B with A, B in Z
    /// with minimal discriminant and isomorphisms f: E -> E' and f^-1: E' -> E.
    /// 
    pub fn isomorphic_curve_over_z(&self) -> (Self, impl Fn(EllipticCurvePoint<QType>) -> EllipticCurvePoint<QType>, impl Fn(EllipticCurvePoint<QType>) -> EllipticCurvePoint<QType>) {
        let (A_num, A_den) = self.base_ring.wrapped_ring().reduce(self.A.val().clone());
        let (B_num, B_den) = self.base_ring.wrapped_ring().reduce(self.B.val().clone());

        let Z = BigInt::RING;
        let Q = QType::singleton();

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

        let u: QEl = Q.from_z_big(&lcm(&BigInt::RING, u_num_A, u_num_B)) / Q.from_z_big(&gcd(&BigInt::RING, u_den_A, u_den_B));
        let u_inv: QEl = Q.one() / &u;
        return (EllipticCurve {
            base_ring: self.base_ring.clone(),
            A: self.A.clone() * u.pow(4),
            B: self.B.clone() * u.pow(6)
        },
            move |P: EllipticCurvePoint<QType>| match P {
                EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(
                    u.pow(2) * x,
                    u.pow(3) * y,
                ),
                EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity
            },
            move |P| match P {
                EllipticCurvePoint::Affine(x, y) => EllipticCurvePoint::Affine(
                    u_inv.pow(2) * x,
                    u_inv.pow(3) * y,
                ),
                EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity
            },
        )
    }

    pub fn torsion_group(&self) -> HashSet<EllipticCurvePoint<QType>> {
        let Z = BigInt::RING;
        let Q = QType::singleton();
        let (E, _f, finv) = self.isomorphic_curve_over_z();
        let disc = E.discriminant();
        // it is a theorem that for all torsion points (x, y), have y^2 | Disc
        let disc = Z.quotient(&disc.val().0, &disc.val().1).unwrap();
        let y_multiple = Z.factor(disc).into_iter()
            .map(|(factor, power)| (factor, power as u32 / 2));
        // note that the divisors of y_multiple are bijective to the cartesian product
        // { 0, ..., e1 } x ... x { 0, ..., en } via e -> product pi^ei
        let possible_y = std::iter::once(BigInt::ZERO).chain(
            multi_cartesian_product(
                y_multiple.map(|(factor, power)| (0..=power).map(move |i| (factor.clone(), i))), 
                |factorization| factorization.iter()
                    .map(|(factor, power)| factor.pow(*power))
                    .product::<RingElWrapper<_>>()
                    .into_val().abs()
            )
        );
        let A = Z.quotient(&E.A.val().0, &E.A.val().1).unwrap();
        let B = Z.quotient(&E.B.val().0, &E.B.val().1).unwrap();

        let mut result: HashSet<EllipticCurvePoint<QType>> = HashSet::new();
        result.insert(EllipticCurvePoint::Infinity);
        for y in possible_y {
            let y2 = Z.mul_ref(&y, &y);
            let B_minus_y2 = Z.sub_ref_fst(&B, y2);
            for x in IntegralCubic::new(&A, &B_minus_y2).find_integral_roots() {
                let point = EllipticCurvePoint::Affine(Q.from_z_big(&x), Q.from_z_big(&y));
                // it is a theorem that the torsion group has order dividing 24
                if E.points_eq(&E.mul_point(&point, &BigInt::from(24)), &EllipticCurvePoint::Infinity) {
                    result.insert(finv(point.clone()));
                    result.insert(finv(E.inv_point(point)));
                }
            }
        }
        assert!(24 % result.len() == 0);
        return result;
    }
}

#[test]
fn test_isomorphic_curve_over_z() {
    let Q = QType::singleton();
    let i = z_hom(&Q);
    let A = i(3).pow(5) / i(2).pow(4);
    let B = i(3).pow(6) / i(2).pow(6);
    let curve = EllipticCurve::new(Q.clone(), A, B);
    let expected_curve = EllipticCurve::new(Q.clone(), i(3), i(1));
    let (actual_curve, f, finv) = curve.isomorphic_curve_over_z();
    assert_eq!(expected_curve, actual_curve);
    let P = EllipticCurvePoint::Affine(i(0), i(1));
    assert!(actual_curve.is_on_curve(&P));
    assert!(curve.is_on_curve(&finv(P.clone())));
    assert!(actual_curve.points_eq(&P, &f(finv(P.clone()))));
}

#[test]
fn test_integral_cubic_eval() {
    let i = z_hom(&BigInt::RING);
    assert_eq!(
        i(1),
        IntegralCubic { p: &i(1), q: &i(1) }.eval(&BigInt::ZERO)
    );
    assert_eq!(
        i(-3),
        IntegralCubic { p: &i(-2), q: &i(1) }.eval(&i(-2))
    )
}

#[test]
fn test_find_integral_roots() {
    let i = z_hom(&BigInt::RING);
    assert_eq!(
        Vec::<BigInt>::new(),
        IntegralCubic { p: &i(1), q: &i(1) }.find_integral_roots().collect::<Vec<_>>()
    );
    assert_eq!(
        vec![i(1)],
        IntegralCubic { p: &i(-2), q: &i(1) }.find_integral_roots().collect::<Vec<_>>()
    );
    assert_eq!(
        vec![i(-3), i(1), i(2)],
        IntegralCubic { p: &i(-7), q: &i(6) }.find_integral_roots().collect::<Vec<_>>()
    );
}

#[test]
fn test_elliptic_curve_point_eq() {
    let Q = QType::singleton();
    let i = z_hom(&Q);
    let E = EllipticCurve::new(Q.clone(), i(0), i(1));
    let P = EllipticCurvePoint::Affine(i(0), i(1));
    let Q = EllipticCurvePoint::Affine(i(0), i(-1));
    assert!(E.is_on_curve(&P));
    assert!(E.is_on_curve(&Q));
    assert_eq!(P, P);
    assert_ne!(P, Q);
}

#[test]
fn test_compute_torsion_group() {
    let Q = QType::singleton();
    let i = z_hom(&Q);
    let E = EllipticCurve::new(Q.clone(), i(0), i(3));
    let mut torsion_group = HashSet::new();
    torsion_group.insert(EllipticCurvePoint::Infinity);
    assert_eq!(torsion_group, E.torsion_group());

    let E = EllipticCurve::new(Q.clone(), i(1), i(0));
    let mut torsion_group = HashSet::new();
    torsion_group.insert(EllipticCurvePoint::Infinity);
    torsion_group.insert(EllipticCurvePoint::Affine(i(0), i(0)));
    assert_eq!(torsion_group, E.torsion_group());

    let E = EllipticCurve::new(Q.clone(), i(0), i(1));
    let mut torsion_group = HashSet::new();
    torsion_group.insert(EllipticCurvePoint::Infinity);
    torsion_group.insert(EllipticCurvePoint::Affine(i(-1), i(0)));
    torsion_group.insert(EllipticCurvePoint::Affine(i(2), i(3)));
    torsion_group.insert(EllipticCurvePoint::Affine(i(2), i(-3)));
    torsion_group.insert(EllipticCurvePoint::Affine(i(0), i(1)));
    torsion_group.insert(EllipticCurvePoint::Affine(i(0), i(-1)));
    assert_eq!(torsion_group, E.torsion_group());
}