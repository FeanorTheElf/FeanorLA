#![allow(non_snake_case)]
use super::super::ring::*;
use super::super::bigint::*;
use super::super::embedding::*;
use super::super::wrapper::*;
use super::super::la::mat::*;
use super::fq::*;
use super::super::combinatorics::iters::*;
use super::ring_ext::*;
use super::poly::*;
use super::fractions::*;
use super::eea::*;
use super::rationals::*;
use super::integer::*;

use std::collections::HashSet;
use std::cmp::Ordering;

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

impl<K> EllipticCurve<K>
    where K: FiniteRing + CanonicalIsomorphismInfo<K>
{

    pub fn points<'a>(&'a self) -> impl 'a + Iterator<Item = EllipticCurvePoint<K>> {
        cartesian_product(elements(&self.base_ring), elements(&self.base_ring))
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

struct IntegralCubic<'a, R: IntegerRing + OrderedRing> {
    p: &'a R::El,
    q: &'a R::El,
    ring: R
}

impl<'a, R: IntegerRing + OrderedRing + EuclideanInfoRing> IntegralCubic<'a, R> {

    fn new(p: &'a R::El, q: &'a R::El, ring: R) -> Self {
        IntegralCubic {
            p, q, ring
        }
    }

    fn eval(&self, x: &R::El) -> R::El {
        self.ring.add_ref(self.ring.add(self.ring.pow(x, 3), self.ring.mul_ref(x, &self.p)), &self.q)
    }

    ///
    /// Calls the given closure on all integral roots of the cubic.
    /// The closure is only called once, even if the root is a multiple roots.
    /// 
    /// The algorithm relies on bisection, and works reliably with very large numbers.
    /// 
    fn calc_integral_roots<F>(&self, mut f: F)
        where F: FnMut(R::El)
    {
        let Z = &self.ring;
        let mut process_potential_root = |x: R::El| if Z.is_zero(&self.eval(&x)) { f(x) };
        let diff_disc = Z.mul(Z.from_z(-12), self.p.clone());
        if Z.cmp(&diff_disc, &Z.zero()) != Ordering::Greater {
            // zero or one maximum/minimum, so there can be at most one root
            process_potential_root(Z.find_zero_floor(|x| self.eval(x), Z.zero()));
        } else {
            let root_size_bound = if Z.abs_cmp(&self.p, &self.q) == std::cmp::Ordering::Less {
                Z.abs(self.q.clone())
            } else {
                Z.abs(self.p.clone())
            };
            let extremum_floor = Z.floor_div(Z.root_floor(&diff_disc, 2), &Z.from_z(6));
            // on the intervals [a0, a1], [b0, b1], [c0, c1], the function is monotonous, 
            // hence has at most one root
            let a0 = Z.neg(root_size_bound.clone());
            let a1 = Z.neg(extremum_floor.clone());
            let b0 = Z.sub_ref_fst(&a1, Z.one());
            let b1 = extremum_floor;
            let c0 = Z.add_ref(Z.one(), &b1);
            let c1 = root_size_bound;
            let a0_val = self.eval(&a0);
            let a1_val = self.eval(&a1);
            let b0_val = self.eval(&b0);
            let b1_val = self.eval(&b1);
            let c0_val = self.eval(&c0);
            let c1_val = self.eval(&c1);
            let zero = Z.zero();
            if Z.cmp(&a0_val, &zero) != Ordering::Greater && Z.cmp(&a1_val, &zero) != Ordering::Less {
                process_potential_root(Z.bisect(|x| self.eval(&x), a0, a1));
            }
            if Z.cmp(&b0_val, &zero) != Ordering::Less && Z.cmp(&b1_val, &zero) != Ordering::Greater {
                process_potential_root(Z.bisect(|x| Z.neg(self.eval(&x)), b0, b1));
            }
            if Z.cmp(&c0_val, &zero) != Ordering::Greater && Z.cmp(&c1_val, &zero) != Ordering::Less {
                process_potential_root(Z.bisect(|x| self.eval(&x), c0, c1));
            }
        }
    }

    fn find_integral_roots(&self) -> impl Iterator<Item = R::El> {
        let mut roots = [None, None, None];
        let mut i = 0;
        self.calc_integral_roots(|root| {
            roots[i] = Some(root);
            i += 1;
        });
        return <[Option<_>; 3] as std::iter::IntoIterator>::into_iter(roots).filter_map(|x| x);
    }
}

impl<QType> EllipticCurve<WrappingRing<QType>>
    where QType: RationalField + HashableElRing + SingletonRing + CanonicalIsomorphismInfo<QType>, 
        QType::UnderlyingIntegers: SingletonRing + UfdInfoRing + EuclideanInfoRing + OrderedRing
{

    ///
    /// Returns an isomorphic elliptic curve E': y^2 = x^3 + Ax + B with A, B in Z
    /// with minimal discriminant (in absolute value) and isomorphisms f: E -> E' 
    /// and f^-1: E' -> E.
    /// 
    pub fn isomorphic_curve_over_z(&self) -> (
        Self, 
        impl Fn(EllipticCurvePoint<WrappingRing<QType>>) -> EllipticCurvePoint<WrappingRing<QType>>, 
        impl Fn(EllipticCurvePoint<WrappingRing<QType>>) -> EllipticCurvePoint<WrappingRing<QType>>
    ) {
        let (A_num, A_den) = (self.base_ring.num(&self.A), self.base_ring.den(&self.A));
        let (B_num, B_den) = (self.base_ring.num(&self.B), self.base_ring.den(&self.B));

        let Z = self.base_ring.underlying_integers();
        let Q = self.base_ring.clone();

        // there seems to be no better way to find the correct u without factoring A resp. B
        // as even for checking whether a number is square-free currently there is no better
        // method known
        let u_num_A = Z.factor(A_den.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow((power as u32 - 1) / 4 + 1))
            .product::<RingElWrapper<QType::UnderlyingIntegers>>();
        let u_num_B = Z.factor(B_den.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow((power as u32 - 1) / 6 + 1))
            .product::<RingElWrapper<QType::UnderlyingIntegers>>();
        
        let u_den_A = Z.factor(A_num.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow(power as u32 / 4))
            .product::<RingElWrapper<QType::UnderlyingIntegers>>();
        let u_den_B = Z.factor(B_num.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow(power as u32 / 6))
            .product::<RingElWrapper<QType::UnderlyingIntegers>>();

        let u: RingElWrapper<QType> = Q.embed(&Z, lcm(&Z, u_num_A, u_num_B)) / Q.embed(&Z, gcd(&Z, u_den_A, u_den_B));
        let u_inv: RingElWrapper<QType> = Q.one() / &u;
        return (EllipticCurve {
            base_ring: Q.clone(),
            A: u.pow(4) * &self.A,
            B: u.pow(6) * &self.B
        },
            move |P: EllipticCurvePoint<WrappingRing<QType>>| match P {
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

    pub fn torsion_group(&self) -> HashSet<EllipticCurvePoint<WrappingRing<QType>>> {
        let Z = self.base_ring.underlying_integers();
        let Q = self.base_ring.clone();

        let (E, _f, finv) = self.isomorphic_curve_over_z();
        let disc = E.discriminant();
        // it is a theorem that for all torsion points (x, y), have y^2 | Disc
        let disc = Z.quotient(&Q.num(&disc), &Q.den(&disc)).unwrap();
        let y_multiple = Z.factor(disc).into_iter()
            .map(|(factor, power)| (factor, power as u32 / 2));
        // note that the divisors of y_multiple are bijective to the cartesian product
        // { 0, ..., e1 } x ... x { 0, ..., en } via e -> product pi^ei
        let possible_y = std::iter::once(Z.zero()).chain(
            multi_cartesian_product(
                y_multiple.map(|(factor, power)| (0..=power).map(move |i| (factor.clone().into_val(), i))), 
                |factorization| factorization.iter()
                    .map(|(factor, power)| factor.pow(*power))
                    .product::<RingElWrapper<_>>()
            )
        );
        let A = Z.quotient(&Q.num(&E.A), &Q.den(&E.A)).unwrap();
        let B = Z.quotient(&Q.num(&E.B), &Q.den(&E.B)).unwrap();

        let mut result: HashSet<EllipticCurvePoint<WrappingRing<QType>>> = HashSet::new();
        result.insert(EllipticCurvePoint::Infinity);
        for y in possible_y {
            let B_minus_y2 = &B - &y * &y;
            for x in IntegralCubic::new(&A, &B_minus_y2, y.ring()).find_integral_roots() {
                let point = EllipticCurvePoint::Affine(Q.embed(&Z, x), Q.embed(&Z, y.clone()));
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
    let Q = FieldOfFractions::<BigIntRing>::singleton().bind_ring_by_value();
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
        IntegralCubic { p: &i(1), q: &i(1), ring: &BigInt::RING }.eval(&BigInt::ZERO)
    );
    assert_eq!(
        i(-3),
        IntegralCubic { p: &i(-2), q: &i(1), ring: &BigInt::RING }.eval(&i(-2))
    )
}

#[test]
fn test_find_integral_roots() {
    let i = z_hom(&BigInt::RING);
    assert_eq!(
        Vec::<BigInt>::new(),
        IntegralCubic { p: &i(1), q: &i(1), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
    println!("Done first");
    assert_eq!(
        vec![i(1)],
        IntegralCubic { p: &i(-2), q: &i(1), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
    println!("Done second");
    assert_eq!(
        vec![i(-3), i(1), i(2)],
        IntegralCubic { p: &i(-7), q: &i(6), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
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

#[test]
fn test_compute_torsion_group() {
    let Q = FieldOfFractions::<BigIntRing>::singleton().bind_ring_by_value();
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