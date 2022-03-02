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

impl<K: Ring> EllipticCurve<K> {

    pub fn new(base_ring: K, A: K::El, B: K::El) -> Self {
        let result = EllipticCurve {
            base_ring, A, B
        };
        assert!(!result.base_ring.is_zero(&result.discriminant()));
        return result;
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

impl std::hash::Hash for EllipticCurvePoint<QType> {

    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let Q = QType::singleton();
        match self {
            EllipticCurvePoint::Infinity => {
                true.hash(state);
            },
            EllipticCurvePoint::Affine(x, y) => {
                false.hash(state);
                let (x_num, x_den) = Q.reduce(x.clone());
                let (y_num, y_den) = Q.reduce(y.clone());
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

    fn eval(&self, x: &BigInt) -> BigInt {
        BigInt::RING.add_ref(BigInt::RING.add(BigInt::RING.pow(x, 3), BigInt::RING.mul_ref(x, &self.p)), &self.q)
    }

    fn eval_diff(&self, x: &BigInt) -> BigInt {
        BigInt::RING.add_ref(BigInt::RING.mul_ref(x, x) * 3, &self.p)
    }

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
                dbg!(&b0, &b1, &b0_val, &b1_val);
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

    pub fn torsion_group(&self) -> HashSet<EllipticCurvePoint<QType>> {
        let Z = BigInt::RING;
        let (E, _f, finv) = self.isomorphic_curve_over_z();
        let (num, den) = E.discriminant();
        let disc = Z.quotient(&num, &den).unwrap();
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
        let A = Z.quotient(&E.A.0, &E.A.1).unwrap();
        let B = Z.quotient(&E.B.0, &E.B.1).unwrap();
        for y in possible_y {
            let y2 = Z.mul_ref(&y, &y);
            let potential_x1 = BigInt::find_zero_floor(
                |x| Z.add_ref(Z.add(x.pow(3), Z.mul_ref(&x, &A)), &B), 
                BigInt::ZERO,
            );
        }
        unimplemented!() 
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

#[test]
fn test_integral_cubic_eval() {
    let i = BigInt::RING.z_embedding(i64::RING);
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
    let i = BigInt::RING.z_embedding(i64::RING);
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