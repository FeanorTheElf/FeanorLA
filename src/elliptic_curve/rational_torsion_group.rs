use super::*;
use super::super::eea::*;
use super::super::rational::*;
use super::super::integer::*;

use std::collections::HashSet;
use std::cmp::Ordering;

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
    where QType: RationalField + HashableElRing + SingletonRing, 
        QType::BaseRing: SingletonRing + UfdInfoRing + EuclideanInfoRing
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
        let (A_num, A_den) = (self.A.num().clone_ring(), self.A.den().clone_ring());
        let (B_num, B_den) = (self.B.num().clone_ring(), self.B.den().clone_ring());

        let Z = self.base_field().base_ring().clone_ring();
        let Q = self.base_field().clone();
        let i = self.base_field().embedding();

        // there seems to be no better way to find the correct u without factoring A resp. B
        // as even for checking whether a number is square-free currently there is no better
        // method known
        let u_num_A = Z.factor(A_den.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow((power as u32 - 1) / 4 + 1))
            .product::<RingElWrapper<QType::BaseRing>>();
        let u_num_B = Z.factor(B_den.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow((power as u32 - 1) / 6 + 1))
            .product::<RingElWrapper<QType::BaseRing>>();
        
        let u_den_A = Z.factor(A_num.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow(power as u32 / 4))
            .product::<RingElWrapper<QType::BaseRing>>();
        let u_den_B = Z.factor(B_num.clone()).into_iter()
            .map(|(factor, power)| factor.into_val().pow(power as u32 / 6))
            .product::<RingElWrapper<QType::BaseRing>>();

        let u: RingElWrapper<QType> = i(lcm(&Z, u_num_A, u_num_B)) / i(gcd(&Z, u_den_A, u_den_B));
        let u_inv: RingElWrapper<QType> = Q.one() / &u;
        return (EllipticCurve {
            base_field: Q.clone(),
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
        let Z = self.base_field().base_ring().clone_ring();
        let Q = self.base_field().clone();
        let i = self.base_field().embedding();

        let (E, _f, finv) = self.isomorphic_curve_over_z();
        let disc = E.discriminant();
        // it is a theorem that for all torsion points (x, y), have y^2 | Disc
        let disc = (disc.num() / disc.den()).clone_ring();
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
        let A = (E.A.num() / E.A.den()).clone_ring();
        let B = (E.B.num() / E.B.den()).clone_ring();

        let mut result: HashSet<EllipticCurvePoint<WrappingRing<QType>>> = HashSet::new();
        result.insert(EllipticCurvePoint::Infinity);
        for y in possible_y {
            let B_minus_y2 = &B - &y * &y;
            for x in IntegralCubic::new(A.val(), B_minus_y2.val(), y.parent_ring()).find_integral_roots() {
                let x = Z.from(x);
                let point = EllipticCurvePoint::Affine(i(x), i(y.clone()));
                // it is a theorem that the torsion group has order dividing 24
                if E.points_eq(&E.mul_point::<QType>(&point, &BigInt::from(24), &Q), &EllipticCurvePoint::Infinity) {
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
    let Q = FractionFieldImpl::<BigIntRing>::singleton().bind_ring_by_value();
    let i = z_hom(&Q);
    let A = i(3).pow(5) / i(2).pow(4);
    let B = i(3).pow(6) / i(2).pow(6);
    let curve = EllipticCurve::new(Q.clone(), A, B);
    let expected_curve = EllipticCurve::new(Q.clone(), i(3), i(1));
    let (actual_curve, f, finv) = curve.isomorphic_curve_over_z();
    assert_eq!(expected_curve, actual_curve);
    let P = EllipticCurvePoint::Affine(i(0), i(1));
    assert!(actual_curve.is_on_curve(&P, &Q));
    assert!(curve.is_on_curve(&finv(P.clone()), &Q));
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
    assert_eq!(
        vec![i(1)],
        IntegralCubic { p: &i(-2), q: &i(1), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
    assert_eq!(
        vec![i(-3), i(1), i(2)],
        IntegralCubic { p: &i(-7), q: &i(6), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
}

#[test]
fn test_compute_torsion_group() {
    let Q = FractionFieldImpl::<BigIntRing>::singleton().bind_ring_by_value();
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