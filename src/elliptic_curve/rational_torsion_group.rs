use super::*;
use super::super::eea::*;
use super::super::rational::*;
use super::super::integer::roots;

use std::collections::HashSet;

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
            for x in roots::IntegralCubic::new(A.val(), B_minus_y2.val(), y.parent_ring()).find_integral_roots() {
                let x = Z.from(x);
                let point = EllipticCurvePoint::Affine(i(x), i(y.clone()));
                // it is a theorem that the torsion group has order dividing 24
                if E.points_eq(&E.mul_point::<QType>(&point, &StdInt::from(24), &Q), &EllipticCurvePoint::Infinity) {
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
    let Q = WrappingRing::new(FractionFieldImpl::<BigIntSOORing>::singleton());
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
fn test_compute_torsion_group() {
    let Q = WrappingRing::new(FractionFieldImpl::<BigIntSOORing>::singleton());
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