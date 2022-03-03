use super::super::super::ring::*;
use super::super::super::embedding::*;
use super::uni_var::*;
use super::multi_var::*;

fn rising_power_poly<R>(ring: &PolyRing<R>, n: usize) -> <PolyRing<R> as Ring>::El
    where R: DivisibilityInfoRing
{
    let mut result = ring.one();
    let x = ring.unknown();
    for i in 0..n {
        result = ring.mul(result, ring.add_ref(ring.from_z(i as i64), &x));
    }
    return result;
}

fn sumation_poly<R>(ring: &PolyRing<R>, f: <PolyRing<R> as Ring>::El) -> <PolyRing<R> as Ring>::El 
    where R: DivisibilityInfoRing
{
    let mut current = f;
    let mut result = ring.zero();
    let incl = embedding(ring.base_ring(), ring);
    while let Some(n) = ring.deg(&current) {
        let coeff = current.at(n).clone();
        current = ring.sub(current, ring.mul(rising_power_poly(ring, n), incl(coeff.clone())));
        let xn_sum_scaled = rising_power_poly(ring, n + 1);
        let xn_sum = ring.quotient(&xn_sum_scaled, &ring.from_z(n as i64 + 1)).unwrap();
        result = ring.add(result, ring.mul(xn_sum, incl(coeff)));
    }
    return result;
}

pub fn sumation_operator<R>(ring: &MultivariatePolyRing<R>, f: <MultivariatePolyRing<R> as Ring>::El, var: Var) -> <MultivariatePolyRing<R> as Ring>::El
    where R: DivisibilityInfoRing
{
    let (elevated_ring, iso, iso_inv) = ring.elevate_var(var);
    return iso_inv(sumation_poly(&elevated_ring, iso(f)));
}

#[cfg(test)]
use super::super::super::wrapper::*;
#[cfg(test)]
use super::super::super::primitive::*;
#[cfg(test)]
use super::super::rat::*;

#[test]
fn test_sumation_poly() {
    let ring = PolyRing::adjoint(r64::RING, "X");
    let x = ring.bind(ring.unknown());
    let sum_x_2 = x.clone() * (x.clone() + 1) / 2;
    assert_eq!(sum_x_2, ring.bind(sumation_poly(&ring, ring.unknown())));

    let sum_x_3 = x.clone() * (x.clone() + 1) * (x.clone() * 2 + 1) / 6;
    assert_eq!(sum_x_3, ring.bind(sumation_poly(&ring, (x.clone() * x.clone()).val().clone())));
}
