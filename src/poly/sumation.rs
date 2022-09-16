use super::super::ring::*;
use super::*;

fn rising_power_poly<P>(ring: &P, n: usize) -> El<P>
    where P: PolyRing
{
    let mut result = ring.one();
    let x = ring.unknown();
    for i in 0..n {
        result = ring.mul(result, ring.add_ref(ring.from_z(i as i64), &x));
    }
    return result;
}

pub fn sumation_poly<P>(ring: &P, f: El<P>) -> El<P> 
    where P: PolyRing + DivisibilityInfoRing
{
    assert!(ring.is_divisibility_computable().can_use());
    let mut current = f;
    let mut result = ring.zero();
    let incl = ring.embedding();
    while let Some(n) = ring.deg(&current) {
        let coeff = ring.coefficient_at(&current, n).clone();
        current = ring.sub(current, ring.mul(rising_power_poly(ring, n), incl(coeff.clone())));
        let xn_sum_scaled = rising_power_poly(ring, n + 1);
        let xn_sum = ring.quotient(&xn_sum_scaled, &ring.from_z(n as i64 + 1)).unwrap();
        result = ring.add(result, ring.mul(xn_sum, incl(coeff)));
    }
    return result;
}

#[cfg(test)]
use super::super::rational::*;

#[test]
fn test_sumation_poly() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(r64::RING, "X"));
    let x = ring.unknown();
    let sum_x_2 = x.clone() * (x.clone() + 1) / 2;
    assert_eq!(sum_x_2, ring.from(sumation_poly(ring.wrapped_ring(), ring.unknown().into_val())));

    let sum_x_3 = x.clone() * (x.clone() + 1) * (x.clone() * 2 + 1) / 6;
    assert_eq!(sum_x_3, ring.from(sumation_poly(ring.wrapped_ring(), (&x * &x).into_val())));
}
