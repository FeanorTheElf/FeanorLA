use super::super::alg::*;
use super::bigint::*;
use super::rat::*;
use super::primality::*;

///
/// Returns an integer n such that p^n divides x and p^(n + 1) does not.
/// 
/// # Termination & Noetherian rings
/// 
/// Returns None if the valuation is infinite. Note that if the input
/// ring is not noetherian, this function might loop forever.
/// 
/// # Complexity
/// 
/// The algorithm runs in O(log(v)^2) where v is the p-adic valuation of x.
/// 
pub fn p_adic_valuation<R>(ring: &R, x: R::El, p: &R::El) -> Option<u32>
    where R: DivisibilityInformationRing
{
    assert!(ring.is_integral());
    // using integrality and noetherianity, it is easy to show
    // that in all other cases, there exists some n with p^n does
    // not divide x
    if ring.is_unit(p) || ring.is_zero(&x) {
        return None;
    }
    let mut current = x;
    let mut valuation = 0;
    while let Some(new_current) = ring.quotient(&current, p) {
        valuation += 1;
        current = new_current;
        let mut power_p = ring.mul_ref(p, p);
        let mut valuation_delta = 2;
        while let Some(new_current) = ring.quotient(&current, &power_p) {
            current = new_current;
            valuation += valuation_delta;
            let new_power_p = ring.mul_ref(&power_p, &power_p);
            power_p = new_power_p;
            valuation_delta = valuation_delta * 2;
        }
    }
    return Some(valuation);
}

pub fn p_adic_abs(x: r64, p: i64) -> f64 {
    if let Some(num_val) = p_adic_valuation(&i64::RING, x.num(), &p) {
        if let Some(den_val) = p_adic_valuation(&i64::RING, x.den(), &p) {
            (p as f64).pow(den_val) / (p as f64).pow(num_val)
        } else {
            f64::INFINITY
        }
    } else {
        0.
    }
}

#[test]
fn test_p_adic_valuation() {
    assert_eq!(Some(3), p_adic_valuation(&i64::RING, 7 * 7 * 7 * 5 * 2 * 2 * 3, &7));
    assert_eq!(None, p_adic_valuation(&i64::RING, 0, &7));
}

#[test]
fn test_p_adic_abs() {
    assert_eq!(1./9., p_adic_abs(r64::new(36, 1), 3));
    assert_eq!(0., p_adic_abs(r64::ZERO, 3));
}