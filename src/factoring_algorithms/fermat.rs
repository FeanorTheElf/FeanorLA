use super::super::bigint::*;
use super::super::integer::*;

///
/// Uses the Fermat method to find nontrivial factors of a given number.
/// Use only for numbers that have two odd factors of approximately
/// equal size. The running time is in O(d^2/sqrt(n)), where d is
/// the difference of both factors.
/// 
/// # Panics
/// 
/// Will panic (after a potentially very long computation) if n does not
/// have two odd factors.
/// 
pub fn fermat(n: &BigInt) -> BigInt
{
    let mut m = BigInt::RING.root_floor(&(n.clone() - 1), 2) + 1;

    // We try to find m such that m^2 - n = d^2 for some integer d, then
    // (m - d)(m + d) = n.
    // the smallest possible m is sqrt(n), for the largest have m^2 = n + d^2;
    // Using sqrt(1 + x^2) - 1 <= x^2 / 2 we get that the number of possible m is
    // <= d^2 / 2 / sqrt(n)
    loop {
        assert!(m < *n);
        let m_squared = m.pow(2);
        let d = BigInt::RING.root_floor(&(m_squared.clone() - n.clone()), 2);
        if m_squared - d.pow(2) == *n {
            return m - d;
        }
        m += 1;
    }
}

#[test]
fn test_fermat_factorization() {
    assert_eq!(BigInt::from(9851), fermat(&BigInt::from(99997501)));
}