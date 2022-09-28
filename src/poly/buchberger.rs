use super::super::prelude::*;
use super::monomial_order::*;
use super::*;

use std::borrow::Cow;

type Monomial = Vec<usize>;

fn monomial_divide(lhs: &Monomial, rhs: &Monomial) -> Option<Monomial> {
    if rhs.len() > lhs.len() {
        return None;
    } else if lhs.iter().zip(rhs.iter()).any(|(l, r)| l < r) {
        return None;
    }
    Some((0..lhs.len()).map(|i| lhs[i] - rhs.get(i).unwrap_or(&0)).collect())
}

pub fn multi_poly_div<R, M>(ring: &R, mut lhs: El<R>, rhs: &El<R>, order: M) -> El<R>
    where R: MultiPolyRing, M: MonomialOrder + Copy
{
    assert!(ring.base_ring().is_field().can_use());
    let (rhs_lm, rhs_lc) = ring.lt(rhs, order).unwrap();
    while let Some((lhs_lm, lhs_lc)) = ring.lt(&lhs, order) {
        if let Some(m) = monomial_divide(&lhs_lm, &rhs_lm) {
            let factor = ring.base_ring().div(lhs_lc.clone(), rhs_lc);
            lhs = ring.sub_shifted(lhs, rhs, &m, &factor);   
        } else {
            return lhs;
        }
    }
    return lhs;
}

#[cfg(test)]
use super::super::rational::r64;

#[test]
fn test_multi_poly_divide() {
    let ring = MultivariatePolyRing::new(r64::RING, vec!["X", "Y"]);
    let ring = WrappingRing::new(&ring);
    let x= ring.as_poly("X");
    let y = ring.as_poly("Y");
    let f = &y * x.pow(5); // x^5 y
    let g = &y * x.pow(2) + y.pow(2) * &x * 2 + 1; // x^2 y + 2x y^2 + 1
    let h = -x.pow(3) + x.pow(2) * &y * 2 - &x * y.pow(2) * 4 - x.pow(2) * y.pow(4) * 8;
    // x^3 (x^2 y + 2x y^2 + 1) - 2x^4 y^2 - x^3
    // -2x^2 y (x^2 y + 2x y^2 + 1) + 4x^3 y^3 + 2x^2 y
    // 4x y^2 (x^2 y + 2x y^2 + 1) - 8x^2 y^4 - 4x y^2
    assert_eq!(h, ring.from(multi_poly_div(ring.wrapped_ring(), f.into_val(), g.val(), Lex {})))

}