use super::super::prelude::*;
use super::monomial_order::*;
use super::*;
use super::super::rational::*;

use std::collections::BinaryHeap;
use std::borrow::Cow;

type Monomial = Vec<usize>;

fn monomial_divide(lhs: &Monomial, rhs: &Monomial) -> Option<Monomial> {
    if rhs.len() > lhs.len() {
        return None;
    } else if lhs.iter().zip(rhs.iter()).any(|(l, r)| l < r) {
        return None;
    }
    Some(normalized_monomial((0..lhs.len()).map(|i| lhs[i] - rhs.get(i).unwrap_or(&0))))
}

fn monomial_correlation(lhs: &Monomial, rhs: &Monomial) -> i64 {
    let inner_prod = |l: &Monomial, r: &Monomial| l.iter().zip(r.iter()).map(|(a, b)| a * b).sum::<usize>() as i64;
    (inner_prod(lhs, rhs) * inner_prod(lhs, rhs)) / (inner_prod(lhs, lhs) * inner_prod(rhs, rhs))
}

///
/// Computes minimal a, b such that `a * lhs == b * rhs`
/// 
fn monomial_lcm_combine(lhs: &Monomial, rhs: &Monomial) -> (Monomial, Monomial) {
    let lcm = (0..usize::max(lhs.len(), rhs.len()))
        .map(|i| usize::max(
            *lhs.get(i).unwrap_or(&0), 
            *rhs.get(i).unwrap_or(&0)
        ));
    (
        normalized_monomial(lcm.clone().enumerate().map(|(i, x)| x - lhs.get(i).unwrap_or(&0))),
        normalized_monomial(lcm.enumerate().map(|(i, x)| x - rhs.get(i).unwrap_or(&0))),
    )
}

pub fn multi_poly_div<R, M>(ring: &R, mut lhs: El<R>, rhs: &El<R>, order: M) -> El<R>
    where R: MultiPolyRing, M: MonomialOrder + Copy
{
    assert!(ring.base_ring().is_field().can_use());
    let (rhs_lm, rhs_lc) = ring.lt(rhs, order).unwrap();
    while let Some((lhs_lm, lhs_lc)) = ring.lt(&lhs, order) {
        if let Some(m) = monomial_divide(&lhs_lm, &rhs_lm) {
            let factor = ring.base_ring().div(lhs_lc.clone(), rhs_lc);
            lhs = ring.sub_scaled(lhs, rhs, &m, &factor);   
        } else {
            return lhs;
        }
    }
    return lhs;
}

pub fn buchberger<R, M>(ring: &R, mut basis: Vec<El<R>>, order: M) -> Vec<El<R>>
    where R: MultiPolyRing, M: MonomialOrder + Copy
{
    let mut ij_pairs = BinaryHeap::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            ij_pairs.push((monomial_correlation(&ring.lm(&basis[i], order).unwrap(), &ring.lm(&basis[j], order).unwrap()), i, j));
        }
    }
    while let Some((_, i, j)) = ij_pairs.pop() {
        let (lhs_lm, lhs_lc) = ring.lt(&basis[i], order).unwrap();
        let (rhs_lm, rhs_lc) = ring.lt(&basis[j], order).unwrap();
        let (a, b) = monomial_lcm_combine(&lhs_lm, &rhs_lm);
        let S = ring.sub(
            ring.scale(basis[i].clone(), &a, rhs_lc),
            ring.scale(basis[j].clone(), &b, lhs_lc)
        );
        let mut result = S;
        for b in &basis {
            result = multi_poly_div(ring, result, b, order);
        }
        if !ring.is_zero(&result) {
            let m = ring.lm(&result, order).unwrap();
            println!("{}", ring.display(&result));
            ij_pairs.extend((0..basis.len()).map(|i| (monomial_correlation(&ring.lm(&basis[i], order).unwrap(), &m), i, basis.len())));
            basis.push(result);
        }
    }
    return basis;
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

#[test]
fn test_groebner() {
    let ring = MultivariatePolyRing::new(r64::RING, vec!["X", "Y"]);
    let ring = WrappingRing::new(&ring);
    let x= ring.as_poly("X");
    let y = ring.as_poly("Y");
;
    let f = &y * x.pow(5); // x^5 y
    let g = &y * x.pow(2) + y.pow(2) * &x * 2 + 1; // x^2 y + 2x y^2 + 1
    let gb = buchberger(ring.wrapped_ring(), vec![f.into_val(), g.into_val()], Lex {});
    for h in &gb {
        println!("{}", ring.wrapped_ring().display(h))
    }
}