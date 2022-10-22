use super::super::prelude::*;
use super::monomial_order::*;
use super::*;

type Monomial = Vec<usize>;

fn monomial_divide(lhs: &Monomial, rhs: &Monomial) -> Option<Monomial> {
    if rhs.len() > lhs.len() {
        return None;
    } else if lhs.iter().zip(rhs.iter()).any(|(l, r)| l < r) {
        return None;
    }
    Some(normalized_monomial((0..lhs.len()).map(|i| lhs[i] - rhs.get(i).unwrap_or(&0))))
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
            let factor = ring.base_ring().neg(ring.base_ring().div(lhs_lc.clone(), rhs_lc));
            lhs = ring.add_scaled(lhs, rhs, &m, &factor);   
        } else {
            return lhs;
        }
    }
    return lhs;
}

fn reduce_basis<R, M>(ring: &R, basis: &mut Vec<El<R>>, order: M)
    where R: MultiPolyRing, R::BaseRing: DivisibilityInfoRing, M: MonomialOrder + Copy
{
    for i in 0..basis.len() {
        let mut current = std::mem::replace(&mut basis[i], ring.zero());
        for j in 0..basis.len() {
            if !ring.is_zero(&basis[j]) {
                current = multi_poly_div(ring, current, &basis[j], order);
            }
        }
        basis[i] = ring.reduced_poly(&current);
    }
}

pub fn buchberger<R, M>(ring: &R, mut basis: Vec<El<R>>, order: M) -> Vec<El<R>>
    where R: MultiPolyRing, R::BaseRing: DivisibilityInfoRing, M: MonomialOrder + Copy
{
    let mut ij_pairs = Vec::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            ij_pairs.push((i, j));
        }
    }
    while let Some((i, j)) = ij_pairs.pop() {
        if ring.is_zero(&basis[i]) || ring.is_zero(&basis[j]) {
            continue;
        }
        let (lhs_lm, lhs_lc) = ring.lt(&basis[i], order).unwrap();
        let (rhs_lm, rhs_lc) = ring.lt(&basis[j], order).unwrap();
        let (a, b) = monomial_lcm_combine(&lhs_lm, &rhs_lm);
        let s_poly = ring.sub(
            ring.scale(basis[i].clone(), &a, rhs_lc),
            ring.scale(basis[j].clone(), &b, lhs_lc)
        );
        let mut result = s_poly;
        for b in &basis {
            if !ring.is_zero(b) {
                result = multi_poly_div(ring, result, b, order);
            }
        }
        if !ring.is_zero(&result) {
            ij_pairs.extend((0..basis.len()).map(|i| (i, basis.len())));
            result = ring.reduced_poly(&result);
            basis.push(result);
            reduce_basis(ring, &mut basis, order);
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
    assert_eq!(h, ring.wrap(multi_poly_div(ring.wrapped_ring(), f.into_val(), g.val(), Lex {})))
}

#[test]
fn test_groebner() {
    let ring = MultivariatePolyRing::new(r64::RING, vec!["X", "Y"]);
    let ring = WrappingRing::new(&ring);
    let x= ring.as_poly("X");
    let y = ring.as_poly("Y");

    let f = &y * x.pow(5);
    let g = &y * x.pow(2) + y.pow(2) * &x + 1;
    let gb = buchberger(ring.wrapped_ring(), vec![f.into_val(), g.into_val()], Lex {});
    assert!(gb.iter().any(|x| ring.wrapped_ring().is_one(x)));
}