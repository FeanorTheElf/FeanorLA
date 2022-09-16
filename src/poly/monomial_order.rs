use std::cmp::Ordering;

///
/// We assume the same contract as for the keys in the elements of
/// `MultivariatePolyRing`. In particular, the monomial should have no
/// trailing zeros.
/// 
type Monomial = Vec<usize>;

pub trait MonomialOrder: Clone {

    fn cmp(&self, lhs: &Monomial, rhs: &Monomial) -> Ordering;
}

#[derive(Clone, Copy)]
pub struct Lex;

impl MonomialOrder for Lex {
    
    fn cmp(&self, lhs: &Monomial, rhs: &Monomial) -> Ordering {
        for i in 0..usize::min(lhs.len(), rhs.len()) {
            if lhs[i] < rhs[i] {
                return Ordering::Less;
            } else if lhs[i] > rhs[i] {
                return Ordering::Greater;
            }
        }
        // this works, since we assume that lhs and rhs have no trailing zeros
        return usize::cmp(&lhs.len(), &rhs.len());
    }
}

#[derive(Clone, Copy)]
pub struct GrLex;

impl MonomialOrder for GrLex {
    
    fn cmp(&self, lhs: &Monomial, rhs: &Monomial) -> Ordering {
        match usize::cmp(&lhs.iter().sum(), &rhs.iter().sum()) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => Lex.cmp(lhs, rhs)
        }
    }
}

#[derive(Clone, Copy)]
pub struct GrRevLex;

impl MonomialOrder for GrRevLex {
    
    fn cmp(&self, lhs: &Monomial, rhs: &Monomial) -> Ordering {
        match usize::cmp(&lhs.iter().sum(), &rhs.iter().sum()) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            // this works, since we assume that lhs and rhs have no trailing zeros
            Ordering::Equal => match usize::cmp(&lhs.len(), &rhs.len()) {
                Ordering::Less => Ordering::Greater,
                Ordering::Greater => Ordering::Less,
                Ordering::Equal => {
                    for (l, r) in lhs.iter().rev().zip(rhs.iter().rev()) {
                        if l < r {
                            return Ordering::Greater;
                        } else if l > r {
                            return Ordering::Less;
                        }
                    }
                    return Ordering::Equal;
                }
            }
        }
    }
}

#[test]
fn test_lex() {
    let mut monomials = vec![
        vec![1], vec![2], vec![1, 1], vec![1, 0, 1], vec![0, 2], vec![0, 1, 1], vec![0, 0, 2]
    ];
    monomials.sort_by(|a, b| Lex.cmp(a, b));
    assert_eq!(vec![
        vec![0, 0, 2], vec![0, 1, 1], vec![0, 2], vec![1], vec![1, 0, 1], vec![1, 1], vec![2]
    ], monomials);
}

#[test]
fn test_grlex() {
    let mut monomials = vec![
        vec![1], vec![2], vec![1, 1], vec![1, 0, 1], vec![0, 2], vec![0, 1, 1], vec![0, 0, 2]
    ];
    monomials.sort_by(|a, b| GrLex.cmp(a, b));
    assert_eq!(vec![
        vec![1], vec![0, 0, 2], vec![0, 1, 1], vec![0, 2], vec![1, 0, 1], vec![1, 1], vec![2]
    ], monomials);
}

#[test]
fn test_grrevlex() {
    let mut monomials = vec![
        vec![1], vec![2], vec![1, 1], vec![1, 0, 1], vec![0, 2], vec![0, 1, 1], vec![0, 0, 2]
    ];
    monomials.sort_by(|a, b| GrRevLex.cmp(a, b));
    assert_eq!(vec![
        vec![1], vec![0, 0, 2], vec![0, 1, 1], vec![1, 0, 1], vec![0, 2], vec![1, 1], vec![2]
    ], monomials);
}
