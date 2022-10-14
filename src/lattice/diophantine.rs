#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::ring::*;
use super::super::wrapper::*;
use super::super::integer::*;

use std::cmp::min;

pub fn diophantine_solve<R, M, N, K, V>(partial_snf: Matrix<M, RingElWrapper<R>>, left_inv: Matrix<N, RingElWrapper<R>>, right_inv: Matrix<K, RingElWrapper<R>>, b: Vector<V, RingElWrapper<R>>) -> Option<Vector<MatrixCol<RingElWrapper<R>, MatrixOwned<RingElWrapper<R>>>, RingElWrapper<R>>>
    where R: IntegerRing, M: MatrixView<RingElWrapper<R>>, N: MatrixView<RingElWrapper<R>>, K: MatrixView<RingElWrapper<R>>, V: VectorView<RingElWrapper<R>>
{
    let ring = partial_snf.at(0, 0).ring().clone();
    let b = left_inv.mul(Matrix::col_vec(b), &ring).compute().into_col_vec();
    let mut x = Vector::zero_ring(partial_snf.col_count(), &ring).into_owned();
    for i in 0..min(partial_snf.row_count(), partial_snf.col_count()) {
        *x.at_mut(i) = ring.quotient(b.at(i), partial_snf.at(i, i))?;
    }
    for i in partial_snf.col_count()..partial_snf.row_count() {
        if !b.at(i).is_zero() {
            return None;
        }
    }
    let x = right_inv.mul(Matrix::col_vec(x), &ring).compute().into_col_vec();
    return Some(x);
}

#[cfg(test)]
use super::super::primitive::*;
#[cfg(test)]
use super::smith::*;

#[test]
fn test_diophantine() {
    let ring = i64::WRAPPED_RING;
    let mut A = Matrix::map(Matrix::from_array([[15, 10], [6, 7]]), ring.wrapping_embedding());
    let b = Vector::map(Vector::from_array([195, 87]), ring.wrapping_embedding());
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert_eq!(x.unwrap(), Vector::from_array([11, 3]));
}

#[test]
fn test_diophantine_no_solution() {
    let ring = i64::WRAPPED_RING;
    let i = ring.wrapping_embedding();
    let mut A = Matrix::from_array([[i(2), i(-2)]]);
    let b = Vector::from_array([i(1)]);
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert!(x.is_none());
}

#[test]
fn test_diophantine_no_solution_three_dim() {
    let ring = i64::WRAPPED_RING;
    #[rustfmt::skip]
    let mut A = Matrix::map(Matrix::from_array([[1, 2, 0], 
                                                [1, 0, 2]]), ring.wrapping_embedding());

    let b = Vector::map(Vector::from_array([2, 1]), ring.wrapping_embedding());
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert!(x.is_none());
}

#[test]
fn test_diophantine_three_dim() {
    let ring = i64::WRAPPED_RING;
    #[rustfmt::skip]
    let mut A = Matrix::map(Matrix::from_array([[1, 2, 0], 
                                                [1, 0, 2]]), ring.wrapping_embedding());

    let b = Vector::map(Vector::from_array([2, 4]), ring.wrapping_embedding());
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert_eq!(x.unwrap(), Vector::from_array([4, -1, 0]));
}

#[test]
fn test_diophantine_unnecessary_conditions() {
    let ring = i64::WRAPPED_RING;
    #[rustfmt::skip]
    let mut A = Matrix::map(Matrix::from_array([[1, 2, 0], 
                                                [1, 2, 0], 
                                                [1, 2, 0], 
                                                [1, 0, 2]]), ring.wrapping_embedding());

    let b = Vector::map(Vector::from_array([2, 2, 2, 4]), ring.wrapping_embedding());
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert_eq!(x.unwrap(), Vector::from_array([2, 0, 1]));
}

#[test]
fn test_diophantine_no_rational_solutions() {
    let ring = i64::WRAPPED_RING;
    #[rustfmt::skip]
    let mut A = Matrix::map(Matrix::from_array([[1, 2, 0], 
                                                [1, 2, 0], 
                                                [1, 0, 2]]), ring.wrapping_embedding());

    let b = Vector::map(Vector::from_array([2, 3, 4]), ring.wrapping_embedding());
    let (iL, iR) = partial_snf(&mut A);
    let x = diophantine_solve(A.as_ref(), iL.as_ref(), iR.as_ref(), b.as_ref());
    assert!(x.is_none());
}
