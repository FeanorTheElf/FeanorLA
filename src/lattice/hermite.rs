#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::eea::*;
use super::super::ring::*;
use super::super::wrapper::*;
use super::super::embedding::*;
use super::super::integer::*;

fn row_echelon_form<R, M, F>(A: &mut Matrix<M, RingElWrapper<R>>, mut transformed_rows: F)
    where R: EuclideanInfoRing, M: MatrixViewMut<RingElWrapper<R>>, F: FnMut(usize, usize, [RingElWrapper<R>; 4])
{
    let mut i = 0;
    let mut j = 0;
    while i < A.row_count() && j < A.col_count() {
        let mut last_row = 1;
        while let Some(index) = find_transform_row(A.submatrix(i.., j..), last_row) {
            last_row = index + 1;
            let a = A.at(i, j).clone();
            let b = A.at(index + i, j).clone();
            let (s, t, d) = eea(&a.ring(), a.clone(), b.clone());
            let transform = [s, t, -b / &d, a / &d];
            // taking only the submatrix here is not necessary, but increases performance
            A.submatrix_mut(i.., j..).transform_two_dims_left(0, index, &transform, &d.ring());
            transformed_rows(i, index + i, transform);
        }
        if !A.at(i, j).is_zero() {
            i += 1;
        }
        j += 1;
    }
}

pub fn row_hnf<R, M>(A: &mut Matrix<M, RingElWrapper<R>>) -> Matrix<MatrixOwned<RingElWrapper<R>>, RingElWrapper<R>>
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>
{
    let ring = A.at(0, 0).ring().clone();
    let incl = z_hom(&ring);
    let mut U = Matrix::identity_ring(A.row_count(), A.row_count(), &ring).to_owned();
    let mut transform_U = |i, j, [a, b, c, d]: [RingElWrapper<R>; 4]| {
        let inv_det = (&a * &d - &b * &c).inv();
        let col_transform = [d * &inv_det, -b * &inv_det, -c * &inv_det, a * &inv_det];
        U.transform_two_dims_right(i, j, &col_transform, &ring);
    };
    row_echelon_form(A, &mut transform_U);
    println!("{}", A);
    let mut i = 0;
    let mut j = 0;
    while i < A.row_count() && j < A.col_count() {
        if A.at(i, j).is_zero() {
            j += 1;
            continue;
        }
        let pivot = A.at(i, j).clone().abs();
        let mut pivot_sign = A.at(i, j).sgn();
        for k in 0..i {
            let factor = A.at(k, j).clone().floor_div(&pivot);
            let transform = [incl(1), -factor * &pivot_sign, incl(0), pivot_sign.clone()];
            A.submatrix_mut(..=i, j..).transform_two_dims_left(k, i, &transform, &ring);
            transform_U(k, i, transform);
            println!("{}, {}, {}", A, i, j);
            pivot_sign = ring.one();
        }
        i += 1;
        j += 1;
    }
    return U;
}

fn find_transform_row<R, M>(A: Matrix<M, RingElWrapper<R>>, start_col: usize) -> Option<usize>
    where R: EuclideanInfoRing, M: MatrixView<RingElWrapper<R>>
{
    let pivot = A.at(0, 0);
    for i in start_col..A.row_count() {
        if !A.at(i, 0).divides(pivot) {
            return Some(i);
        }
    }
    return None;
}

#[cfg(test)]
use super::super::primitive::*;

#[test]
fn test_row_triangular() {
    let ring = i64::RING.bind_ring_by_value();
    let mut A = Matrix::map(Matrix::from_array([[ 6, 5 ], [ 4, 10 ]]), ring.wrapping_embedding());
    let expected = Matrix::map(Matrix::from_array([[ 2, -5 ], [ 0, 20 ]]), ring.wrapping_embedding());
    row_echelon_form(&mut A, |_, _, _| {});
    assert_eq!(expected, A);

    let mut A = Matrix::map(Matrix::from_array([[ 6, 5, 9 ], [ 4, 10, 11 ], [ 7, 5, 2 ]]), ring.wrapping_embedding());
    let expected = Matrix::map(Matrix::from_array([[ 1, 20, 8 ], [ 0, 5, -12 ], [ 0, 0, -63 ]]), ring.wrapping_embedding());
    row_echelon_form(&mut A, |_, _, _| {});
    println!("{}", A);
    assert_eq!(expected, A);
}

#[test]
fn test_row_hnf() {
    let ring = i64::RING.bind_ring_by_value();
    let A = Matrix::map(Matrix::from_array([[ 6, 5, 9 ], [ 4, 10, 11 ], [ 7, 5, 2 ]]), ring.wrapping_embedding());
    let expected = Matrix::map(Matrix::from_array([[ 1, 0, 56 ], [ 0, 5, 51 ], [ 0, 0, 63 ]]), ring.wrapping_embedding());
    let mut actual = A.clone();
    let U = row_hnf(&mut actual);
    assert_eq!(expected, actual);
    assert_eq!(A, U.mul(actual, &ring));
}