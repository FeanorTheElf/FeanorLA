#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::eea::*;
use super::super::wrapper::*;
use super::super::embedding::*;
use super::super::integer::*;

pub fn partial_snf<R, M>(A: &mut Matrix<M, RingElWrapper<R>>) -> (Matrix<MatrixOwned<RingElWrapper<R>>, RingElWrapper<R>>, Matrix<MatrixOwned<RingElWrapper<R>>, RingElWrapper<R>>)
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>
{
    let ring = A.at(0, 0).ring();
    let mut iL = Matrix::identity_ring(A.row_count(), A.row_count(), &ring);
    let mut iR = Matrix::identity_ring(A.col_count(), A.col_count(), &ring);
    partial_smith(A.as_mut(), iL.as_mut(), iR.as_mut(), 0);
    return (iL, iR);
}

///
/// Transforms the matrix A into diagonal form and
/// changes L, R so that L' * A' * R' = L * A * R
/// and |det L'| = |det L|, |det R'| = |det R| holds
/// Instead of L and R, this function works on their
/// inverses iL and iR.
///
fn partial_smith<R, M, N, K>(
    mut A: Matrix<M, RingElWrapper<R>>, 
    mut iL: Matrix<N, RingElWrapper<R>>, 
    mut iR: Matrix<K, RingElWrapper<R>>, 
    pivot: usize
) 
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>, K: MatrixViewMut<RingElWrapper<R>>
{
    let ring = A.at(0, 0).ring().clone();
    if pivot == A.row_count() || pivot == A.col_count() {
        return;
    }
    let is_zero_matrix = swap_pivot_entry_if_zero(
        A.as_mut(), iL.as_mut(), iR.as_mut(), pivot
    );
    // pivot must be != 0
    if is_zero_matrix {
        return;
    }
    // pivot must divide all entries on pivot row and pivot column
    let mut changed = true;
    while changed {
        changed = transform_pivot_gcd_col(A.as_mut(), iL.as_mut(), pivot, &ring) || 
            transform_pivot_gcd_row(A.as_mut(), iR.as_mut(), pivot, &ring);
    }
    // eliminate the entries on pivot row and pivot column
    eliminate_col(A.as_mut(), iL.as_mut(), pivot, &ring);
    eliminate_row(A.as_mut(), iR.as_mut(), pivot, &ring);
    partial_smith(A, iL, iR, pivot + 1);
}

fn eliminate_col<R, M, N>(mut A: Matrix<M, RingElWrapper<R>>, mut iL: Matrix<N, RingElWrapper<R>>, pivot: usize, ring: &WrappingRing<R>)
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>
{
    let i = z_hom(ring);
    for row in (pivot + 1)..A.row_count() {
        let transform = [i(1), i(0), -A.at(row, pivot) / A.at(pivot, pivot), i(1)];
        A.transform_two_dims_left(pivot, row, &transform, ring);
        iL.transform_two_dims_left(pivot, row, &transform, ring);
    }
}

fn eliminate_row<R, M, N>(mut A: Matrix<M, RingElWrapper<R>>, mut iR: Matrix<N, RingElWrapper<R>>, pivot: usize, ring: &WrappingRing<R>)
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>
{
    let i = z_hom(ring);
    for col in (pivot + 1)..A.col_count() {
        let transform = [i(1), -A.at(pivot, col) / A.at(pivot, pivot), i(0), i(1)];
        A.transform_two_dims_right(pivot, col, &transform, ring);
        iR.transform_two_dims_right(pivot, col, &transform, ring);
    }
}

fn transform_pivot_gcd_col<R, M, N>(
    mut A: Matrix<M, RingElWrapper<R>>, 
    mut iL: Matrix<N, RingElWrapper<R>>, 
    pivot: usize,
    ring: &WrappingRing<R>
) -> bool
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>
{
    let pivot_row = pivot;
    let pivot_col = pivot;
    let mut current =
        find_smallest_gcd_entry_in_pivot_col(
            A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count()), ring
        );
    if current == 0 {
        return false;
    }
    while current != 0 {
        let (a, b) = (A.at(pivot_row, pivot_col).clone(), A.at(pivot_row + current, pivot_col).clone());
        let (s, t, d) = signed_eea(a.clone(), b.clone(), ring);
        let transform = [s, t, -b / &d, a / d];
        A.transform_two_dims_left(pivot_row, pivot_row + current, &transform, ring);
        iL.transform_two_dims_left(pivot_row, pivot_row + current, &transform, ring);
        current =
            find_smallest_gcd_entry_in_pivot_col(
                A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count()), ring
            );
    }
    return true;
}

fn transform_pivot_gcd_row<R, M, N>(
    mut A: Matrix<M, RingElWrapper<R>>, 
    mut iR: Matrix<N, RingElWrapper<R>>, 
    pivot: usize,
    ring: &WrappingRing<R>
) -> bool 
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>
{
    let pivot_row = pivot;
    let pivot_col = pivot;
    let mut current =
        find_smallest_gcd_entry_in_pivot_row(
            A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count()), ring
        );
    if current == 0 {
        return false;
    }
    while current != 0 {
        let (a, b) = (A.at(pivot_row, pivot_col).clone(), A.at(pivot_row, pivot_col + current).clone());
        let (s, t, d) = signed_eea(a.clone(), b.clone(), ring);
        let transform = [s, -b / &d, t, a / d];
        A.transform_two_dims_right(pivot_col, pivot_col + current, &transform, ring);
        iR.transform_two_dims_left(pivot_col, pivot_col + current, &transform, ring);
        current =
            find_smallest_gcd_entry_in_pivot_row(
                A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count()), ring
            );
    }
    return true;
}

fn find_min_by_key<T, I, F, C>(mut it: I, mut f: F) -> Option<T>
where
    I: Iterator<Item = T>,
    F: FnMut(&T) -> C,
    C: Ord
{
    let mut result: T = it.next()?;
    let mut current_val: C = f(&result);
    for item in it {
        let value = f(&item);
        if value < current_val {
            result = item;
            current_val = value;
        }
    }
    return Some(result);
}

fn find_smallest_gcd_entry_in_pivot_row<R, M>(A: Matrix<M, RingElWrapper<R>>, ring: &WrappingRing<R>) -> usize 
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>
{
    find_min_by_key(0..A.col_count(), |col: &usize| 
        signed_gcd(A.at(0, 0).clone(), A.at(0, *col).clone(), ring)
    ).unwrap()
}

fn find_smallest_gcd_entry_in_pivot_col<R, M>(A: Matrix<M, RingElWrapper<R>>, ring: &WrappingRing<R>) -> usize
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>
{
    find_min_by_key(0..A.row_count(), |row: &usize| 
        signed_gcd(A.at(0, 0).clone(), A.at(*row, 0).clone(), ring)
    ).unwrap()
}

fn swap_pivot_entry_if_zero<R, M, N, K>(
    mut A: Matrix<M, RingElWrapper<R>>,
    mut iL: Matrix<N, RingElWrapper<R>>,
    mut iR: Matrix<K, RingElWrapper<R>>,
    pivot: usize,
) -> bool 
    where R: IntegerRing, M: MatrixViewMut<RingElWrapper<R>>, N: MatrixViewMut<RingElWrapper<R>>, K: MatrixViewMut<RingElWrapper<R>>
{
    let pivot_row = pivot;
    let pivot_col = pivot;
    if let Some((row, col)) =
        find_not_zero(A.submatrix_mut(pivot_row.., pivot_col..))
    {
        A.swap_rows(pivot_row, row + pivot_row);
        iL.swap_rows(pivot_row, row + pivot_row);
        A.swap_cols(pivot_col, col + pivot_col);
        iR.swap_cols(pivot_col, col + pivot_col);
        return false;
    } else {
        return true;
    }
}

fn find_not_zero<R, M>(mat: Matrix<M, RingElWrapper<R>>) -> Option<(usize, usize)> 
    where R: IntegerRing, M: MatrixView<RingElWrapper<R>>
{
    for row in 0..mat.row_count() {
        for col in 0..mat.col_count() {
            if !mat.at(row, col).is_zero() {
                return Some((row, col));
            }
        }
    }
    return None;
}

#[cfg(test)]
use super::super::primitive::*;

#[test]
fn test_partial_smith_4x2() {
    let ring = i64::WRAPPED_RING;
    #[rustfmt::skip]
    let mut A = Matrix::map(Matrix::from_array([[ 1,  1, -1, -1],
                                                [ 1,  0, -1,  0]]), ring.wrapping_embedding());
    let A_copy = A.clone();
    let mut iL = Matrix::identity_ring(2, 2, &ring);
    let mut iR = Matrix::identity_ring(4, 4, &ring);
    partial_smith(A.as_mut(), iL.as_mut(), iR.as_mut(), 0);

    assert_eq!(A, iL.mul(A_copy, &ring).mul(iR, &ring));
}