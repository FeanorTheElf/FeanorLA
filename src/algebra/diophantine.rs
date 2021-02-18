#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::alg::*;
use super::eea::{signed_eea, signed_gcd};

type Item = i32;

pub fn diophantine_solve<M, V>(
    A: Matrix<M, Item>, 
    b: Vector<V, Item>
) -> Option<Vector<VectorOwned<Item>, Item>> 
    where M: MatrixView<Item>, V: VectorView<Item>
{
    let mut smith_A = A.to_owned();
    let mut iL = Matrix::identity(A.row_count(), A.row_count());
    let mut iR = Matrix::identity(A.col_count(), A.col_count());
    partial_smith(
        smith_A.as_mut(),
        iL.as_mut(),
        iR.as_mut(),
        0,
    );
    // x is solution of (L * smith_A) x = b, get result through r := R^-1 * x
    let mut x = Vector::zero(A.col_count());
    let c = iL * Matrix::col_vec(b);
    for i in 0..usize::min(x.len(), A.row_count()) {
        let entry = *smith_A.at(i, i);
        if entry == 0 && *c.at(i, 0) != 0 {
            return None;
        } else if entry != 0 && *c.at(i, 0) % entry != 0 {
            return None;
        } else if entry != 0 {
            *x.at_mut(i) = *c.at(i, 0) / entry;
        }
    }
    return Some((iR * Matrix::col_vec(x)).copy_vec());
}

///
/// Transforms the matrix A into diagonal form and
/// changes L, R so that L' * A' * R' = L * A * R
/// and |det L'| = |det L|, |det R'| = |det R| holds
/// Instead of L and R, this function works on their
/// inverses iL and iR.
///
pub fn partial_smith<M, N, K>(
    mut A: Matrix<M, Item>, 
    mut iL: Matrix<N, Item>, 
    mut iR: Matrix<K, Item>, 
    pivot: usize
) 
    where M: MatrixViewMut<Item>, N: MatrixViewMut<Item>, K: MatrixViewMut<Item>
{
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
        changed = transform_pivot_gcd_col(A.as_mut(), iL.as_mut(), pivot) || 
            transform_pivot_gcd_row(A.as_mut(), iR.as_mut(), pivot);
    }
    // eliminate the entries on pivot row and pivot column
    eliminate_col(A.as_mut(), iL.as_mut(), pivot);
    eliminate_row(A.as_mut(), iR.as_mut(), pivot);
    partial_smith(A, iL, iR, pivot + 1);
}

fn eliminate_col<M, N>(mut A: Matrix<M, Item>, mut iL: Matrix<N, Item>, pivot: usize)
    where M: MatrixViewMut<Item>, N: MatrixViewMut<Item>
{
    for row in (pivot + 1)..A.row_count() {
        let transform = [1, 0, -A.at(row, pivot) / A.at(pivot, pivot), 1];
        A.transform_two_dims_left(pivot, row, &transform);
        iL.transform_two_dims_left(pivot, row, &transform);
    }
}

fn eliminate_row<M, N>(mut A: Matrix<M, Item>, mut iR: Matrix<N, Item>, pivot: usize)
    where M: MatrixViewMut<Item>, N: MatrixViewMut<Item>
{
    for col in (pivot + 1)..A.col_count() {
        let transform = [1, -A.at(pivot, col) / A.at(pivot, pivot), 0, 1];
        A.transform_two_dims_right(pivot, col, &transform);
        iR.transform_two_dims_left(pivot, col, &transform);
    }
}

fn transform_pivot_gcd_col<M, N, T>(
    mut A: Matrix<M, T>, 
    mut iL: Matrix<N, T>, 
    pivot: usize
) -> bool
    where M: MatrixViewMut<T>, N: MatrixViewMut<T>, T: Integer
{
    let pivot_row = pivot;
    let pivot_col = pivot;
    let mut current =
        find_smallest_gcd_entry_in_pivot_col(
            A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count())
        );
    if current == 0 {
        return false;
    }
    while current != 0 {
        let (a, b) = (A.at(pivot_row, pivot_col).clone(), A.at(pivot_row + current, pivot_col).clone());
        let (s, t, _) = signed_eea(a.clone(), b.clone());
        let gcd = s.clone() * a.clone() + t.clone() * b.clone();
        let transform = [s, t, -b / gcd.clone(), a / gcd];
        A.transform_two_dims_left(pivot_row, pivot_row + current, &transform);
        iL.transform_two_dims_left(pivot_row, pivot_row + current, &transform);
        current =
            find_smallest_gcd_entry_in_pivot_col(
                A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count())
            );
    }
    return true;
}

fn transform_pivot_gcd_row<M, N, T>(
    mut A: Matrix<M, T>, 
    mut iR: Matrix<N, T>, 
    pivot: usize
) -> bool 
    where M: MatrixViewMut<T>, N: MatrixViewMut<T>, T: Integer
{
    let pivot_row = pivot;
    let pivot_col = pivot;
    let mut current =
        find_smallest_gcd_entry_in_pivot_row(
            A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count())
        );
    if current == 0 {
        return false;
    }
    while current != 0 {
        let (a, b) = (A.at(pivot_row, pivot_col).clone(), A.at(pivot_row, pivot_col + current).clone());
        let (s, t, _) = signed_eea(a.clone(), b.clone());
        let gcd = s.clone() * a.clone() + t.clone() * b.clone();
        let transform = [s, -b / gcd.clone(), t, a / gcd];
        A.transform_two_dims_right(pivot_col, pivot_col + current, &transform);
        iR.transform_two_dims_left(pivot_col, pivot_col + current, &transform);
        current =
            find_smallest_gcd_entry_in_pivot_row(
                A.submatrix_mut(pivot..A.row_count(), pivot..A.col_count())
            );
    }
    return true;
}

fn find_min<T, I, F, Int>(mut it: I, mut f: F) -> Option<T>
where
    I: Iterator<Item = T>,
    F: FnMut(&T) -> Int,
    Int: Integer
{
    let mut result: T = it.next()?;
    let mut current_val: Int = f(&result);
    for item in it {
        let value = f(&item);
        if value < current_val {
            result = item;
            current_val = value;
        }
    }
    return Some(result);
}

fn find_smallest_gcd_entry_in_pivot_row<M, T>(A: Matrix<M, T>) -> usize 
    where M: MatrixViewMut<T>, T: Integer
{
    find_min(0..A.col_count(), |col: &usize| 
        signed_gcd(A.at(0, 0).clone(), A.at(0, *col).clone())
    ).unwrap()
}

fn find_smallest_gcd_entry_in_pivot_col<M, T>(A: Matrix<M, T>) -> usize
    where M: MatrixViewMut<T>, T: Integer
{
    find_min(0..A.row_count(), |row: &usize| 
        signed_gcd(A.at(0, 0).clone(), A.at(*row, 0).clone())
    ).unwrap()
}

fn swap_pivot_entry_if_zero<M, N, K, T>(
    mut A: Matrix<M, T>,
    mut iL: Matrix<N, T>,
    mut iR: Matrix<K, T>,
    pivot: usize,
) -> bool 
    where M: MatrixViewMut<T>, N: MatrixViewMut<T>, K: MatrixViewMut<T>, T: Integer
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

fn find_not_zero<M, T>(mat: Matrix<M, T>) -> Option<(usize, usize)> 
    where M: MatrixView<T>, T: Integer
{
    for row in 0..mat.row_count() {
        for col in 0..mat.col_count() {
            if *mat.at(row, col) != T::zero() {
                return Some((row, col));
            }
        }
    }
    return None;
}

#[test]
fn test_diophantine() {
    let A = Matrix::from_array([[15, 10], [6, 7]]);
    let b = Vector::from_array([195, 87]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert_eq!(Vector::from_array([11, 3]), x.unwrap());
}

#[test]
fn test_diophantine_no_solution() {
    let A = Matrix::from_array([[2, -2]]);
    let b = Vector::from_array([1]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert!(x.is_none());
}

#[test]
fn test_diophantine_no_solution_three_dim() {
    #[rustfmt::skip]
    let A = Matrix::from_array([[1, 2, 0], 
                                [1, 0, 2]]);

    let b = Vector::from_array([2, 1]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert!(x.is_none());
}

#[test]
fn test_diophantine_three_dim() {
    #[rustfmt::skip]
    let A = Matrix::from_array([[1, 2, 0], 
                                [1, 0, 2]]);

    let b = Vector::from_array([2, 4]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert_eq!(Vector::from_array([4, -1, 0]), x.unwrap());
}

#[test]
fn test_diophantine_unnecessary_conditions() {
    #[rustfmt::skip]
    let A = Matrix::from_array([[1, 2, 0], 
                                [1, 2, 0], 
                                [1, 2, 0], 
                                [1, 0, 2]]);

    let b = Vector::from_array([2, 2, 2, 4]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert_eq!(Vector::from_array([4, -1, 0]), x.unwrap());
}

#[test]
fn test_diophantine_no_rational_solutions() {
    #[rustfmt::skip]
    let A = Matrix::from_array([[1, 2, 0], 
                                [1, 2, 0], 
                                [1, 0, 2]]);

    let b = Vector::from_array([2, 3, 4]);
    let x = diophantine_solve(A.as_ref(), b.as_ref());
    assert!(x.is_none());
}
