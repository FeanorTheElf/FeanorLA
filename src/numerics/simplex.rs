use super::super::la::mat::*;
use super::super::alg::*;
use std::ops::MulAssign;
use std::vec::Vec;

type BasicVars = Box<[usize]>;

#[derive(Debug, PartialEq)]
pub struct SystemUnbounded;

///
/// Optimize c^T x with x >= 0 and Ax=b
/// table: (0 | c^T)
///        (b |  A )
/// 
/// Returns Err if problem is unbounded,
/// otherwise table (0 | c^T)
///                 (b |  A )
/// with
///
fn simplex<M, T>(mut table: Matrix<M, T>, basic_vars: &mut BasicVars) -> Result<(), SystemUnbounded> 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone + 'static
{
    while let Some(pivot_col) = find_pivot_col(table.row(0)) {
        pivot(table.as_mut(), pivot_col, basic_vars)?;
    }
    return Ok(());
}

/*
 * Find solution of Ax = b with x >= 0
 * table: (b | A)
 */
pub fn solve<M, T>(table: Matrix<M, T>) -> Option<Vec<T>> 
    where M: MatrixView<T>, T: FieldEl + PartialOrd + Clone + 'static
{
    let (mut matrix, mut basic_vars) = add_artificials(table.as_ref());
    simplex(matrix.as_mut(), &mut basic_vars).unwrap();
    let solution = extract_solution(matrix, &basic_vars);

    let result: Vec<T> = Vec::from(&solution[0..(table.col_count() - 1)]);
    if is_solution(&result, table) {
        return Some(result);
    } else {
        return None;
    }
}

fn is_solution<M, T>(vars: &[T], table: Matrix<M, T>) -> bool 
    where M: MatrixView<T>, T: FieldEl + PartialOrd + Clone
{
    assert_eq!(
        vars.len() + 1,
        table.col_count(),
        "Expected one variable for each column except the first, got {} variables and {} columns",
        vars.len(),
        table.col_count()
    );
    for row_index in 0..table.row_count() {
        let mut current: T = T::zero();
        for var_index in 0..vars.len() {
            current += vars[var_index].clone() * table.at(row_index, var_index + 1).clone();
        }
        if current != *table.at(row_index, 0) {
            return false;
        }
    }
    return true;
}

fn extract_solution<M, T>(table: Matrix<M, T>, basic_vars: &BasicVars) -> Box<[T]> 
    where M: MatrixView<T>, T: FieldEl + Clone
{
    assert_eq!(basic_vars.len(), table.row_count() - 1);
    let mut result: Box<[T]> = {
        let mut vec = Vec::new();
        vec.resize(table.col_count(), T::zero());
        vec.into_boxed_slice()
    };
    for row_index in 0..basic_vars.len() {
        result[basic_vars[row_index] - 1] = table.at(row_index + 1, 0).clone();
        debug_assert!(*table.at(row_index + 1, basic_vars[row_index]) == T::one());
    }
    result[table.col_count() - 1] = -table.at(0, 0).clone();
    return result;
}

fn pivot<M, T>(
    table: Matrix<M, T>,
    pivot_col_index: usize,
    basic_vars: &mut BasicVars,
) -> Result<(), SystemUnbounded> 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone
{
    let pivot_row_index: usize = find_pivot_row(table.as_ref(), pivot_col_index)?;
    basic_vars[pivot_row_index - 1] = pivot_col_index;
    eliminate(table, pivot_row_index, pivot_col_index);
    return Ok(());
}

fn eliminate<M, T>(mut table: Matrix<M, T>, row_index: usize, col_index: usize) 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone
{
    let pivot_value: T = table.at(row_index, col_index).clone();
    assert!(pivot_value != T::zero());
    table.row_mut(row_index).mul_assign(T::one() / pivot_value);
    for target_row_index in 0..table.row_count() {
        if target_row_index < row_index {
            let factor = table.at(target_row_index, col_index).clone();
            table.transform_two_dims_left(target_row_index, row_index, &[T::one(), -factor, T::zero(), T::one()], &StaticRing::<T>::RING);
        } else if target_row_index > row_index {
            let factor = table.at(target_row_index, col_index).clone();
            table.transform_two_dims_left(row_index, target_row_index, &[T::one(), T::zero(), -factor, T::one()], &StaticRing::<T>::RING);
        }
    }
}

fn find_pivot_row<M, T>(table: Matrix<M, T>, pivot_col_index: usize) -> Result<usize, SystemUnbounded> 
    where M: MatrixView<T>, T: FieldEl + PartialOrd + Clone
{
    let last_col: usize = table.col_count() - 1;
    let mut current_min: Option<(usize, T)> = None;
    for row_index in 1..table.row_count() {
        if *table.at(row_index, pivot_col_index) > T::zero() {
            let row_value = table.at(row_index, last_col).clone() / table.at(row_index, pivot_col_index).clone();
            if current_min.as_ref().map_or(true, |(_index, min)| *min > row_value) {
                current_min = Some((row_index, row_value));
            }
        }
    }
    if let Some((result, _value)) = current_min {
        return Ok(result);
    } else {
        return Err(SystemUnbounded);
    }
}

fn find_pivot_col<V, T>(row: Vector<V, T>) -> Option<usize> 
    where V: VectorView<T>, T: FieldEl + PartialOrd
{
    for i in 0..row.len() {
        if *row.at(i) > T::zero() {
            return Some(i);
        }
    }
    return None;
}

fn add_artificials<M, T>(table: Matrix<M, T>) -> (Matrix<MatrixOwned<T>, T>, BasicVars) 
    where M: MatrixView<T>, T: FieldEl + Clone + PartialOrd, T: 'static
{
    let rows = table.row_count() + 1;
    let cols = table.col_count() + table.row_count();
    let mut basic_vars = {
        let mut vec = Vec::new();
        vec.resize(table.row_count(), 0);
        vec.into_boxed_slice()
    };
    let mut result: Matrix<MatrixOwned<T>, T> = Matrix::zero(rows, cols);
    for row_index in 1..rows {
        for col_index in 0..table.col_count() {
            *result.at_mut(row_index, col_index) = table.at(row_index - 1, col_index).clone();
        }
        if *result.at(row_index, 0) < T::zero() {
            result.row_mut(row_index).mul_assign(-T::one());
        }
        result.transform_two_dims_left(0, row_index, &[T::one(), T::one(), T::zero(), T::one()], &StaticRing::<T>::RING);
    }
    for row_index in 1..rows {
        let basic_var_col = table.col_count() + row_index - 1;
        *result.at_mut(row_index, basic_var_col) = T::one();
        basic_vars[row_index - 1] = basic_var_col;
    }
    *result.at_mut(0, 0) = T::zero();
    return (result, basic_vars);
}

#[cfg(test)]
use super::macros::ApproxEq;

#[test]
fn test_simplex_no_artificials() {
    let mut basic_vars: Box<[usize]> = Box::new([2, 3]);

    #[rustfmt::skip]
	let mut m = Matrix::from_array([[3.0, 4.0, 0.0, 0.0, 0.0,  2.0],
	                                [2.0, 1.0, 1.0, 0.0, 10.0, 3.0],
                                    [5.0, 3.0, 0.0, 1.0, 15.0, 2.0]]);

    assert_eq!(Ok(()), simplex(m.as_mut(), &mut basic_vars));

    #[rustfmt::skip]
	assert_approx_eq!(&Matrix::from_array([[-3.6666, 0.0, 0.0, -1.3333, -20.0, -0.6666],
	                                       [0.3333,  0.0, 1.0, -0.3333, 5.0,   2.3333 ], 
                                           [1.6666,  1.0, 0.0, 0.3333,  5.0,   0.6666 ]]), &m, 0.001);

    assert_eq!(&[2, 1], &*basic_vars);
}

#[test]
fn test_extract_solution() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[-11.0, 0.0, 0.0, -4.0, -60.0, -2.0],  
	                            [1.0,   0.0, 1.0, -1.0, 15.0,  7.0],  
				                [5.0,   1.0, 0.0, 1.0,  10.0,  2.0]]);
    let basic_vars: Box<[usize]> = Box::new([2, 1]);
    let solution = extract_solution(m.as_ref(), &basic_vars);
    assert_eq!(&[5.0, 1.0, 0.0, 0.0, 0.0, 11.0], &*solution);
}

#[test]
fn test_is_solution() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[1.0, 0.0, 1.0, -1.0, 15.0,  7.0],  
				                [5.0, 1.0, 0.0, 1.0,  10.0,  2.0]]);
    assert_eq!(true, is_solution(&[5.0, 1.0, 0.0, 0.0, 0.0], m.as_ref()));

    #[rustfmt::skip]
	let m = Matrix::from_array([[9.0, 3.0, 0.0, 0.0, 0.0,  2.0],
	                            [2.0, 1.0, 1.0, 0.0, 10.0, 3.0],
							    [5.0, 3.0, 0.0, 1.0, 15.0, 2.0]]);
    assert_eq!(true, is_solution(&[3.0, -1.0, -4.0, 0.0, 0.0], m.as_ref()));
    assert_eq!(false, is_solution(&[4.0, -1.0, -4.0, 0.0, 0.0], m.as_ref()));
}

#[test]
fn test_add_artificials() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[6.0,  5.0,  4.0], 
								[-9.0, 8.0,  7.0], 
								[12.0, 11.0, 10.0]]);
    let (result, basic_vars) = add_artificials(m.as_ref());

    #[rustfmt::skip]
	assert_eq!(&Matrix::from_array([[0.0,  8.0,  7.0,  0.0, 0.0, 0.0], 
	                                [6.0,  5.0,  4.0,  1.0, 0.0, 0.0], 
				                    [9.0,  -8.0, -7.0, 0.0, 1.0, 0.0], 
				                    [12.0, 11.0, 10.0, 0.0, 0.0, 1.0]]), &result);
    assert_eq!(&[3, 4, 5], &*basic_vars);
}

#[test]
fn test_solve() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[-1.0, -1.0, 0.0,  1.0, 0.0, 0.0], 
								[4.0,  1.0,  1.0,  0.0, 1.0, 0.0], 
								[0.0,  1.0,  -1.0, 0.0, 0.0, 1.0]]);
    let solution = solve(m.as_ref());
    assert!(is_solution(&solution.unwrap()[0..5], m.as_ref()));
}

#[test]
fn test_solve_zero_vec_solution() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[0.0, 1.0, 1.0, -1.0, 0.0],
	                            [0.0, 1.0, 0.0, -1.0, -1.0]]);
    assert_eq!(&[0.0, 0.0, 0.0, 0.0], &*solve(m.as_ref()).unwrap());
}

#[test]
fn test_impossible_system_solve() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[1.0, 1.0,  -1.0],
	                            [1.0, -1.0, 1.0]]);
    assert_eq!(None, solve(m.as_ref()));
}
