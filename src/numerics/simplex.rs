use super::super::la::mat::*;
use super::super::alg::*;
use std::ops::MulAssign;
use std::vec::Vec;

type BasicVars = Box<[usize]>;

#[derive(Debug, PartialEq)]
pub struct SystemUnbounded;

const CONSTANTS_COL: usize = 0;
const OBJECTIVE_ROW: usize = 0;

///
/// Minimize c^T x with x >= 0 and Ax = b
/// table: (0 | c^T)
///        (b |  A )
/// 
/// Returns Err if problem is unbounded,
/// otherwise table (0 | c^T)
///                 (b |  A )
/// together with the basic variables, such that
/// setting all variables not occuring in basic_vars
/// to zero and the variables occuring in basic_vars to
/// the value required for satisfying Ax = b is the
/// optimal solution.
/// 
/// Requires that for the i-th row except for the first one,
/// the basic_vars[i]-th column contains the (i + 1)-th unit
/// vector. The variables corresponding to these columns
/// are called basic variables, and are computed by linear
/// elimination after setting the other variables to zero.
/// 
/// The algorithm itself proceeds by changing the set of basic
/// variables and adjusting the table correspondingly.
///
fn simplex<M, T>(mut table: Matrix<M, T>, basic_vars: &mut BasicVars) -> Result<(), SystemUnbounded> 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone + 'static
{
    assert_eq!(table.row_count(), basic_vars.len() + 1);
    while let Some(pivot_col) = find_pivot_col(table.row(OBJECTIVE_ROW)) {
        pivot(table.as_mut(), pivot_col, basic_vars)?;
    }
    return Ok(());
}

///
/// Find solution of Ax = b with x >= 0
/// table: (b | A)
///
pub fn solve<M, T>(table: Matrix<M, T>) -> Option<Vector<VectorOwned<T>, T>> 
    where M: MatrixView<T>, T: FieldEl + PartialOrd + Clone + 'static + std::fmt::Display
{
    let (mut matrix, mut basic_vars) = add_artificials(table.as_ref());
    simplex(matrix.as_mut(), &mut basic_vars).unwrap();
    let (solution, _value) = extract_solution(matrix, &basic_vars);
    if is_solution(solution.as_ref(), table) {
        return Some(solution);
    } else {
        return None;
    }
}

fn is_solution<V, M, T>(vars: Vector<V, T>, table: Matrix<M, T>) -> bool 
    where V: VectorView<T>, M: MatrixView<T>, T: FieldEl + PartialOrd + Clone
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
        if current != *table.at(row_index, CONSTANTS_COL) {
            return false;
        }
    }
    return true;
}

fn extract_solution<M, T>(table: Matrix<M, T>, basic_vars: &BasicVars) -> (Vector<VectorOwned<T>, T>, T)
    where M: MatrixView<T>, T: FieldEl + Clone
{
    assert_eq!(basic_vars.len(), table.row_count() - 1);
    let mut result = Vector::zero(table.col_count() - 1).into_owned();
    for row_index in 0..basic_vars.len() {
        *result.at_mut(basic_vars[row_index] - 1) = table.at(row_index + 1, CONSTANTS_COL).clone();
        debug_assert!(*table.at(row_index + 1, basic_vars[row_index]) == T::one());
    }
    return (result, -table.at(OBJECTIVE_ROW, CONSTANTS_COL).clone());
}

fn pivot<M, T>(
    table: Matrix<M, T>,
    pivot_col_index: usize,
    basic_vars: &mut BasicVars,
) -> Result<(), SystemUnbounded> 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone
{
    // the variable with pivot column index replaces another
    // basic variable, i.e. we set the variable with index pivot
    // column index to zero
    let pivot_row_index: usize = find_pivot_row(table.as_ref(), pivot_col_index)?;
    basic_vars[pivot_row_index - 1] = pivot_col_index;
    eliminate(table, pivot_row_index, pivot_col_index);
    return Ok(());
}

///
/// Subtracts multiples of the row with given index from all other
/// rows, such that the only non-zero entry in the column with given
/// index is in this row.
/// 
/// Obviously requires that the entry at row_index, col_index is nonzero.
/// 
fn eliminate<M, T>(mut table: Matrix<M, T>, row_index: usize, col_index: usize) 
    where M: MatrixMutRowIter<T>, T: FieldEl + PartialOrd + Clone
{
    let pivot_value: T = table.at(row_index, col_index).clone();
    assert!(pivot_value != T::zero());
    table.row_mut(row_index).mul_assign(T::one() / pivot_value);
    for target_row_index in 0..table.row_count() {
        if target_row_index < row_index {
            let factor = table.at(target_row_index, col_index).clone();
            table.transform_two_dims_left(target_row_index, row_index, &[T::one(), -factor, T::zero(), T::one()], &T::RING);
        } else if target_row_index > row_index {
            let factor = table.at(target_row_index, col_index).clone();
            table.transform_two_dims_left(row_index, target_row_index, &[T::one(), T::zero(), -factor, T::one()], &T::RING);
        }
    }
}

///
/// Finds the row in the table that corresponds to the
/// constraint that prevents decreasing the variable with given
/// index most.
/// 
fn find_pivot_row<M, T>(table: Matrix<M, T>, pivot_col_index: usize) -> Result<usize, SystemUnbounded> 
    where M: MatrixView<T>, T: FieldEl + PartialOrd + Clone
{
    let mut current_min: Option<(usize, T)> = None;
    for row_index in 1..table.row_count() {
        if *table.at(row_index, pivot_col_index) > T::zero() {
            let row_value = table.at(row_index, CONSTANTS_COL).clone() / table.at(row_index, pivot_col_index).clone();
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

///
/// Finds the index of a variable such that decreasing it will
/// decrease the value of the objective function, given by 
/// c^T x where c is the given vector.
/// 
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

///
/// Creates a simplex tableau that contains the variables from the original
/// table without objective function (i.e. just (b | A)), and an additional
/// variable for each constraint, representing the deviation from that constraint.
/// Apart from that, the resulting table contains an objective function that
/// corresponds to minimizing the deviation from the constraints.
/// 
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
    let mut result: Matrix<MatrixOwned<T>, T> = Matrix::zero(rows, cols).into_owned();
    for row_index in 1..rows {
        for col_index in 0..table.col_count() {
            *result.at_mut(row_index, col_index) = table.at(row_index - 1, col_index).clone();
        }
        if *result.at(row_index, 0) < T::zero() {
            result.row_mut(row_index).mul_assign(-T::one());
        }
        result.transform_two_dims_left(0, row_index, &[T::one(), T::one(), T::zero(), T::one()], &T::RING);
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
    assert_eq!(Vector::from_array([5.0, 1.0, 0.0, 0.0, 0.0, 11.0]), solution);
}

#[test]
fn test_is_solution() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[1.0, 0.0, 1.0, -1.0, 15.0,  7.0],  
				                [5.0, 1.0, 0.0, 1.0,  10.0,  2.0]]);
    assert_eq!(true, is_solution(Vector::from_array([5.0, 1.0, 0.0, 0.0, 0.0]), m.as_ref()));

    #[rustfmt::skip]
	let m = Matrix::from_array([[9.0, 3.0, 0.0, 0.0, 0.0,  2.0],
	                            [2.0, 1.0, 1.0, 0.0, 10.0, 3.0],
							    [5.0, 3.0, 0.0, 1.0, 15.0, 2.0]]);
    assert_eq!(true, is_solution(Vector::from_array([3.0, -1.0, -4.0, 0.0, 0.0]), m.as_ref()));
    assert_eq!(false, is_solution(Vector::from_array([4.0, -1.0, -4.0, 0.0, 0.0]), m.as_ref()));
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
    assert!(is_solution(solution.unwrap(), m.as_ref()));
}

#[test]
fn test_solve_zero_vec_solution() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[0.0, 1.0, 1.0, -1.0, 0.0],
	                            [0.0, 1.0, 0.0, -1.0, -1.0]]);
    assert_eq!(Vector::from_array([0.0, 0.0, 0.0, 0.0]), solve(m.as_ref()).unwrap());
}

#[test]
fn test_impossible_system_solve() {
    #[rustfmt::skip]
	let m = Matrix::from_array([[1.0, 1.0,  -1.0],
	                            [1.0, -1.0, 1.0]]);
    assert_eq!(None, solve(m.as_ref()));
}

#[test]
fn test_fractional_solution() {
    #[rustfmt::skip]
    let m = Matrix::from_array([
        [                1.,                 1.,                 1.],
        [0.3333333333333333,                 1.,                 0.],
        [0.6666666666666666,                 0.,                 1.],
        [                1.,                 1.,                 1.]
    ]);
    assert!(solve(m.as_ref()).is_some());
}