use super::indexed::*;
use super::matrix::Matrix;
use std::ops::MulAssign;
use std::vec::Vec;

type Tableau = Matrix<f64>;
type TableauRow = [f64];
type BasicVars = Box<[usize]>;

#[derive(Debug, PartialEq)]
pub struct SystemUnbounded;

/*
 * Optimize c^T x with x >= 0 and Ax=b
 * table: (0 | c^T)
 *        (b |  A )
 *
 * Returns Err if problem is unbounded,
 * otherwise table (0 | c^T)
 *                 (b |  A )
 * with
 */
pub fn simplex(table: &mut Tableau, basic_vars: &mut BasicVars) -> Result<(), SystemUnbounded> {
    while let Some(pivot_col) = find_pivot_col(&table[0]) {
        pivot(table, pivot_col, basic_vars)?;
    }
    return Ok(());
}

/*
 * Find solution of Ax=b with x >= 0
 * table: (b | A)
 */
#[allow(unused)]
pub fn solve(table: &Tableau) -> Option<Vec<f64>> {
    let (mut matrix, mut basic_vars) = add_artificials(table);
    simplex(&mut matrix, &mut basic_vars).unwrap();
    let solution = extract_solution(&matrix, &basic_vars);

    let result: Vec<f64> = Vec::from(&solution[0..(table.cols() - 1)]);
    if is_solution(&result, table) {
        return Some(result);
    } else {
        return None;
    }
}

fn is_solution(vars: &[f64], table: &Tableau) -> bool {
    assert_eq!(
        vars.len() + 1,
        table.cols(),
        "Expected one variable for each column except the first, got {} variables and {} columns",
        vars.len(),
        table.cols()
    );
    for row_index in 0..table.rows() {
        let mut current: f64 = 0.0;
        for var_index in 0..vars.len() {
            current += vars[var_index] * table[row_index][var_index + 1];
        }
        if current != table[row_index][0] {
            return false;
        }
    }
    return true;
}

fn extract_solution(table: &Tableau, basic_vars: &BasicVars) -> Box<[f64]> {
    assert_eq!(basic_vars.len(), table.rows() - 1);
    let mut result: Box<[f64]> = {
        let mut vec = Vec::new();
        vec.resize(table.cols(), 0.0);
        vec.into_boxed_slice()
    };
    for row_index in 0..basic_vars.len() {
        result[basic_vars[row_index] - 1] = table[row_index + 1][0];
        debug_assert!(table[row_index + 1][basic_vars[row_index]] == 1.0);
    }
    result[table.cols() - 1] = -table[0][0];
    return result;
}

fn pivot(
    table: &mut Tableau,
    pivot_col_index: usize,
    basic_vars: &mut BasicVars,
) -> Result<(), SystemUnbounded> {
    let pivot_row_index: usize = find_pivot_row(table, pivot_col_index)?;
    basic_vars[pivot_row_index - 1] = pivot_col_index;
    eliminate(table, pivot_row_index, pivot_col_index);
    return Ok(());
}

fn eliminate(table: &mut Tableau, row_index: usize, col_index: usize) {
    let pivot_value: f64 = table[row_index][col_index];
    assert!(pivot_value != 0.0);
    table.get_mut(row_index).mul_assign(1.0 / pivot_value);
    for target_row_index in 0..table.rows() {
        if target_row_index != row_index {
            let factor = table[target_row_index][col_index];
            let (mut target_row, base_row) = table.get_rows(target_row_index, row_index);
            target_row.add_product(base_row.as_const(), -factor);
        }
    }
}

fn find_pivot_row(table: &Tableau, pivot_col_index: usize) -> Result<usize, SystemUnbounded> {
    let last_col: usize = table.cols() - 1;
    let mut current_min: Option<(usize, f64)> = None;
    for row_index in 1..table.rows() {
        if table[row_index][pivot_col_index] > 0.0 {
            let row_value = table[row_index][last_col] / table[row_index][pivot_col_index];
            if current_min.map_or(true, |(_index, min)| min > row_value) {
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

fn find_pivot_col(row: &TableauRow) -> Option<usize> {
    for i in 0..row.len() {
        if row[i] > 0.0 {
            return Some(i);
        }
    }
    return None;
}

fn add_artificials(table: &Tableau) -> (Tableau, BasicVars) {
    let rows = table.rows() + 1;
    let cols = table.cols() + table.rows();
    let mut basic_vars = {
        let mut vec = Vec::new();
        vec.resize(table.rows(), 0);
        vec.into_boxed_slice()
    };
    let mut result: Tableau = Tableau::zero(rows, cols);
    for row_index in 1..rows {
        for col_index in 0..table.cols() {
            result[row_index][col_index] = table[row_index - 1][col_index];
        }
        if result[row_index][0] < 0.0 {
            result.get_mut(row_index).mul_assign(-1.0);
        }
        let (mut dst_row, cnt_row) = result.get_rows(0, row_index);
        dst_row.add_product(cnt_row.as_const(), 1.0);
    }
    for row_index in 1..rows {
        let basic_var_col = table.cols() + row_index - 1;
        result[row_index][basic_var_col] = 1.0;
        basic_vars[row_index - 1] = basic_var_col;
    }
    result[0][0] = 0.0;
    return (result, basic_vars);
}

#[test]
fn test_simplex_no_artificials() {
    let mut basic_vars: Box<[usize]> = Box::new([2, 3]);

    #[rustfmt::skip]
	let mut m = Matrix::new(Box::new([3.0, 4.0, 0.0, 0.0, 0.0,  2.0,
	                                  2.0, 1.0, 1.0, 0.0, 10.0, 3.0,
                                      5.0, 3.0, 0.0, 1.0, 15.0, 2.0]), 3);

    assert_eq!(Ok(()), simplex(&mut m, &mut basic_vars));

    #[rustfmt::skip]
	assert_approx_eq!(&[-3.6666, 0.0, 0.0, -1.3333, -20.0, -0.6666,
	                    0.3333,  0.0, 1.0, -0.3333, 5.0,   2.3333, 
                        1.6666,  1.0, 0.0, 0.3333,  5.0,   0.6666], m.data(), 0.001);

    assert_eq!(&[2, 1], &*basic_vars);
}

#[test]
fn test_extract_solution() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([-11.0, 0.0, 0.0, -4.0, -60.0, -2.0,  
	                              1.0,   0.0, 1.0, -1.0, 15.0,  7.0,  
				                  5.0,   1.0, 0.0, 1.0,  10.0,  2.0]), 3);
    let basic_vars: Box<[usize]> = Box::new([2, 1]);
    let solution = extract_solution(&m, &basic_vars);
    assert_eq!(&[5.0, 1.0, 0.0, 0.0, 0.0, 11.0], &*solution);
}

#[test]
fn test_is_solution() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([1.0, 0.0, 1.0, -1.0, 15.0,  7.0,  
				                  5.0, 1.0, 0.0, 1.0,  10.0,  2.0]), 2);
    assert_eq!(true, is_solution(&[5.0, 1.0, 0.0, 0.0, 0.0], &m));

    #[rustfmt::skip]
	let m = Matrix::new(Box::new([9.0, 3.0, 0.0, 0.0, 0.0,  2.0,
	                              2.0, 1.0, 1.0, 0.0, 10.0, 3.0,
							      5.0, 3.0, 0.0, 1.0, 15.0, 2.0]), 3);
    assert_eq!(true, is_solution(&[3.0, -1.0, -4.0, 0.0, 0.0], &m));
    assert_eq!(false, is_solution(&[4.0, -1.0, -4.0, 0.0, 0.0], &m));
}

#[test]
fn test_add_artificials() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([6.0,  5.0,  4.0, 
								  -9.0, 8.0,  7.0, 
								  12.0, 11.0, 10.0]), 3);
    let (result, basic_vars) = add_artificials(&m);

    #[rustfmt::skip]
	assert_eq!(&[0.0,  8.0,  7.0,  0.0, 0.0, 0.0, 
	             6.0,  5.0,  4.0,  1.0, 0.0, 0.0, 
				 9.0,  -8.0, -7.0, 0.0, 1.0, 0.0, 
				 12.0, 11.0, 10.0, 0.0, 0.0, 1.0, ], result.data());
    assert_eq!(&[3, 4, 5], &*basic_vars);
}

#[test]
fn test_solve() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([-1.0, -1.0, 0.0,  1.0, 0.0, 0.0, 
								  4.0,  1.0,  1.0,  0.0, 1.0, 0.0, 
								  0.0,  1.0,  -1.0, 0.0, 0.0, 1.0]), 3);
    let solution = solve(&m);
    assert!(is_solution(&solution.unwrap()[0..5], &m));
}

#[test]
fn test_solve_zero_vec_solution() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([0.0, 1.0, 1.0, -1.0, 0.0,
	                              0.0, 1.0, 0.0, -1.0, -1.0]), 2);
    assert_eq!(&[0.0, 0.0, 0.0, 0.0], &*solve(&m).unwrap());
}

#[test]
fn test_impossible_system_solve() {
    #[rustfmt::skip]
	let m = Matrix::new(Box::new([1.0, 1.0,  -1.0,
	                              1.0, -1.0, 1.0]), 2);
    assert_eq!(None, solve(&m));
}
