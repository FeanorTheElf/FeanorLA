use super::mat::*;
use super::super::prelude::*;

use std::cmp::{min, max};

pub trait MatrixSolve<M: MatrixViewMut<Self::El>>: Ring {
    ///
    /// If you want to solve linear equations, prefer [`Matrix::solve_right_modifying()`].
    /// This is an implemenation tool to use the correct specializations for 
    /// different rings and matrix views.
    /// 
    /// Solves the linear equation AX = B where A is an square matrix.
    /// This is done by transforming A and B, and after this function
    /// successfully terminated, A will be in strict upper triangle form and
    /// B will contain the solution.
    /// 
    /// If A is not invertible, this function should terminate as soon
    /// as this has been detected. As part of the error object, the smallest
    /// integer i is returned for which the first i columns are linearly
    /// dependent.
    /// 
    /// The default implementation has complexity O(n^2(n + m)) where 
    /// self is nxn and rhs is nxm.
    /// 
    fn solve_right_modifying<N: MatrixViewMut<Self::El>>(&self, a: Matrix<M, Self::El>, b: &mut Matrix<N, Self::El>) -> Result<(), usize>;

    ///
    /// If you want to compute a kernel basis, prefer [`Matrix::right_kernel_base_modifying()`].
    /// This is an implemenation tool to use the correct specializations for 
    /// different rings and matrix views.
    /// 
    /// Calculates a base of the (right)-kernel of this matrix, or returns None 
    /// if this kernel is trivial. These are all vectors x such that A * x = 0.
    /// 
    /// The default implementation has complexity O(n m min(n, m)) where self is nxm.
    /// 
    fn right_kernel_base_modifying(&self, a: Matrix<M, Self::El>) -> Option<Matrix<MatrixOwned<Self::El>, Self::El>>;

    ///
    /// If you want to compute a matrix rank, prefer [`Matrix::rank()`].
    /// This is an implemenation tool to use the correct specializations for 
    /// different rings and matrix views.
    /// 
    /// Calculates the rank of this matrix.
    /// 
    /// The default implementation has complexity O(n m min(n, m)) where self is nxm.
    /// 
    fn matrix_rank_modifying(&self, a: Matrix<M, Self::El>) -> usize;

    ///
    /// If you want to solve linear equations, prefer [`Matrix::find_any_solution_modifying()`].
    /// This is an implemenation tool to use the correct specializations for 
    /// different rings and matrix views.
    /// 
    /// Finds a solution to the inhomogeneous linear equation AX = B and returns it,
    /// or None if no solution exists.
    /// The returned solution is not unique, but you can find all solutions to the system
    /// by considering the returned solution plus the kernel of A (use [`right_kernel_base_modifying()`]).
    /// 
    fn find_any_solution_modifying<N: MatrixView<Self::El>>(&self, a: Matrix<M, Self::El>, b: Matrix<N, Self::El>) -> Result<Matrix<Submatrix<MatrixOwned<El<Self>>, El<Self>>, El<Self>>, ()>;
}

impl<R: Ring, M: MatrixViewMut<Self::El>> MatrixSolve<M> for R {

    default fn solve_right_modifying<N: MatrixViewMut<Self::El>>(&self, mut lhs: Matrix<M, Self::El>, rhs: &mut Matrix<N, Self::El>) -> Result<(), usize> {
        assert!(self.is_field().can_use());
        assert_eq!(lhs.row_count(), lhs.col_count());
        assert_eq!(lhs.row_count(), rhs.row_count());
        lhs.gaussion_elimination_half(
            |row, a, rhs| rhs.submatrix_mut(row..=row, 0..).scale(&a, self), 
            |i, j, rhs| rhs.swap_rows(i, j), 
            |dst, a, src, rhs| rhs.transform_two_dims_left(
                src, dst, &[self.one(), self.zero(), self.neg(a), self.one()], self
            ), 
            rhs, 
            self)?;

        lhs.solve_strict_triangular(rhs, self);

        return Ok(());
    }
    
    default fn right_kernel_base_modifying(&self, mut a: Matrix<M, Self::El>) -> Option<Matrix<MatrixOwned<Self::El>, Self::El>> {
        assert!(self.is_field().can_use());

        // the approach is to transform the matrix in upper triangle form, 
        // so ( U | R ) with an upper triangle matrix U and a nonsquare 
        // rest matrix R. Then the kernel base matrix
        // is given by ( -inv(U)*R )
        //             (     I     )
        // we just have to watch out if the left side is singular, then swap cols

        let mut col_swaps: Vec<(usize, usize)> = Vec::new();
        let non_zero_row_count = upper_trapezoid_form(a.as_mut(), self, |c1, c2| col_swaps.push((c1, c2)), (&mut (), |_, _, _| {}, |_, _, _| {}, |_, _, _, _| {}));

        // now a is in upper triangle form
        let effective_matrix = a.submatrix_mut(0..non_zero_row_count, ..);
        if effective_matrix.row_count() >= effective_matrix.col_count() {
            return None;
        }
        let upper_part = effective_matrix.row_count();
        let lower_part = effective_matrix.col_count() - effective_matrix.row_count();
        let mut result = Matrix::zero_ring(upper_part + lower_part, lower_part, self).into_owned();

        // set to identity in the lower part
        for i in 0..lower_part {
            *result.at_mut(upper_part + i, i) = self.neg(self.one());
        }
        
        // set the interesting upper part
        let mut result_upper_part = result.submatrix_mut(..upper_part, ..);
        result_upper_part.add_assign(effective_matrix.submatrix(.., upper_part..), self);
        effective_matrix.submatrix(.., ..upper_part).solve_strict_triangular(
            &mut result_upper_part, self
        );

        // and now perform the swaps
        for (row1, row2) in col_swaps.iter().rev() {
            result.swap_rows(*row1, *row2);
        }

        return Some(result);
    }

    default fn matrix_rank_modifying(&self, a: Matrix<M, Self::El>) -> usize {
        assert!(self.is_field().can_use());
        upper_trapezoid_form(a, self, |_, _| {}, (&mut (), |_, _, _| {}, |_, _, _| {}, |_, _, _, _| {}))
    }

    default fn find_any_solution_modifying<N>(&self, mut a: Matrix<M, Self::El>, b: Matrix<N, Self::El>) -> Result<Matrix<Submatrix<MatrixOwned<El<Self>>, El<Self>>, El<Self>>, ()> 
        where N: MatrixView<Self::El>
    {
        assert!(self.is_field().can_use());
        assert!(a.row_count() == b.row_count());
        let mut col_swaps: Vec<(usize, usize)> = Vec::new();

        let mut result = Matrix::zero_ring(max(b.row_count(), a.col_count()), b.col_count(), self).into_owned();
        result.submatrix_mut(..b.row_count(), ..).assign(b.as_ref());

        let gaussian_elim_ops = (
            &mut result.submatrix_mut(..b.row_count(), ..),
            |row, a, rhs: &mut Matrix<Submatrix<&mut MatrixOwned<El<Self>>, El<Self>>, El<Self>>| rhs.submatrix_mut(row..=row, 0..).scale(&a, self), 
            |i, j, rhs: &mut Matrix<Submatrix<&mut MatrixOwned<El<Self>>, El<Self>>, El<Self>>| rhs.swap_rows(i, j), 
            |dst, a, src, rhs: &mut Matrix<Submatrix<&mut MatrixOwned<El<Self>>, El<Self>>, El<Self>>| rhs.transform_two_dims_left(
                src, dst, &[self.one(), self.zero(), self.neg(a), self.one()], self
            ), 
        );
        let non_zero_row_count = upper_trapezoid_form(a.as_mut(), self, |c1, c2| col_swaps.push((c1, c2)), gaussian_elim_ops);

        if !result.submatrix(non_zero_row_count..b.row_count(), ..).eq(Matrix::zero_ring(b.row_count() - non_zero_row_count, b.col_count(), self), self) {
            return Err(());
        }
        a.submatrix(..non_zero_row_count, ..non_zero_row_count).solve_strict_triangular(&mut result.submatrix_mut(..non_zero_row_count, ..), self);
        
        return Ok(result.into_submatrix(..a.col_count(), ..));
    }
}

///
/// Using row operations and column swaps, transforms a matrix A into
/// what we call "upper trapezoid form", which means that except for
/// zero rows at the bottom, the matrix is in upper triangle form with
/// non-zero entries on the diagonal. In other words, it looks like
/// ```text
/// [ 1 * * ... * ... * ]
/// [   1 * ... * ... * ]
///   ...
/// [           1 ... * ]
/// [                   ]
///   ...
/// [                   ]
/// ```
/// where `*` means any entry.
/// 
/// Returns the number of non-zero rows after the transformation.
/// 
fn upper_trapezoid_form<R, M, F, G1, G2, G3, S>(mut a: Matrix<M, El<R>>, ring: &R, mut col_swap: F, gaussian_elim_data: (&mut S, G1, G2, G3)) -> usize 
    where F: FnMut(usize, usize), 
        R: Ring, 
        M: MatrixViewMut<El<R>>,
        G1: Copy + FnMut(usize, El<R>, &mut S), 
        G2: Copy + FnMut(usize, usize, &mut S), 
        G3: Copy + FnMut(usize, El<R>, usize, &mut S),
{
    assert!(ring.is_field().can_use());
    let mut i = 0;
    loop {
        let mut current_submatrix = a.submatrix_mut(i.., i..);
        let gaussian_elim_result = current_submatrix.gaussion_elimination_half(
            gaussian_elim_data.1, gaussian_elim_data.2, gaussian_elim_data.3, gaussian_elim_data.0, ring
        );
        let (col1, col2) = if let Err(null_col) = gaussian_elim_result {
            if let Some(other_col) = find_non_null_column(
                a.submatrix((null_col + i).., (null_col + i)..), ring
            ) {
                // swap columns
                (null_col + i, null_col + i + other_col)
            } else {
                // we have a whole 0-rectangle in the lower right corner, 
                // so we really are in upper triangle form
                return null_col + i;
            }
        } else {
            // upper triangle form is reached
            return min(a.col_count(), a.row_count());
        };
        col_swap(col1, col2);
        a.swap_cols(col1, col2);
        i = col1;
    }
}

fn find_non_null_column<R, N, T>(matrix: Matrix<N, T>, ring: &R) -> Option<usize> 
    where N: MatrixView<T>, T: std::fmt::Debug + Clone, R: Ring<El = T>
{
    for i in 0..matrix.col_count() {
        for j in 0..matrix.row_count() {
            if !ring.is_zero(matrix.at(j, i)) {
                return Some(i);
            }
        }
    }
    return None;
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>, T: std::fmt::Debug + Clone
{
    ///
    /// Transforms this matrix into strict upper triangle form by gaussian elimination.
    /// If the matrix is singular, the process is aborted as soon as this is detected,
    /// and an index i is returned such that the first i columns of this matrix are in
    /// upper triangle form, with the (i,i)-th entry being zero.
    /// 
    /// The passed functions are called whenever a gaussian elimination step is made,
    /// so that the caller can adjust other data as well (i.e. a potential 
    /// right-hand side of a linear equation, or a determinant)
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n m min(n,m))
    /// 
    pub fn gaussion_elimination_half<R, F, G, H, S>(
        &mut self, 
        mut mul_row: F, 
        mut swap_rows: G, 
        mut sub_row: H, 
        state: &mut S,
        ring: &R
    ) -> Result<(), usize>
        where F: FnMut(usize, T, &mut S), 
            G: FnMut(usize, usize, &mut S), 
            H: FnMut(usize, T, usize, &mut S),
            R: Ring<El = T>
    {
        assert!(ring.is_field().can_use());

        for i in 0..min(self.col_count(), self.row_count()) {
            // pivot
            if ring.is_zero(self.at(i, i)) {
                let mut has_swapped = false;
                for j in (i + 1)..self.row_count() {
                    if !ring.is_zero(self.at(j, i)) {
                        has_swapped = true;
                        self.swap_rows(i, j);
                        swap_rows(i, j, state);
                        break;
                    }
                }
                if !has_swapped {
                    return Err(i);
                }
            }

            // normalize pivot to 1
            let inverse = ring.div(ring.one(), self.at(i, i));
            self.submatrix_mut(i..=i, i..).scale(&inverse, ring);
            mul_row(i, inverse, state);

            // eliminate
            for j in (i + 1)..self.row_count() {
                let transform = [
                    ring.one(), 
                    ring.zero(), 
                    ring.neg(self.at(j, i).clone()), 
                    ring.one()
                ];
                sub_row(j, self.at(j, i).clone(), i, state);
                self.transform_two_dims_left(i, j, &transform, ring);
            }
        }
        return Ok(());
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: std::fmt::Debug + Clone
{
    ///
    /// Expects self to be a square, strict upper triangular matrix (i.e. 1 on the diagonal)
    /// and assigns to rhs the solution of the equation self * X = rhs
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn solve_strict_triangular<N, R>(&self, rhs: &mut Matrix<N, T>, ring: &R)
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        assert!(ring.is_field().can_use());
        assert_eq!(self.row_count(), self.col_count());
        assert_eq!(self.row_count(), rhs.row_count());

        // check upper triangle
        #[cfg(debug)] {
            for i in 0..self.row_count() {
                assert!(ring.is_one(self.at(i, i)));
                for j in (0..i) {
                    assert!(ring.is_zero(self.at(i, j)));
                }
            }
        }

        // now self is in upper triangle form
        for col in 0..rhs.col_count() {
            for row in (0..self.row_count()).rev() {
                for i in (row + 1)..self.row_count() {
                    let d = ring.mul_ref(self.at(row, i), rhs.at(i, col));
                    ring.add_assign(rhs.at_mut(row, col), ring.neg(d));
                }
            }
        }
    }

    ///
    /// Solves the linear equation AX = B where A is an square matrix.
    /// This is done by transforming A into echelon form.
    /// 
    /// If A is not invertible, this function will terminate as soon
    /// as this has been detected. Then A will consist of a left part
    /// in upper triangle form and some right part, where the (i-1)-th
    /// column and the i-th column both have at most the first (i-1) entries
    /// unequal to zero. i is returned as part of the error object. Therefore,
    /// i is the smallest integer so that the first i columns are linearly
    /// dependent
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn solve_right<R, N>(self, rhs: &mut Matrix<N, T>, ring: &R) -> Result<(), usize>
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        <R as MatrixSolve<MatrixOwned<T>>>::solve_right_modifying(ring, self.into_owned(), rhs)
    }

    ///
    /// Calculates a base of the kernel of this matrix, or returns None 
    /// if this kernel is trivial.
    /// 
    /// Complexity O(n^3) where self is nxn
    /// 
    pub fn right_kernel_base<R>(self, ring: &R) -> Option<Matrix<MatrixOwned<T>, T>> 
        where R: Ring<El = T>
    {
        <R as MatrixSolve<MatrixOwned<T>>>::right_kernel_base_modifying(ring, self.into_owned())
    }

    ///
    /// Computes the rank of this matrix.
    /// 
    /// Complexity O(n^3) where self is nxn
    /// 
    pub fn rank<R>(self, ring: &R) -> usize 
        where R: Ring<El = T>
    {
        <R as MatrixSolve<MatrixOwned<T>>>::matrix_rank_modifying(ring, self.into_owned())
    }

    ///
    /// Finds a solution to the inhomogeneous equation AX = B, or returns `None` if not such
    /// solution exists. As opposed to [`solve_right()`], this also works for non-invertible
    /// (and non-square) matrices, but then the returned solution is non-unique. No guarantee
    /// is given which of the possible solutions is returned. You can find all solutions by
    /// combining the result with [`right_kernel_base()`].
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn find_any_solution<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Result<Matrix<Submatrix<MatrixOwned<T>, T>, T>, ()>
        where N: MatrixView<T>, R: Ring<El = T>
    {
        <R as MatrixSolve<MatrixOwned<T>>>::find_any_solution_modifying(ring, self.into_owned(), rhs)
    }
}

#[cfg(test)]
use super::super::rational::primitive_rational::*;

#[test]
fn test_invert_matrix() {
    let a = Matrix::from_array([[1., 2.], [4., 8.]]);
    let mut a_inv = Matrix::identity(2, 2);
    assert!(a.solve_right(&mut a_inv, &f32::RING).is_err());

    let b = Matrix::from_array([[1., 2.], [4., 4.]]);
    let mut b_inv = Matrix::identity(2, 2);
    b.solve_right(&mut b_inv, &f32::RING).unwrap();

    assert_eq!(1., *b_inv.at(1, 0));
    assert_eq!(-0.25, *b_inv.at(1, 1));
    assert_eq!(-1., *b_inv.at(0, 0));
    assert_eq!(0.5, *b_inv.at(0, 1));
}

#[test]
fn test_kernel_base() {
    #[rustfmt::skip]
    let a = Matrix::from_array([[1., 3., 1., 3.], 
                                [2., 6., 1., 2.], 
                                [0., 0., 1., 1.]]);

    let b = Matrix::from_array([[3.], [-1.], [0.], [0.]]);
    assert_eq!(b, a.right_kernel_base(&f32::RING).unwrap());

    #[rustfmt::skip]
    let a = Matrix::from_array([[1., 3., 1., 3.], 
                                [2., 6., 1., 5.], 
                                [0., 0., 1., 1.]]);
                                    
    #[rustfmt::skip]
    let b = Matrix::from_array([[3.,  2.], 
                                [-1., 0.], 
                                [0.,  1.], 
                                [0., -1.]]);
    assert_eq!(b, a.right_kernel_base(&f32::RING).unwrap());
}

#[test]
fn test_find_any_solution() {
    #[rustfmt::skip]
    let a = Matrix::from_array([[1.,  2.],
                                [3.,  0.],
                                [-2., 1.]]);
    let b = Matrix::col_vec(Vector::from_array([-3., 3., -4.]));
    let expected = Matrix::col_vec(Vector::from_array([1., -2.]));
    assert_eq!(expected, a.find_any_solution(b, &f32::RING).unwrap());

    let i = r64::RING.embedding();
    #[rustfmt::skip]
    let a = Matrix::map(Matrix::from_array([[1,  2,  1],
                                            [3,  0,  2],
                                            [4,  2,  3],
                                            [-2, 1,  0]]), i);
    let b = Matrix::col_vec(Vector::from_array([-2, 5, 3, -4]).map(i).compute());
    let expected = Matrix::col_vec(Vector::from_array([1, -2, 1]).map(i).compute());
    assert_eq!(expected, a.find_any_solution(b, &r64::RING).unwrap());

    #[rustfmt::skip]
    let a = Matrix::from_array([[1., 2., 3., 1.],
                                [3., 0., 0., 2.]]);
    let b = Matrix::col_vec(Vector::from_array([-3., 3.]));
    let expected = Matrix::col_vec(Vector::from_array([1., -2., 0., 0.]));
    assert_eq!(expected, a.find_any_solution(b, &f32::RING).unwrap());
    
    #[rustfmt::skip]
    let a = Matrix::map(Matrix::from_array([[1,  2,  1],
                                            [3,  0,  2],
                                            [4,  2,  3],
                                            [-2, 1,  0]]), i);
    let b = Matrix::col_vec(Vector::from_array([-2, 5, 2, -1]).map(i).compute());
    assert_eq!(Err(()), a.find_any_solution(b, &r64::RING));
}
