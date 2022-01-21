use super::mat::*;
use super::super::alg::*;

pub trait MatrixSolve: Ring {
    ///
    /// Solves the linear equation AX = B where A is an square matrix.
    /// This is done by transforming A and B, and after this function
    /// successfully terminated, A will be in strict upper triangle form and
    /// B will contain the solution.
    /// 
    /// If A is not invertible, this function will terminate as soon
    /// as this has been detected. As part of the error object, the smallest
    /// integer i is returned for which the first i columns are linearly
    /// dependent
    /// 
    /// The default implementation has complexity O(n^2(n + m)) where 
    /// self is nxn and rhs is nxm.
    /// 
    fn solve_linear_equation<M: MatrixView<Self::El>, N: MatrixViewMut<Self::El>>(&self, a: Matrix<M, Self::El>, b: &mut Matrix<N, Self::El>) -> Result<(), usize>;
}

impl<R: Ring> MatrixSolve for R {

    default fn solve_linear_equation<M: MatrixView<Self::El>, N: MatrixViewMut<Self::El>>(&self, a: Matrix<M, Self::El>, b: &mut Matrix<N, Self::El>) -> Result<(), usize> {
        a.into_owned().solve_modifying(b, self)
    }
}

pub trait MatrixKernelBase: Ring {
    ///
    /// Calculates a base of the (right)-kernel of this matrix, or returns None 
    /// if this kernel is trivial. These are all vectors x such that A * x = 0.
    /// 
    /// The default implementation has complexity O(n m min(n, m)) where self is nxm.
    /// 
    fn calc_matrix_kernel_space<M: MatrixView<Self::El>>(&self, a: Matrix<M, Self::El>) -> Option<Matrix<MatrixOwned<Self::El>, Self::El>>;
}

impl<R: Ring> MatrixKernelBase for R {
    
    default fn calc_matrix_kernel_space<M: MatrixView<Self::El>>(&self, a: Matrix<M, Self::El>) -> Option<Matrix<MatrixOwned<Self::El>, Self::El>> {
        a.into_owned().kernel_base_modifying(self)
    }
}

pub trait MatrixFrobenius: Ring {
    fn calc_matrix_frobenius_norm_square<M: MatrixView<Self::El>>(&self, a: Matrix<M, Self::El>) -> Self::El;

    fn l2_norm_square<V: VectorView<Self::El>>(&self, v: Vector<V, Self::El>) -> Self::El {
        self.calc_matrix_frobenius_norm_square(Matrix::row_vec(v))
    }
}

impl<R: Ring> MatrixFrobenius for R {
    
    default fn calc_matrix_frobenius_norm_square<M: MatrixView<Self::El>>(&self, a: Matrix<M, Self::El>) -> Self::El
    {
        let mut it = a.rows().flat_map(|r| 
            (0..r.len()).map(move |i| self.mul_ref(r.at(i), r.at(i)))
        );
        let initial = it.next().unwrap();
        return it.fold(initial, |a, b| self.add(a, b));
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixView<T>, T: Clone + std::fmt::Debug
{
    ///
    /// Expects self to be a square, strict upper triangular matrix (i.e. 1 on the diagonal)
    /// and assigns to rhs the solution of the equation self * X = rhs
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    fn solve_strict_triangular<N, R>(&self, rhs: &mut Matrix<N, T>, ring: &R)
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
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
                    take_mut::take_or_recover(
                        rhs.at_mut(row, col), 
                        || ring.unspecified_element(), 
                        |v| ring.sub(v, d)
                    );
                }
            }
        }
    }
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
        debug_assert!(ring.is_field());

        for i in 0..std::cmp::min(self.col_count(), self.row_count()) {
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

    ///
    /// Solves the linear equation AX = B where A is an square matrix.
    /// This is done by transforming A and B, and after this function
    /// successfully terminated, A will be in strict upper triangle form and
    /// B will contain the solution.
    /// 
    /// If you do not need that this function requires no copy of the
    /// matrix, prefer solve() instead.
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
    fn solve_modifying<R, N>(&mut self, rhs: &mut Matrix<N, T>, ring: &R) -> Result<(), usize>
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        debug_assert!(ring.is_field());
        assert_eq!(self.row_count(), self.col_count());
        assert_eq!(self.row_count(), rhs.row_count());
        self.gaussion_elimination_half(
            |row, a, rhs| rhs.submatrix_mut(row..=row, 0..).scale(&a, ring), 
            |i, j, rhs| rhs.swap_rows(i, j), 
            |dst, a, src, rhs| rhs.transform_two_dims_left(
                src, dst, &[ring.one(), ring.zero(), ring.neg(a), ring.one()], ring
            ), 
            rhs, 
            ring)?;

        self.solve_strict_triangular(rhs, ring);

        return Ok(());
    }

    ///
    /// Calculates a base of the kernel of this matrix, or returns None 
    /// if this kernel is trivial. Note that this function modifies self, 
    /// so if you can life with an additional copy, prefer
    /// to use kernel_base() instead
    /// 
    fn kernel_base_modifying<R>(&mut self, ring: &R) -> Option<Matrix<MatrixOwned<T>, T>> 
        where R: Ring<El = T>
    {
        debug_assert!(ring.is_field());

        // the approach is to transform the matrix in upper triangle form, 
        // so ( U | R ) with an upper triangle matrix U and a nonsquare 
        // rest matrix R. Then the kernel base matrix
        // is given by ( -inv(U)*R )
        //             (     I     )
        // we just have to watch out if the left side is singular, then swap cols
        let mut col_swaps: Vec<(usize, usize)> = Vec::new();

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

        let mut i = 0;
        let mut non_zero_row_count = self.row_count();
        loop {
            let mut current_submatrix = self.submatrix_mut(i.., i..);
            let gaussian_elim_result = current_submatrix.gaussion_elimination_half(
                |_, _, _| {}, |_, _, _| {}, |_, _, _, _| {}, &mut (), ring
            );
            let (col1, col2) = if let Err(null_col) = gaussian_elim_result {
                if let Some(other_col) = find_non_null_column(
                    self.submatrix((null_col + i).., (null_col + i)..), ring
                ) {
                    // swap columns
                    (null_col + i, null_col + i + other_col)
                } else {
                    // we have a whole 0-rectangle in the lower right corner, 
                    // so we really are in upper triangle form
                    non_zero_row_count = null_col + i;
                    break;
                }
            } else {
                // upper triangle form is reached
                break;
            };
            col_swaps.push((col1, col2));
            self.swap_cols(col1, col2);
            i = col1;
        }

        // now self is in upper triangle form
        let effective_matrix = self.submatrix_mut(0..non_zero_row_count, ..);
        if effective_matrix.row_count() >= effective_matrix.col_count() {
            return None;
        }
        let upper_part = effective_matrix.row_count();
        let lower_part = effective_matrix.col_count() - effective_matrix.row_count();
        let mut result = Matrix::zero_ring(upper_part + lower_part, lower_part, ring).into_owned();

        // set to identity in the lower part
        for i in 0..lower_part {
            *result.at_mut(upper_part + i, i) = ring.neg(ring.one());
        }
        
        // set the interesting upper part
        let mut result_upper_part = result.submatrix_mut(..upper_part, ..);
        result_upper_part.add_assign(effective_matrix.submatrix(.., upper_part..), ring);
        effective_matrix.submatrix(.., ..upper_part).solve_strict_triangular(
            &mut result_upper_part, ring
        );

        // and now perform the swaps
        for (row1, row2) in col_swaps.iter().rev() {
            result.swap_rows(*row1, *row2);
        }

        return Some(result);
    }
}


#[test]
fn test_invert_matrix() {
    let mut a = Matrix::from_array([[1., 2.], [4., 8.]]);
    let mut a_inv = Matrix::identity(2, 2);
    assert!(a.solve_modifying(&mut a_inv, &f32::RING).is_err());

    let mut b = Matrix::from_array([[1., 2.], [4., 4.]]);
    let mut b_inv = Matrix::identity(2, 2);
    b.solve_modifying(&mut b_inv, &f32::RING).unwrap();

    assert_eq!(1., *b_inv.at(1, 0));
    assert_eq!(-0.25, *b_inv.at(1, 1));
    assert_eq!(-1., *b_inv.at(0, 0));
    assert_eq!(0.5, *b_inv.at(0, 1));
}

#[test]
fn test_kernel_base() {
    #[rustfmt::skip]
    let mut a = Matrix::from_array([[1., 3., 1., 3.], 
                                    [2., 6., 1., 2.], 
                                    [0., 0., 1., 1.]]);

    let b = Matrix::from_array([[3.], [-1.], [0.], [0.]]);
    assert_eq!(b, a.kernel_base_modifying(&f32::RING).unwrap());

    #[rustfmt::skip]
    let mut a = Matrix::from_array([[1., 3., 1., 3.], 
                                    [2., 6., 1., 5.], 
                                    [0., 0., 1., 1.]]);
                                    
    #[rustfmt::skip]
    let b = Matrix::from_array([[3.,  2.], 
                                [-1., 0.], 
                                [0.,  1.], 
                                [0., -1.]]);
    assert_eq!(b, a.kernel_base_modifying(&f32::RING).unwrap());
}
