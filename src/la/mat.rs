use super::super::alg::*;

pub use super::matrix_view::*;
pub use super::vector_view::*;
pub use super::vec::*;
pub use super::vector::*;
pub use super::submatrix::*;
pub use super::matrix_owned::*;
pub use super::matrix_vector::*;
pub use super::matrix_row_col::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, Add, Mul, RangeBounds, Bound};

#[derive(Debug, Clone, Copy)]
pub struct Matrix<M, T>
    where M: MatrixView<T>
{
    data: M,
    element: PhantomData<T>
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>
{
    pub fn new(data: M) -> Self {
        Matrix {
            data: data,
            element: PhantomData
        }
    }

    pub fn row_count(&self) -> usize {
        self.data.row_count()
    }

    pub fn col_count(&self) -> usize {
        self.data.col_count()
    }

    pub fn at(&self, row: usize, col: usize) -> &T {
        self.data.at(row, col)
    }

    pub fn as_ref<'a>(&'a self) -> Matrix<MatrixRef<'a, M, T>, T> {
        self.submatrix(.., ..)
    }

    pub fn submatrix<'a, R, S>(
        &'a self,
        rows: R, 
        cols: S
    ) -> Matrix<MatrixRef<'a, M, T>, T> 
        where R: RangeBounds<usize>, S: RangeBounds<usize>
    {
        let rows_begin = match rows.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let rows_end = match rows.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.row_count(),
        };
        let cols_begin = match cols.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let cols_end = match cols.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.col_count(),
        };
        Matrix::new(MatrixRef::new(
            rows_begin, rows_end, cols_begin, cols_end, &self.data
        ))
    }

    pub fn rows(&self) -> MatrixRowIter<T, M> {
        MatrixRowIter::new(&self.data)
    }

    pub fn row(&self, row: usize) -> Vector<MatrixRow<T, M>, T> {
        self.data.assert_row_in_range(row);
        Vector::new(MatrixRow::new(&self.data, row))
    }

    pub fn cols(&self) -> MatrixColIter<T, M> {
        MatrixColIter::new(&self.data)
    }

    pub fn col(&self, col: usize) -> Vector<MatrixCol<T, M>, T> {
        self.data.assert_col_in_range(col);
        Vector::new(MatrixCol::new(&self.data, col))
    }
}

impl<V, T> Matrix<ColumnVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn col_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(ColumnVector::new(vector))
    }
}

impl<V, T> Matrix<RowVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn row_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(RowVector::new(vector))
    }
}

impl<T> Matrix<MatrixOwned<T>, T> {

    pub fn from_fn<F>(rows: usize, cols: usize, f: F) -> Self 
        where F: FnMut(usize, usize) -> T
    {
        Self::new(MatrixOwned::from_fn(rows, cols, f))
    }

    pub fn from_array<const R: usize, const C: usize>(array: [[T; C]; R]) -> Self 
    {
        Self::new(MatrixOwned::from_array(array))
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>
{
    pub fn swap_rows(&mut self, fst: usize, snd: usize) {
        self.data.assert_row_in_range(fst);
        self.data.assert_row_in_range(snd);
        if fst == snd {
            return;
        }
        for col in 0..self.col_count() {
            self.data.swap((fst, col), (snd, col));
        }
    }

    pub fn swap_cols(&mut self, fst: usize, snd: usize) {
        self.data.assert_col_in_range(fst);
        self.data.assert_col_in_range(snd);
        if fst == snd {
            return;
        }
        for row in 0..self.row_count() {
            self.data.swap((row, fst), (row, snd));
        }
    }

    pub fn as_mut<'a>(&'a mut self) -> Matrix<MatrixRefMut<'a, M, T>, T> {
        self.submatrix_mut(.., ..)
    }

    pub fn submatrix_mut<'a, R, S>(
        &'a mut self, 
        rows: R, 
        cols: S
    ) -> Matrix<MatrixRefMut<'a, M, T>, T> 
        where R: RangeBounds<usize>, S: RangeBounds<usize>
    {
        let rows_begin = match rows.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let rows_end = match rows.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.row_count(),
        };
        let cols_begin = match cols.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => x + 1,
            Bound::Unbounded => 0,
        };
        let cols_end = match cols.end_bound() {
            Bound::Included(x) => x + 1,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => self.col_count(),
        };
        Matrix::new(MatrixRefMut::new(
            rows_begin, rows_end, cols_begin, cols_end, &mut self.data
        ))
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixMutRowIter<T>
{
    pub fn row_mut<'a>(
        &'a mut self, 
        row: usize
    ) -> Vector<<M as LifetimeMatrixMutRowIter<'a, T>>::RowRef, T> 
    {
        Vector::new(self.data.get_row_mut(row))
    }

    pub fn rows_mut<'a>(
        &'a mut self
    ) -> impl Iterator<Item = Vector<<M as LifetimeMatrixMutRowIter<'a, T>>::RowRef, T>> 
    {
        self.data.rows_mut().map(|r| Vector::new(r))
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: Clone
{
    pub fn to_owned(&self) -> Matrix<MatrixOwned<T>, T> {
        Matrix::new(
            MatrixOwned::from_fn(self.row_count(), self.col_count(), |i, j| self.at(i, j).clone())
        )
    }

    pub fn copy_vec(&self) -> Vector<VectorOwned<T>, T> {
        assert!(self.col_count() == 1 || self.row_count() == 1);
        if self.col_count() == 1 {
            self.col(0).to_owned()
        } else {
            self.row(0).to_owned()
        }
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixViewMut<T>, T: Clone + std::fmt::Debug
{
    pub fn scale_ring<R>(&mut self, rhs: T, ring: &R) 
        where R: Ring<El = T>
    {
        for i in 0..self.row_count() {
            for j in 0..self.col_count() {
                *self.at_mut(i, j) = ring.mul_ref(self.at(i, j), &rhs);
            }
        }
    }

    ///
    /// Let T be the identity matrix (mxm where this matrix is mxn), in which the entries
    /// [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the values in transform.
    /// This function performs the multiplication A' := T * A, where A is this matrix
    /// 
    /// Complexity O(n)
    ///
    pub fn transform_two_dims_left<R>(&mut self, fst: usize, snd: usize, transform: &[T; 4], ring: &R) 
        where R: Ring<El = T>
    {
        assert!(fst < snd);
        self.data.assert_row_in_range(fst);
        self.data.assert_row_in_range(snd);
        for col in 0..self.col_count() {
            let b = self.at(fst, col).clone();
            *self.at_mut(fst, col) = ring.add(ring.mul_ref(self.at(fst, col), &transform[0]),
                ring.mul_ref(self.at(snd, col), &transform[1]));
            *self.at_mut(snd, col) = ring.add(ring.mul_ref(&b, &transform[2]), 
                ring.mul_ref(self.at(snd, col), &transform[3]));
        }
    }

    ///
    /// Let T be the identity matrix (nxn where this matrix is mxn), in which
    /// the entries [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the
    /// values in transform.
    /// This function performs the multiplication A' := A * T, where A is this matrix
    /// 
    /// Complexity O(m)
    ///
    pub fn transform_two_dims_right<R>(&mut self, fst: usize, snd: usize, transform: &[T; 4], ring: &R) 
        where R: Ring<El = T>
    {
        assert!(fst < snd);
        self.data.assert_col_in_range(fst);
        self.data.assert_col_in_range(snd);
        for row in 0..self.row_count() {
            let b = self.at(row, fst).clone();
            *self.at_mut(row, fst) = ring.add(ring.mul_ref(self.at(row, fst), &transform[0]),
                ring.mul_ref(self.at(row, snd), &transform[2]));
            *self.at_mut(row, snd) = ring.add(ring.mul_ref(&b, &transform[1]),
                ring.mul_ref(self.at(row, snd), &transform[3]));
        }
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixView<T>, T: Clone + std::fmt::Debug
{
    pub fn frobenius_square<R>(&self, ring: &R) -> T 
        where R: Ring<El = T>
    {
        let mut it = self.rows().flat_map(|r| 
            (0..r.len()).map(move |i| ring.mul_ref(r.at(i), r.at(i)))
        );
        let initial = it.next().unwrap();
        it.fold(initial, |a, b| ring.add(a, b))
    }

    ///
    /// Expects self to be a square, strict upper triangular matrix (i.e. 1 on the diagonal)
    /// and assigns to rhs the solution of the equation self * X = rhs
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn solve_strict_triangular<N, R>(&self, rhs: &mut Matrix<N, T>, ring: &R)
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
    where M: MatrixViewMut<T>
{
    pub fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.data.at_mut(row, col)
    }
}

impl<M, N, T, U> PartialEq<Matrix<N, U>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<U>, T: PartialEq<U>
{
    fn eq(&self, rhs: &Matrix<N, U>) -> bool {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                if self.at(row, col) != rhs.at(row, col) {
                    return false;
                }
            }
        }
        return true;
    }
}

impl<M, T> Eq for Matrix<M, T>
    where M: MatrixView<T>, T: Eq
{}

impl<M, N, T> AddAssign<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<T>, T: RingEl
{
    fn add_assign(&mut self, rhs: Matrix<N, T>) {
        self.add_assign_ring(rhs, &StaticRing::<T>::RING)
    }
}

impl<M, N, T> Add<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<T>, T: RingEl
{
    type Output = Matrix<MatrixOwned<T>, T>;

    fn add(self, rhs: Matrix<N, T>) -> Self::Output {
        self.add_ring(rhs, &StaticRing::<T>::RING)
    }
}

impl<M, N, T> Mul<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<T>, T: RingEl
{
    type Output = Matrix<MatrixOwned<T>, T>;

    fn mul(self, rhs: Matrix<N, T>) -> Self::Output {
        self.mul_ring(rhs, &StaticRing::<T>::RING)
    }
}

impl<M, T> Matrix<M, T>
    where T: RingEl, M: MatrixViewMut<T>
{
    pub fn scale(&mut self, rhs: T) {
        self.scale_ring(rhs, &StaticRing::<T>::RING)
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, 
        T: std::fmt::Debug + Clone
{
    pub fn mul_ring<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Matrix<MatrixOwned<T>, T> 
        where R: Ring<El = T>, N: MatrixView<T>
    {
        assert_eq!(self.col_count(), rhs.row_count());
        debug_assert!(self.col_count() > 0);
        Matrix::new(
            MatrixOwned::from_fn(self.row_count(), rhs.col_count(), |i, j| {
                let mut it = (0..self.col_count()).map(|k| 
                    ring.mul_ref(self.at(i, k), rhs.at(k, j))
                );
                let initial = it.next().unwrap();
                it.fold(initial, |a, b| ring.add(a, b))
            })
        )
    }

    pub fn add_ring<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Matrix<MatrixOwned<T>, T> 
        where R: Ring<El = T>, N: MatrixView<T>
    {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        Matrix::new(
            MatrixOwned::from_fn(self.row_count(), self.col_count(), |i, j| {
                ring.add_ref(self.at(i, j).clone(), rhs.at(i, j))
            })
        )
    }

    pub fn eq_ring<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> bool
        where R: Ring<El = T>, N: MatrixView<T>
    {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                if !ring.eq(self.at(row, col), rhs.at(row, col)) {
                    return false;
                }
            }
        }
        return true;
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>, T: std::fmt::Debug + Clone
{
    fn add_assign_ring<R, N>(&mut self, rhs: Matrix<N, T>, ring: &R) 
        where R: Ring<El = T>, N: MatrixView<T>
    {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                take_mut::take_or_recover(
                    self.at_mut(row, col), 
                    || ring.unspecified_element(), 
                    |v| ring.add_ref(v, rhs.at(row, col))
                );
            }
        }
    }
}

impl<T> Matrix<MatrixOwned<T>, T>
    where T: Zero
{
    pub fn zero(rows: usize, cols: usize) -> Self {
        Matrix::new(
            MatrixOwned::from_fn(rows, cols, |_, _| T::zero())
        )
    }
}

impl<T> Matrix<MatrixOwned<T>, T>
    where T: std::fmt::Debug + Clone
{
    pub fn zero_ring<R>(rows: usize, cols: usize, ring: &R) -> Self 
        where R: Ring<El = T>
    {
        Matrix::new(
            MatrixOwned::from_fn(rows, cols, |_, _| ring.zero())
        )
    }

    pub fn identity_ring<R>(rows: usize, cols: usize, ring: &R) -> Self 
        where R: Ring<El = T>
    {
        Matrix::new(
            MatrixOwned::from_fn(rows, cols, |i, j| 
                if i == j { ring.one() } else { ring.zero() }
            )
        )
    }
}

impl<T> Matrix<MatrixOwned<T>, T>
    where T: Zero + One
{
    pub fn identity(rows: usize, cols: usize) -> Self {
        Matrix::new(
            MatrixOwned::from_fn(rows, cols, |i, j| 
                if i == j { T::one() } else { T::zero() }
            )
        )
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
    /// Complexity O(n^3)
    /// 
    fn gaussion_elimination_half<R, F, G, H, S>(
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
        assert!(ring.is_field());

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
            self.submatrix_mut(i..=i, i..).scale_ring(inverse.clone(), ring);
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
    pub fn solve_modifying<R, N>(&mut self, rhs: &mut Matrix<N, T>, ring: &R) -> Result<(), usize>
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        assert!(ring.is_field());
        assert_eq!(self.row_count(), self.col_count());
        assert_eq!(self.row_count(), rhs.row_count());
        
        self.gaussion_elimination_half(
            |row, a, rhs| rhs.submatrix_mut(row..=row, 0..).scale_ring(a, ring), 
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
    pub fn kernel_base_modifying<R>(&mut self, ring: &R) -> Option<Matrix<MatrixOwned<T>, T>> 
        where R: Ring<El = T>
    {
        assert!(ring.is_field());

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
        let mut result = Matrix::zero_ring(upper_part + lower_part, lower_part, ring);

        // set to identity in the lower part
        for i in 0..lower_part {
            *result.at_mut(upper_part + i, i) = ring.neg(ring.one());
        }
        
        // set the interesting upper part
        let mut result_upper_part = result.submatrix_mut(..upper_part, ..);
        result_upper_part.add_assign_ring(effective_matrix.submatrix(.., upper_part..), ring);
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

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: FieldEl + Clone
{
    ///
    /// Solves the linear equation AX = B for a square matrix A.
    /// If A is singular, returns the smallest integer i such that
    /// the first i columns of A are linearly dependent. In this
    /// case, the rhs matrix will be in an unspecified state when
    /// the function terminates.
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn solve<N>(&self, rhs: &mut Matrix<N, T>) -> Result<(), usize>
        where N: MatrixViewMut<T>
    {
        self.to_owned().solve_modifying(rhs, &StaticRing::<T>::RING)
    }

    ///
    /// Inverts a square matrix. If the matrix is singular, instead
    /// returns as error the smallest integer i such that the first
    /// i columns of A are linearly dependent.
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^3)
    /// 
    pub fn invert(&self) -> Result<Matrix<MatrixOwned<T>, T>, usize>
    {
        let mut result = Matrix::identity(self.row_count(), self.col_count());
        self.solve(&mut result)?;
        return Ok(result);
    }

    pub fn kernel_base(&self) -> Option<Matrix<MatrixOwned<T>, T>> {
        self.to_owned().kernel_base_modifying(&StaticRing::<T>::RING)
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: std::fmt::Debug + Clone
{
    ///
    /// Solves the linear equation AX = B for a square matrix A.
    /// If A is singular, returns the smallest integer i such that
    /// the first i columns of A are linearly dependent. In this
    /// case, the rhs matrix will be in an unspecified state when
    /// the function terminates.
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^2(n + m)) where self is nxn and rhs is nxm
    /// 
    pub fn solve_ring<R, N>(&self, rhs: &mut Matrix<N, T>, ring: &R) -> Result<(), usize>
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        assert!(ring.is_field());
        self.to_owned().solve_modifying(rhs,ring)
    }

    ///
    /// Inverts a square matrix. If the matrix is singular, instead
    /// returns as error the smallest integer i such that the first
    /// i columns of A are linearly dependent.
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^3)
    /// 
    pub fn invert_ring<R>(&self, ring: &R) -> Result<Matrix<MatrixOwned<T>, T>, usize>
        where R: Ring<El = T>
    {
        assert!(ring.is_field());
        let mut result = Matrix::identity_ring(self.row_count(), self.col_count(), ring);
        self.solve_ring(&mut result, ring)?;
        return Ok(result);
    }

    pub fn kernel_base_ring<R>(&self, ring: &R) -> Option<Matrix<MatrixOwned<T>, T>> 
        where R: Ring<El = T>
    {
        assert!(ring.is_field());
        self.to_owned().kernel_base_modifying(ring)
    }
}

impl<M, T> std::fmt::Display for Matrix<M, T> 
    where M: MatrixView<T>, T: std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let entries = (0..self.row_count()).flat_map(|row| 
            (0..self.col_count()
        ).map(move |col| 
            format!("{}", self.at(row, col))
        )).collect::<Vec<_>>();

        let width = entries.iter().map(|s| s.chars().count()).max().unwrap();
        writeln!(f, "[")?;
        for row in 0..self.row_count() {
            write!(f, "[")?;
            write!(f, "{:>width$}", entries[row * self.col_count()], width = width)?;
            for col in 1..self.col_count() {
                write!(f, ", {:>width$}", entries[row * self.col_count() + col], width = width)?;
            }
            if row + 1 == self.row_count() {
                writeln!(f, "]")?;
            } else {
                writeln!(f, "],")?;
            }
        }
        write!(f, "]")?;
        return Ok(());
    }
}

#[cfg(test)]
use super::super::algebra::zn::*;

#[test]
fn test_mul_matrix() {
    let a: Matrix<_, i32> = Matrix::from_array([[0, 1], [1, 2]]);
    let res = a.as_ref() * a.as_ref();

    assert_eq!(1, *res.at(0, 0));
    assert_eq!(2, *res.at(0, 1));
    assert_eq!(2, *res.at(1, 0));
    assert_eq!(5, *res.at(1, 1));
}

#[test]
fn test_invert_matrix() {
    let mut a = Matrix::from_array([[1., 2.], [4., 8.]]);
    let mut a_inv = Matrix::identity(2, 2);
    assert!(a.solve_modifying(&mut a_inv, &StaticRing::<f32>::RING).is_err());

    let mut b = Matrix::from_array([[1., 2.], [4., 4.]]);
    let mut b_inv = Matrix::identity(2, 2);
    b.solve_modifying(&mut b_inv, &StaticRing::<f32>::RING).unwrap();

    assert_eq!(1., *b_inv.at(1, 0));
    assert_eq!(-0.25, *b_inv.at(1, 1));
    assert_eq!(-1., *b_inv.at(0, 0));
    assert_eq!(0.5, *b_inv.at(0, 1));
}

#[test]
fn test_row_iter() {
    let a = Matrix::from_fn(4, 4, |i, j| i + 4 * j);
    let b = a.submatrix(1..3, 2..4);
    let mut it = b.rows();
    let r1 = it.next().unwrap();
    let r2 = it.next().unwrap();
    assert!(it.next().is_none());

    assert_eq!(2, r1.len());
    assert_eq!(2, r2.len());
    assert_eq!(9, *r1.at(0));
    assert_eq!(13, *r1.at(1));
    assert_eq!(10, *r2.at(0));
    assert_eq!(14, *r2.at(1));

    assert_eq!(3, *a.row(3).at(0));
}

#[test]
fn test_mul() {
    let a = Matrix::from_array([[1, 2], [0, 1]]);
    let b = Matrix::from_array([[0, 2, 1], [2, 0, 1]]);
    let c = Matrix::from_array([[4, 2, 3], [2, 0, 1]]);

    assert_eq!(c, a * b);
}

#[test]
fn test_kernel_base() {
    #[rustfmt::skip]
    let mut a = Matrix::from_array([[1., 3., 1., 3.], 
                                    [2., 6., 1., 2.], 
                                    [0., 0., 1., 1.]]);

    let b = Matrix::from_array([[3.], [-1.], [0.], [0.]]);
    assert_eq!(b, a.kernel_base_modifying(&StaticRing::<f32>::RING).unwrap());

    #[rustfmt::skip]
    let mut a = Matrix::from_array([[1., 3., 1., 3.], 
                                    [2., 6., 1., 5.], 
                                    [0., 0., 1., 1.]]);
                                    
    #[rustfmt::skip]
    let b = Matrix::from_array([[3.,  2.], 
                                [-1., 0.], 
                                [0.,  1.], 
                                [0., -1.]]);
    assert_eq!(b, a.kernel_base_modifying(&StaticRing::<f32>::RING).unwrap());
}
