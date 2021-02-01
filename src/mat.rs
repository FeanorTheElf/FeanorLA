use super::matrix_view::*;
use super::vector_view::*;
use super::vec::*;
use super::vector::*;
use super::submatrix::*;
use super::matrix_owned::*;
use super::matrix_vector::*;
use super::matrix_row_col::*;
use super::alg::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, Add, Mul, RangeBounds, Bound, MulAssign};

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

    pub fn submatrix<'a, R, S>(&'a self, rows: R, cols: S) -> Matrix<MatrixRef<'a, M, T>, T> 
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
        Matrix::new(MatrixRef::new(rows_begin, rows_end, cols_begin, cols_end, &self.data))
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

    pub fn submatrix_mut<'a, R, S>(&'a mut self, rows: R, cols: S) -> Matrix<MatrixRefMut<'a, M, T>, T> 
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
        Matrix::new(MatrixRefMut::new(rows_begin, rows_end, cols_begin, cols_end, &mut self.data))
    }

    pub fn scale<U>(&mut self, rhs: U) 
        where T: MulAssign<U>, U: Clone
    {
        for i in 0..self.row_count() {
            for j in 0..self.col_count() {
                *self.at_mut(i, j) *= rhs.clone();
            }
        }
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixMutRowIter<T>
{
    pub fn row_mut<'a>(&'a mut self, row: usize) -> Vector<<M as LifetimeMatrixMutRowIter<'a, T>>::RowRef, T> {
        Vector::new(self.data.get_row_mut(row))
    }

    pub fn rows_mut<'a>(&'a mut self) -> impl Iterator<Item = Vector<<M as LifetimeMatrixMutRowIter<'a, T>>::RowRef, T>> {
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
    where M: MatrixViewMut<T>, T: RingEl + Clone
{
    ///
    /// Let T be the identity matrix (mxm where this matrix is mxn), in which the entries
    /// [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the values in transform.
    /// This function performs the multiplication A' := T * A, where A is this matrix
    /// 
    /// Complexity O(n)
    ///
    pub fn transform_two_dims_left(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        assert!(fst < snd);
        self.data.assert_row_in_range(fst);
        self.data.assert_row_in_range(snd);
        for col in 0..self.col_count() {
            let b = self.at(fst, col).clone();
            *self.at_mut(fst, col) = self.at(fst, col).clone() * transform[0].clone() + self.at(snd, col).clone() * transform[1].clone();
            *self.at_mut(snd, col) = b * transform[2].clone() + self.at(snd, col).clone() * transform[3].clone();
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
    pub fn transform_two_dims_right(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        assert!(fst < snd);
        self.data.assert_col_in_range(fst);
        self.data.assert_col_in_range(snd);
        for row in 0..self.row_count() {
            let b = self.at(row, fst).clone();
            *self.at_mut(row, fst) = self.at(row, fst).clone() * transform[0].clone() + self.at(row, snd).clone() * transform[2].clone();
            *self.at_mut(row, snd) = b * transform[1].clone() + self.at(row, snd).clone() * transform[3].clone();
        }
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixView<T>, T: RingEl + Clone
{
    pub fn frobenius_square(&self) -> T {
        let mut it = self.rows().flat_map(|r| (0..r.len()).map(move |i| r.at(i).clone() * r.at(i).clone()));
        let initial = it.next().unwrap();
        it.fold(initial, |a, b| a + b)
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

impl<M, N, T, U> AddAssign<Matrix<N, U>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<U>, T: AddAssign<U>, U: Clone
{
    fn add_assign(&mut self, rhs: Matrix<N, U>) {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                *self.at_mut(row, col) += rhs.at(row, col).clone();
            }
        }
    }
}

impl<M, N, T, U> Add<Matrix<N, U>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<U>, T: Add<U> + Clone, U: Clone
{
    type Output = Matrix<MatrixOwned<T::Output>, T::Output>;

    fn add(self, rhs: Matrix<N, U>) -> Self::Output {
        assert_eq!(self.row_count(), rhs.row_count());
        assert_eq!(self.col_count(), rhs.col_count());
        Matrix::new(
            MatrixOwned::from_fn(self.row_count(), self.col_count(), |i, j| {
                self.at(i, j).clone() + rhs.at(i, j).clone()
            })
        )
    }
}

impl<M, N, T, U, R> Mul<Matrix<N, U>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<U>, T: Mul<U, Output = R> + Clone, U: Clone, R: Add<R, Output = R>
{
    type Output = Matrix<MatrixOwned<R>, R>;

    fn mul(self, rhs: Matrix<N, U>) -> Self::Output {
        assert_eq!(self.col_count(), rhs.row_count());
        debug_assert!(self.col_count() > 0);
        Matrix::new(
            MatrixOwned::from_fn(self.row_count(), rhs.col_count(), |i, j| {
                let mut it = (0..self.col_count()).map(|k| self.at(i, k).clone() * rhs.at(k, j).clone());
                let initial = it.next().unwrap();
                it.fold(initial, |a, b| a + b)
            })
        )
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
    where T: Zero + One
{
    pub fn identity(rows: usize, cols: usize) -> Self {
        Matrix::new(
            MatrixOwned::from_fn(rows, cols, |i, j| if i == j { T::one() } else { T::zero() })
        )
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>, T: FieldEl + Clone
{
    ///
    /// Transforms this matrix into strict upper triangle form by gaussian elimination.
    /// If the matrix is singular, the process is aborted as soon as this is detected,
    /// and an index i is returned such that the first i columns of this matrix are in
    /// upper triangle form, with the (i,i)-th entry being zero.
    /// 
    /// The passed functions are called whenever a gaussian elimination step is made,
    /// so that the caller can adjust other data as well (i.e. a potential right-hand side
    /// of a linear equation, or a determinant)
    /// 
    /// Use not for types that have rounding errors, as the algorithm
    /// can be numerically unstable
    /// 
    /// Complexity O(n^3)
    /// 
    fn gaussion_elimination_half<F, G, H, S>(&mut self, mut mul_row: F, mut swap_rows: G, mut sub_row: H, state: &mut S) -> Result<(), usize>
        where F: FnMut(usize, T, &mut S), G: FnMut(usize, usize, &mut S), H: FnMut(usize, T, usize, &mut S)
    {
        for i in 0..self.col_count() {
            // pivot
            if *self.at(i, i) == T::zero() {
                let mut has_swapped = false;
                for j in (i + 1)..self.row_count() {
                    if *self.at(j, i) != T::zero() {
                        has_swapped = true;
                        self.swap_rows(i, j);
                        swap_rows(i, j, state);
                    }
                }
                if !has_swapped {
                    return Err(i);
                }
            }

            // normalize pivot to 1
            let inverse = T::one() / self.at(i, i).clone();
            self.submatrix_mut(i..=i, i..).scale(inverse.clone());
            mul_row(i, inverse, state);

            // eliminate
            for j in (i + 1)..self.row_count() {
                let transform = [T::one(), T::zero(), -self.at(j, i).clone(), T::one()];
                sub_row(j, self.at(j, i).clone(), i, state);
                self.transform_two_dims_left(i, j, &transform);
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
    pub fn solve_modifying<N>(&mut self, rhs: &mut Matrix<N, T>) -> Result<(), usize>
        where N: MatrixViewMut<T>
    {
        assert_eq!(self.row_count(), self.col_count());
        assert_eq!(self.row_count(), rhs.row_count());
        
        self.gaussion_elimination_half(
            |row, a, rhs| rhs.submatrix_mut(row..=row, 0..).scale(a), 
            |i, j, rhs| rhs.swap_rows(i, j), 
            |dst, a, src, rhs| rhs.transform_two_dims_left(src, dst, &[T::one(), T::zero(), -a, T::one()]), rhs)?;

        // now self is in upper triangle form
        for col in 0..rhs.col_count() {
            for row in (0..self.row_count()).rev() {
                for i in (row + 1)..self.row_count() {
                    let d = self.at(row, i).clone() * rhs.at(i, col).clone();
                    *rhs.at_mut(row, col) -= d;
                }
            }
        }
        return Ok(());
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
        self.to_owned().solve_modifying(rhs)
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
}

#[test]
fn test_mul_matrix() {
    let a = Matrix::from_array([[0, 1], [1, 2]]);
    let big = Matrix::from_fn(2, 2, |_, _| a.clone());
    let res = big.as_ref() * big.as_ref();

    assert_eq!(2, *res.at(0, 0).at(0, 0));
    assert_eq!(10, *res.at(0, 1).at(1, 1));
}

#[test]
fn test_invert_matrix() {
    let mut a = Matrix::from_array([[1., 2.], [4., 8.]]);
    let mut a_inv = Matrix::identity(2, 2);
    assert!(a.solve_modifying(&mut a_inv).is_err());

    let mut b = Matrix::from_array([[1., 2.], [4., 4.]]);
    let mut b_inv = Matrix::identity(2, 2);
    b.solve_modifying(&mut b_inv).unwrap();

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