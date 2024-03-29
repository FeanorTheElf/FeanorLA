use super::super::ring::*;
use super::super::primitive::*;
use super::mat_fn::*;
use super::ops::*;

pub use super::mat_fn::MatFn;
pub use super::matrix_view::*;
pub use super::matrix_view::constant_value_matrix::*;
pub use super::matrix_view::matrix_transpose::*;
pub use super::matrix_view::matrix_vector::*;
pub use super::matrix_view::diagonal_matrix::*;
pub use super::matrix_view::submatrix::*;
pub use super::vec::*;
pub use super::vector_view::*;
pub use super::vector_view::matrix_row_col::*;
pub use super::vector_view::matrix_diagonal::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, SubAssign, Mul, RangeBounds, Bound};

#[derive(Debug)]
pub struct Matrix<M, T>
    where M: MatrixView<T>
{
    data: M,
    element: PhantomData<T>
}

impl<M, T> Copy for Matrix<M, T>
    where M: MatrixView<T> + Copy
{}

impl<M, T> Clone for Matrix<M, T>
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        Self::new(self.data.clone())
    }
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

    pub fn as_ref<'a>(&'a self) -> Matrix<&'a M, T> {
        Matrix::new(&self.data)
    }

    pub fn transpose(self) -> Matrix<MatrixTranspose<T, M>, T> {
        Matrix::new(MatrixTranspose::new(self.data))
    }

    pub fn into_submatrix<R, S>(self, rows: R, cols: S) -> Matrix<Submatrix<M, T>, T> 
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
        Matrix::new(Submatrix::new(
            rows_begin, rows_end, cols_begin, cols_end, self.data
        ))
    }

    pub fn submatrix<'a, R, S>(&'a self, rows: R, cols: S) -> Matrix<Submatrix<&'a M, T>, T>
        where R: RangeBounds<usize>, S: RangeBounds<usize>
    {
        self.as_ref().into_submatrix(rows, cols)
    }

    pub fn rows(&self) -> std::iter::Map<MatrixRowIter<T, &M>, ToVector<MatrixRow<T, &M>, T>> {
        MatrixRowIter::new(&self.data).map(ToVector::INSTANCE)
    }

    pub fn into_row(self, row: usize) -> Vector<MatrixRow<T, M>, T> {
        self.data.assert_row_in_range(row);
        Vector::new(MatrixRow::new(self.data, row))
    }

    pub fn row(&self, row: usize) -> Vector<MatrixRow<T, &M>, T> {
        self.as_ref().into_row(row)
    }

    pub fn cols(&self) -> std::iter::Map<MatrixColIter<T, &M>, ToVector<MatrixCol<T, &M>, T>> {
        MatrixColIter::new(&self.data).map(ToVector::INSTANCE)
    }

    pub fn into_col(self, col: usize) -> Vector<MatrixCol<T, M>, T> {
        self.data.assert_col_in_range(col);
        Vector::new(MatrixCol::new(self.data, col))
    }

    pub fn col(&self, col: usize) -> Vector<MatrixCol<T, &M>, T> {
        self.as_ref().into_col(col)
    }

    pub fn diag(&self) -> Vector<MatrixDiagonal<&M, T>, T> {
        self.nonmain_diag(0)
    }

    pub fn nonmain_diag(&self, diag_index: i64) -> Vector<MatrixDiagonal<&M, T>, T> {
        self.as_ref().into_nonmain_diag(diag_index)
    }

    pub fn into_nonmain_diag(self, diag_index: i64) -> Vector<MatrixDiagonal<M, T>, T> {
        Vector::new(MatrixDiagonal::new(self.data, diag_index))
    }

    pub fn into_row_vec(self) -> Vector<MatrixRow<T, M>, T> {
        assert!(self.row_count() == 1);
        Vector::new(MatrixRow::new(self.data, 0))
    }

    pub fn into_col_vec(self) -> Vector<MatrixCol<T, M>, T> {
        assert!(self.col_count() == 1);
        Vector::new(MatrixCol::new(self.data, 0))
    }

    pub fn as_el(&self) -> &T {
        assert!(self.col_count() == 1 && self.row_count() == 1);
        return self.at(0, 0);
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>
{
    pub fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.data.at_mut(row, col)
    }

    pub fn col_mut(&mut self, col: usize) -> Vector<MatrixCol<T, &mut M>, T> {
        self.data.assert_col_in_range(col);
        Vector::new(MatrixCol::new(&mut self.data, col))
    }
}

impl<V, T> Matrix<ColumnVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn col_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(ColumnVector::new(vector.into_raw_data()))
    }
}

impl<V, T> Matrix<RowVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn row_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(RowVector::new(vector.into_raw_data()))
    }
}

impl<V, T> Matrix<DiagonalMatrix<V, T>, T>
    where V: VectorView<T>
{
    pub fn diag_matrix_ring<R: Ring<El = T>>(vector: Vector<V, T>, ring: &R) -> Self {
        Matrix::nonmain_diag_matrix_ring(vector, 0, ring)
    }

    pub fn nonmain_diag_matrix_ring<R: Ring<El = T>>(vector: Vector<V, T>, diag_index: i64, ring: &R) -> Self {
        Matrix::new(DiagonalMatrix::new(vector.into_raw_data(), diag_index, ring.zero()))
    }
}

impl<V, T> Matrix<DiagonalMatrix<V, T>, T>
    where V: VectorView<T>, T: RingEl
{
    pub fn diag_matrix(vector: Vector<V, T>) -> Self {
        Matrix::nonmain_diag_matrix(vector, 0)
    }

    pub fn nonmain_diag_matrix(vector: Vector<V, T>, diag_index: i64) -> Self {
        Matrix::nonmain_diag_matrix_ring(vector, diag_index, &T::RING)
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

    // Sadly, we cannot implement the trait From as it would conflict
    // with the default implementation From<T> for T
    pub fn from<U, M>(matrix: Matrix<M, U>) -> Self
        where M: MatrixView<U>, T: From<U>, U: Clone
    {
        Self::map(matrix, T::from)
    }

    pub fn map<U, M, F>(matrix: Matrix<M, U>, f: F) -> Self
        where M: MatrixView<U>, F: FnMut(U) -> T, U: Clone
    {
        let row_count = matrix.row_count();
        let col_count = matrix.col_count();
        let mut data_it = matrix.into_owned().data.into_data_iter().map(f);
        let result = Self::from_fn(row_count, col_count, |_, _| data_it.next().unwrap());
        debug_assert!(data_it.next().is_none());
        return result;
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

    pub fn as_mut<'a>(&'a mut self) -> Matrix<&'a mut M, T> {
        Matrix::new(&mut self.data)
    }

    pub fn submatrix_mut<'a, R, S>(&'a mut self, rows: R, cols: S) -> Matrix<Submatrix<&'a mut M, T>, T> 
        where R: RangeBounds<usize>, S: RangeBounds<usize>
    {
        Matrix::new(&mut self.data).into_submatrix(rows, cols)
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: Clone
{
    pub fn into_owned(self) -> Matrix<MatrixOwned<T>, T> {
        Matrix::new(self.data.into_owned())
    }

    pub fn copy_vec(&self) -> Vector<VectorOwned<T>, T> {
        assert!(self.col_count() == 1 || self.row_count() == 1);
        if self.col_count() == 1 {
            self.as_ref().col(0).into_owned()
        } else {
            self.as_ref().row(0).into_owned()
        }
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>, T: Clone
{
    pub fn assign<N>(&mut self, rhs: N) 
        where N: MatFn<T>
    {
        rhs.assign_to(self)
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: Clone + std::fmt::Debug
{
    pub fn mul<R, N>(self, rhs: Matrix<N, T>, ring: R) -> MatrixProd<R, M, N>
        where R: Ring<El = T>, N: MatrixView<T>
    {
        assert_eq!(self.col_count(), rhs.row_count());
        MatrixProd::new(self, rhs, ring)
    }

    pub fn kronecker<R, N>(self, rhs: Matrix<N, T>, ring: R) -> MatrixKronecker<R, M, N>
        where R: Ring<El = T>, N: MatrixView<T>
    {
        MatrixKronecker::new(self, rhs, ring)
    }

    pub fn eq<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> bool
        where R: Ring<El = T>, N: MatrixView<T>
    {
        <R as MatrixEq<M, N>>::eq_matrix(ring, self.data, rhs.data)
    }

    pub fn frobenius_norm_square<R>(self, ring: &R) -> R::El
        where R: Ring<El = T>
    {
        <R as MatrixFrobenius<M>>::calc_matrix_frobenius_norm_square(ring, self.data)
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixViewMut<T>, T: Clone + std::fmt::Debug
{
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

    pub fn add_assign<R, N>(&mut self, rhs: N, ring: &R)
        where R: Ring<El = T>, N: MatFn<T>
    {
        rhs.add_to(self, ring);
    }

    pub fn sub_assign<R, N>(&mut self, rhs: N, ring: &R)
        where R: Ring<El = T>, N: MatFn<T>
    {
        MatrixNeg::new(rhs, ring).add_to(self, ring)
    }

    pub fn scale<R>(&mut self, rhs: &T, ring: &R)
        where R: Ring<El = T>
    {
        <R as MatrixScale<M>>::scale_matrix(ring, rhs, &mut self.data);
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

impl<M, N, T> AddAssign<N> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatFn<T>, T: RingEl
{
    fn add_assign(&mut self, rhs: N) {
        self.add_assign(rhs, &T::RING)
    }
}

impl<M, N, T> SubAssign<N> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatFn<T>, T: RingEl
{
    fn sub_assign(&mut self, rhs: N) {
        self.sub_assign(rhs, &T::RING)
    }
}

impl<M, N, T> Mul<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<T>, T: RingEl
{
    type Output = MatrixProd<T::RingType, M, N>;

    fn mul(self, rhs: Matrix<N, T>) -> Self::Output {
        self.mul(rhs, T::RING)
    }
}

impl<T> Matrix<MatrixConstant<T>, T>
    where T: Zero
{
    pub fn zero(rows: usize, cols: usize) -> Self {
        Matrix::new(
            MatrixConstant::new(rows, cols, T::zero())
        )
    }
}

impl<T> Matrix<MatrixConstant<T>, T>
    where T: std::fmt::Debug + Clone
{
    pub fn zero_ring<R>(rows: usize, cols: usize, ring: &R) -> Self 
        where R: Ring<El = T>
    {
        Matrix::new(
            MatrixConstant::new(rows, cols, ring.zero())
        )
    }
}

impl<T> Matrix<MatrixOwned<T>, T>
    where T: std::fmt::Debug + Clone
{
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
    where M: MatrixView<T>
{
    fn format_raw<F>(&self, mut format_entry: F, f: &mut std::fmt::Formatter) -> std::fmt::Result 
        where F: FnMut(&T) -> Result<String, std::fmt::Error>
    {
        let entries = (0..self.row_count())
            .flat_map(|row| (0..self.col_count()).map(move |col| (row, col)))
            .map(|(r, c)| format_entry(self.at(r, c)))
            .collect::<Result<Vec<_>, _>>()?;

        let width = entries.iter().map(|s| s.chars().count()).max().unwrap_or(0);
        writeln!(f, "[")?;
        for row in 0..self.row_count() {
            write!(f, "[")?;
            if self.col_count() > 0 {
                write!(f, "{:>width$}", entries[row * self.col_count()], width = width)?;
            }
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

impl<M, T> std::fmt::Display for Matrix<M, T> 
    where M: MatrixView<T>, T: std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.format_raw(|x| Ok(format!("{}", x)), f)
    }
}

pub struct DisplayMatrix<M, T, R>
    where M: MatrixView<T>, R: Ring<El = T>, T: Clone + std::fmt::Debug
{
    matrix: Matrix<M, T>,
    ring: R
}

impl<M, T, R> std::fmt::Display for DisplayMatrix<M, T, R> 
    where M: MatrixView<T>, R: Ring<El = T>, T: Clone + std::fmt::Debug
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.matrix.format_raw(|x| Ok(format!("{}", self.ring.display(x))), f)
    }
}

impl<M, T> Matrix<M, T> 
    where M: MatrixView<T>, T: std::fmt::Debug + Clone
{
    pub fn display<'a, 'b, R>(&'b self, ring: &'a R) -> DisplayMatrix<&'b M, T, &'a R> 
        where R: Ring<El = T>
    {
        DisplayMatrix {
            matrix: self.as_ref(),
            ring: ring
        }
    }
}

#[test]
fn test_mul_matrix() {
    let a: Matrix<_, i32> = Matrix::from_array([[0, 1], [1, 2]]);
    let res = (a.as_ref() * a.as_ref()).compute();

    assert_eq!(1, *res.at(0, 0));
    assert_eq!(2, *res.at(0, 1));
    assert_eq!(2, *res.at(1, 0));
    assert_eq!(5, *res.at(1, 1));
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

    assert_eq!(c, (a * b).compute());
}

#[test]
fn test_zero_sized_matrix() {
    let a: Matrix<_, i32> = Matrix::from_array([[], [], []]);
    assert_eq!(3, a.row_count());
    assert_eq!(0, a.col_count());

    let b = a.as_ref();
    assert_eq!(3, b.row_count());
    assert_eq!(0, b.col_count());

    let c = b.submatrix(..2, ..);
    assert_eq!(2, c.row_count());
    assert_eq!(0, c.col_count());

    assert_eq!(Matrix::<_, i32>::zero(3, 0), a);

    let d: Matrix<_, i32> = Matrix::zero(0, 5).into_owned();

    assert_eq!(d.clone(), d);

    assert_eq!(1, d.submatrix(.., 2..3).col_count());
    assert_eq!(a.submatrix(..0, ..), d.submatrix(.., 5..));

    assert_eq!(3, a.rows().count());
    assert_eq!(0, d.rows().count());

    assert_eq!(0, a.row(0).len());
    assert_eq!(0, d.col(1).len());
}

#[test]
fn test_diagonal() {
    let a = Matrix::from_array([[1, 2], [3, 4], [5, 6]]);
    assert_eq!(2, a.diag().len());
    assert_eq!(1, a.nonmain_diag(1).len());
    assert_eq!(2, a.nonmain_diag(-1).len());
    assert_eq!(1, a.nonmain_diag(-2).len());

    let b = Matrix::from_array([[1, 0], [0, 2]]);
    assert_eq!(b, Matrix::diag_matrix(Vector::from_array([1, 2])));

    let c = Matrix::from_array([[0, 0], [3, 0], [0, 6]]);
    assert_eq!(c, Matrix::nonmain_diag_matrix(a.nonmain_diag(-1), -1));
}