use super::super::alg::*;
use super::ops::*;

pub use super::matrix_view::*;
pub use super::vector_view::*;
pub use super::vec::*;
pub use super::submatrix::*;
pub use super::vector::*;

use super::matrix_transpose::*;
use super::diagonal::*;
use super::matrix_vector::*;
use super::matrix_row_col::*;
use super::constant::*;

use std::marker::PhantomData;
use std::ops::{AddAssign, Sub, SubAssign, MulAssign, Add, Mul, RangeBounds, Bound};

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

    pub fn as_ref<'a>(&'a self) -> Matrix<MatrixRef<'a, M, T>, T> {
        self.submatrix(.., ..)
    }

    pub fn transpose(self) -> Matrix<MatrixTranspose<T, M>, T> {
        Matrix::new(MatrixTranspose::new(self.data))
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

    pub fn rows(&self) -> MatrixRowIter<T, MatrixRef<M, T>> {
        MatrixRowIter::new(self.as_ref().data)
    }

    pub fn row(&self, row: usize) -> Vector<MatrixRow<T, MatrixRef<M, T>>, T> {
        self.data.assert_row_in_range(row);
        Vector::new(MatrixRow::new(self.as_ref().data, row))
    }

    pub fn cols(&self) -> MatrixColIter<T, MatrixRef<M, T>> {
        MatrixColIter::new(self.as_ref().data)
    }

    pub fn col(&self, col: usize) -> Vector<MatrixCol<T, MatrixRef<M, T>>, T> {
        self.data.assert_col_in_range(col);
        Vector::new(MatrixCol::new(self.as_ref().data, col))
    }

    pub fn diag(&self) -> Vector<MatrixDiagonal<MatrixRef<M, T>, T>, T> {
        self.nonmain_diag(0)
    }

    pub fn nonmain_diag(&self, diag_index: i64) -> Vector<MatrixDiagonal<MatrixRef<M, T>, T>, T> {
        Vector::new(MatrixDiagonal::new(self.as_ref().data, diag_index))
    }

    pub fn into_row_vec(self) -> Vector<MatrixRow<T, M>, T> {
        assert!(self.row_count() == 1);
        Vector::new(MatrixRow::new(self.data, 0))
    }

    pub fn into_col_vec(self) -> Vector<MatrixCol<T, M>, T> {
        assert!(self.col_count() == 1);
        Vector::new(MatrixCol::new(self.data, 0))
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixViewMut<T>
{
    pub fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.data.at_mut(row, col)
    }

    pub fn col_mut(&mut self, col: usize) -> Vector<MatrixCol<T, MatrixRefMut<M, T>>, T> {
        self.data.assert_col_in_range(col);
        Vector::new(MatrixCol::new(self.as_mut().data, col))
    }
}

impl<V, T> Matrix<ColumnVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn col_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(vector.as_column_vector())
    }
}

impl<V, T> Matrix<RowVector<V, T>, T>
    where V: VectorView<T>
{
    pub fn row_vec(vector: Vector<V, T>) -> Self {
        Matrix::new(vector.as_row_vector())
    }
}

impl<V, T> Matrix<DiagonalMatrix<V, T>, T>
    where V: VectorView<T>
{
    pub fn diag_matrix_ring<R: Ring<El = T>>(vector: Vector<V, T>, ring: &R) -> Self {
        Matrix::nonmain_diag_matrix_ring(vector, 0, ring)
    }

    pub fn nonmain_diag_matrix_ring<R: Ring<El = T>>(vector: Vector<V, T>, diag_index: i64, ring: &R) -> Self {
        Matrix::new(vector.as_diag_matrix(diag_index, ring.zero()))
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
    pub fn assign<N>(&mut self, rhs: Matrix<N, T>) 
        where N: MatrixView<T>
    {
        <T as MatrixAssign<M, N>>::assign_matrix(&mut self.data, rhs.data)
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>, T: Clone + std::fmt::Debug
{
    pub fn mul<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Matrix<MatrixOwned<T>, T>
        where R: Ring<El = T>, N: MatrixView<T>
    {
        Matrix::new(<R as MatrixMul<M, N>>::mul_matrix(ring, self.data, rhs.data))
    }

    pub fn add<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Matrix<MatrixOwned<T>, T>
        where R: Ring<El = T>, N: MatrixView<T>
    {
        let mut result = self.into_owned();
        result.add_assign(rhs, ring);
        return result;
    }


    pub fn sub<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> Matrix<MatrixOwned<T>, T>
        where R: Ring<El = T>, N: MatrixView<T>
    {
        let mut result = self.into_owned();
        result.sub_assign(rhs, ring);
        return result;
    }

    pub fn eq<R, N>(self, rhs: Matrix<N, T>, ring: &R) -> bool
        where R: Ring<El = T>, N: MatrixView<T>
    {
        <R as MatrixEq<M, N>>::eq_matrix(ring, self.data, rhs.data)
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

    pub fn add_assign<R, N>(&mut self, rhs: Matrix<N, T>, ring: &R)
        where R: Ring<El = T>, N: MatrixView<T>
    {
        <R as MatrixAddAssign<M, N>>::add_assign_matrix(ring, &mut self.data, rhs.data);
    }

    pub fn sub_assign<R, N>(&mut self, rhs: Matrix<N, T>, ring: &R)
        where R: Ring<El = T>, N: MatrixView<T>
    {
        <R as MatrixAddAssign<M, N>>::sub_assign_matrix(ring, &mut self.data, rhs.data);
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

impl<M, N, T> AddAssign<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<T>, T: RingEl
{
    fn add_assign(&mut self, rhs: Matrix<N, T>) {
        self.add_assign(rhs, &T::RING)
    }
}

impl<M, N, T> SubAssign<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<T>, T: RingEl
{
    fn sub_assign(&mut self, rhs: Matrix<N, T>) {
        self.sub_assign(rhs, &T::RING)
    }
}

impl<M, N, T> Add<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixViewMut<T>, N: MatrixView<T>, T: RingEl
{
    type Output = Matrix<MatrixOwned<T>, T>;

    fn add(self, rhs: Matrix<N, T>) -> Self::Output {
        self.add(rhs, &T::RING)
    }
}

impl<M, N, T> Sub<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<T>, T: RingEl
{
    type Output = Matrix<MatrixOwned<T>, T>;

    fn sub(self, rhs: Matrix<N, T>) -> Self::Output {
        self.sub(rhs, &T::RING)
    }
}

impl<M, N, T> Mul<Matrix<N, T>> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<T>, T: RingEl
{
    type Output = Matrix<MatrixOwned<T>, T>;

    fn mul(self, rhs: Matrix<N, T>) -> Self::Output {
        self.mul(rhs, &T::RING)
    }
}

impl<M, T> MulAssign<T> for Matrix<M, T>
    where T: RingEl, M: MatrixViewMut<T>
{
    fn mul_assign(&mut self, rhs: T) {
        self.scale(&rhs, &T::RING)
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
        self.matrix.format_raw(|x| Ok(format!("{}", display_ring_el(&self.ring, x))), f)
    }
}


impl<M, T> Matrix<M, T> 
    where M: MatrixView<T>, T: std::fmt::Debug + Clone
{
    pub fn display<'a, 'b, R>(&'b self, ring: &'a R) -> DisplayMatrix<MatrixRef<'b, M, T>, T, &'a R> 
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
    let res = a.as_ref() * a.as_ref();

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

    assert_eq!(c, a * b);
}

#[test]
fn test_zero_sized_matrix() {
    let mut a: Matrix<_, i32> = Matrix::from_array([[], [], []]);
    assert_eq!(3, a.row_count());
    assert_eq!(0, a.col_count());

    let b = a.as_ref();
    assert_eq!(3, b.row_count());
    assert_eq!(0, b.col_count());

    let c = b.submatrix(..2, ..);
    assert_eq!(2, c.row_count());
    assert_eq!(0, c.col_count());

    assert_eq!(Matrix::<_, i32>::zero(3, 0), a);

    let mut d: Matrix<_, i32> = Matrix::zero(0, 5).into_owned();

    assert_eq!(d.clone(), d);

    assert_eq!(1, d.submatrix(.., 2..3).col_count());
    assert_eq!(a.submatrix(..0, ..), d.submatrix(.., 5..));

    assert_eq!(3, a.rows().count());
    assert_eq!(3, a.rows_mut().count());
    assert_eq!(0, d.rows().count());
    assert_eq!(0, d.rows_mut().count());

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