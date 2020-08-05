use super::arith::*;
use super::indexed::{Indexed, IndexedMut};
use std::mem::swap;
use std::ops::{
    Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, RangeBounds, Bound, Range,
    Sub, SubAssign,
};
use super::vector::*;

pub trait MatrixView<T> {
    fn rows(&self) -> usize;
    fn cols(&self) -> usize;
    fn at(&self, row: usize, col: usize) -> &T;

    fn assert_row_in_range(&self, row: usize) {
        assert!(row < self.rows(), "Row index {} out of bounds in matrix with {} rows", row, self.rows());
    }

    fn assert_col_in_range(&self, col: usize) {
        assert!(col < self.cols(), "Column index {} out of bounds in matrix with {} columns", col, self.cols());
    }
}

fn assert_sizes_match<T, M, N>(a: &M, b: &N)
    where M: MatrixView<T>, N: MatrixView<T>
{
    assert!(a.rows() == b.rows() && a.cols() == b.cols(),
        "Expected sizes of added matrices to match, but got {}x{} and {}x{}",
        a.rows(), a.cols(), b.rows(), b.cols()
    );
}

pub trait MatrixViewMut<T>: MatrixView<T> {
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T;
}

///
/// Represents a mxn matrix with elements of type T. Typical matrix operations
/// are not optimized much, so this type is only suitable for small matrices.
/// Instead, the focus lies on a convenient interface and a generic design.
///
#[derive(Debug, Clone)]
pub struct Matrix<T> {
    rows: usize,
    data: Box<[T]>,
}

#[derive(Debug, Clone, Copy)]
pub struct MatrixRef<'a, T> {
    rows_begin: usize,
    rows_end: usize,
    cols_begin: usize,
    cols_end: usize,
    matrix: &'a Matrix<T>,
}

#[derive(Debug)]
pub struct MatrixRefMut<'a, T> {
    rows: Range<usize>,
    cols: Range<usize>,
    matrix: &'a mut Matrix<T>,
}

impl<T> MatrixView<T> for Matrix<T> {
    
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.data.len() / self.rows
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        &self.data[row * self.rows + col]
    }
}

impl<T> MatrixViewMut<T> for Matrix<T> {
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        &mut self.data[row * self.rows + col]
    }
}

impl<'a, T> MatrixView<T> for MatrixRef<'a, T> {
    
    fn rows(&self) -> usize {
        self.rows_end - self.rows_begin
    }

    fn cols(&self) -> usize {
        self.cols_end - self.cols_begin
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.matrix.at(row + self.rows_begin, col + self.cols_begin)
    }
}

impl<'a, T> MatrixView<T> for MatrixRefMut<'a, T> {
    
    fn rows(&self) -> usize {
        self.rows.start
    }

    fn cols(&self) -> usize {
        self.cols.start
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.matrix.at(row + self.rows.start, col + self.cols.start)
    }
}

impl<'a, T> MatrixViewMut<T> for MatrixRefMut<'a, T> {
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.matrix.at_mut(row + self.rows.start, col + self.cols.start)
    }
}

impl<T> Matrix<T> {

    pub fn new(data: Box<[T]>, rows: usize) -> Matrix<T> {
        assert!(data.len() > 0, "Cannot create matrix with zero elements");
        assert!(
            data.len() % rows == 0,
            "Data length must be a multiple of row count, but got {} and {}",
            data.len(),
            rows
        );
        Matrix {
            data: data,
            rows: rows
        }
    }

    pub fn from_array<const R: usize, const C: usize>(mut data: [[T; C]; R]) -> Matrix<T> {
        let data = (0..R).flat_map(|row| (0..C).map(|col| data[row][col])).collect::<Vec<T>>().into_boxed_slice();
        Self::new(data, R)
    }

    pub fn from_nocopy<U>(value: Matrix<U>) -> Self
    where
        T: From<U>,
    {
        let rows = value.rows();
        let data: Vec<T> = value
            .data
            .into_vec()
            .into_iter()
            .map(|d| T::from(d))
            .collect();
        return Matrix {
            data: data.into_boxed_slice(),
            rows: rows,
        };
    }

    pub fn as_ref(&self) -> MatrixRef<T> {
        MatrixRef {
            rows_begin: 0,
            rows_end: self.rows(),
            cols_begin: 0,
            cols_end: self.cols(),
            matrix: self
        }
    }

    pub fn as_mut(&mut self) -> MatrixRefMut<T> {
        MatrixRefMut {
            rows: 0..self.rows(),
            cols: 0..self.cols(),
            matrix: self
        }
    }

    pub fn get_rows<'b>(&'b mut self, fst: usize, snd: usize) -> (VectorRefMut<'b, T>, VectorRefMut<'b, T>) {
        self.assert_row_in_range(fst);
        self.assert_row_in_range(snd);
        assert!(
            fst != snd,
            "When borrowing two rows, their indices must be different, got {}",
            fst
        );

        let cols = self.cols();
        if fst < snd {
            let part: &mut [T] = &mut self.data[(fst * cols)..((snd + 1) * cols)];
            let (fst_row, rest) = part.split_at_mut(cols);
            let snd_row_start = rest.len() - cols;
            return (
                VectorRefMut { data: fst_row },
                VectorRefMut {
                    data: &mut rest[snd_row_start..],
                },
            );
        } else {
            let part: &mut [T] = &mut self.data[(snd * cols)..((fst + 1) * cols)];
            let (snd_row, rest) = part.split_at_mut(cols);
            let fst_row_start = rest.len() - cols;
            return (
                VectorRefMut {
                    data: &mut rest[fst_row_start..],
                },
                VectorRefMut { data: snd_row },
            );
        }
    }
}

impl<T> Matrix<T>
    where T: Clone 
{
    pub fn copy_of<M: MatrixView<T>>(value: M) -> Matrix<T> {
        let data = (0..value.rows()).flat_map(|row| (0..value.cols()).map(|col| value.at(row, col).clone())).collect::<Vec<T>>().into_boxed_slice();
        Self::new(data, value.rows())
    }
}

impl<'a, T> MatrixRefMut<'a, T> {
    pub fn into_const(self) -> MatrixRef<'a, T> {
        MatrixRef {
            rows_begin: self.rows.start,
            rows_end: self.rows.end,
            cols_begin: self.cols.start,
            cols_end: self.cols.end,
            matrix: self.matrix
        }
    }
}

impl<T> Matrix<T>
where
    T: Zero + One + Clone,
{
    pub fn identity(size: usize) -> Matrix<T> {
        let mut result = Matrix::<T>::zero(size, size);
        for i in 0..size {
            *result.at_mut(i, i) = T::one();
        }
        return result;
    }
}

impl<T> Matrix<T>
where
    T: Zero + Clone,
{
    pub fn zero(rows: usize, cols: usize) -> Matrix<T> {
        let mut data: Vec<T> = Vec::new();
        data.resize(rows * cols, T::zero());
        return Matrix::new(data.into_boxed_slice(), rows);
    }
}

impl<'a, T, U: 'a> From<MatrixRef<'a, U>> for Matrix<T>
where
    T: From<&'a U>,
{
    fn from(value: MatrixRef<'a, U>) -> Self {
        let data: Vec<T> = (value.rows_begin..value.rows_end)
            .flat_map(|row| {
                (value.cols_begin..value.cols_end)
                    .map(move |col| &value.matrix.at(row, col))
                    .map(|d| T::from(d))
            })
            .collect();
        return Matrix {
            data: data.into_boxed_slice(),
            rows: value.rows(),
        };
    }
}

impl<'a, T, U: 'a> From<MatrixRefMut<'a, U>> for Matrix<T>
where
    T: From<&'a U>,
{
    fn from(value: MatrixRefMut<'a, U>) -> Self {
        Matrix::from(value.into_const())
    }
}

impl<'a, T> MatrixRef<'a, T>
where
    T: Clone
        + PartialEq<T>
        + AddAssign<T>
        + MulAssign<T>
        + Neg<Output = T>
        + Zero
        + One
        + Div<T, Output = T>
        + Mul<T, Output = T>,
    <T as Div<T>>::Output: Clone,
{
    ///
    /// Calculates the inverse of this matrix. Use only for small matrices, this
    /// is just simple gaussian elimination, and is neither very performant nor
    /// numerically stable!
    ///
    pub fn invert(&self) -> Result<Matrix<T>, ()> {
        assert_eq!(self.rows(), self.cols());
        let n: usize = self.rows();

        let mut result: Matrix<T> = Matrix::identity(n);
        let mut work: Matrix<T> = Matrix::copy_of(*self);

        // just simple gaussian elimination
        for i in 0..n {
            let not_null_index = (i..n).find(|r| work[*r][i] != T::zero()).ok_or(())?;
            if not_null_index != i {
                result.get_mut((.., ..)).swap_rows(i, not_null_index);
                work.get_mut((.., ..)).swap_rows(i, not_null_index);
            }
            for j in (i + 1)..n {
                let (row1, mut row2) = work.get_rows(i, j);
                let factor = -row2[i].clone() / row1[i].clone();
                row2.add_product(row1.as_const(), factor.clone());
                let (res1, mut res2) = result.get_rows(i, j);
                res2.add_product(res1.as_const(), factor);
            }
        }
        // now we have an upper triangle matrix
        for i in 1..n {
            for j in 0..i {
                let (row1, mut row2) = work.get_rows(i, j);
                let factor = -row2[i].clone() / row1[i].clone();
                row2.add_product(row1.as_const(), factor.clone());
                let (res1, mut res2) = result.get_rows(i, j);
                res2.add_product(res1.as_const(), factor);
            }
        }
        // now we have a diagonal matrix
        for i in 0..n {
            result.get_mut(i).mul_assign(T::one() / work[i][i].clone());
        }
        return Ok(result);
    }
}

impl<'a, T> MatrixRefMut<'a, T> {
    
    pub fn get_rows<'b>(&'b mut self, fst: usize, snd: usize) -> (VectorRefMut<'b, T>, VectorRefMut<'b, T>) {
        self.assert_row_in_range(fst);
        self.assert_row_in_range(snd);
        let (mut fst_row, mut snd_row) = self.matrix.get_rows(fst, snd);
        return (fst_row.into_subrange(self.cols.start..self.cols.end), snd_row.into_subrange(self.cols.start..self.cols.end))
    }

    pub fn swap_rows(&mut self, fst: usize, snd: usize) {
        if fst != snd {
            let cols = self.cols();
            let (mut fst_row, mut snd_row) = self.get_rows(fst, snd);
            for col in 0..cols {
                swap(&mut fst_row[col], &mut snd_row[col]);
            }
        }
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        self.assert_row_in_range(index);
        &self.data[(index * self.rows)..(index * self.rows + self.rows)]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.assert_row_in_range(index);
        &mut self.data[(index * self.rows)..(index * self.rows + self.rows)]
    }
}

impl<'a, T> Index<usize> for MatrixRef<'a, T> {
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        self.assert_row_in_range(index);
        &self.matrix[index + self.rows_begin][self.cols_begin..self.cols_end]
    }
}

impl<'a, T> Index<usize> for MatrixRefMut<'a, T> {
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        self.assert_row_in_range(index);
        &self.matrix[index + self.rows.start][self.cols.start..self.cols.end]
    }
}

impl<'a, T> IndexMut<usize> for MatrixRefMut<'a, T> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.assert_row_in_range(index);
        &mut self.matrix[index + self.rows.start][self.cols.start..self.cols.end]
    }
}

fn get_lower_index<R: RangeBounds<usize>>(range: &R, size: usize) -> usize {
    match range.start_bound() {
        Bound::Excluded(i) => *i + 1,
        Bound::Included(i) => *i,
        Bound::Unbounded => 0
    }
}

fn get_upper_index<R: RangeBounds<usize>>(range: &R, size: usize) -> usize {
    match range.end_bound() {
        Bound::Excluded(i) => {
            assert!(*i <= size, "Index {} out of range for size {}", i - 1, size);
            *i
        },
        Bound::Included(i) => {
            assert!(*i < size, "Index {} out of range for size {}", i, size);
            *i + 1
        },
        Bound::Unbounded => size
    }
}

impl<'b, T: 'b, R: RangeBounds<usize>, S: RangeBounds<usize>> Indexed<'b, (R, S)> for Matrix<T> {
    type Output = MatrixRef<'b, T>;

    fn get(&'b self, (rows, cols): (R, S)) -> Self::Output {
        MatrixRef {
            rows_begin: get_lower_index(&rows, self.rows()),
            rows_end: get_upper_index(&rows, self.rows()),
            cols_begin: get_lower_index(&cols, self.cols()),
            cols_end: get_upper_index(&cols, self.cols()),
            matrix: self
        }
    }
}

impl<'b, T: 'b, R: RangeBounds<usize>, S: RangeBounds<usize>> IndexedMut<'b, (R, S)> for Matrix<T> {
    type Output = MatrixRefMut<'b, T>;

    fn get_mut(&'b mut self, (rows, cols): (R, S)) -> Self::Output {
        MatrixRefMut {
            rows: get_lower_index(&rows, self.rows())..get_upper_index(&rows, self.rows()),
            cols: get_lower_index(&cols, self.cols())..get_upper_index(&cols, self.cols()),
            matrix: self
        }
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>, S: RangeBounds<usize>> Indexed<'b, (R, S)> for MatrixRef<'a, T> {
    type Output = MatrixRef<'b, T>;

    fn get(&'b self, (rows, cols): (R, S)) -> Self::Output {
        MatrixRef {
            rows_begin: get_lower_index(&rows, self.rows()) + self.rows_begin,
            rows_end: get_upper_index(&rows, self.rows()) + self.rows_begin,
            cols_begin: get_lower_index(&cols, self.cols()) + self.cols_begin,
            cols_end: get_upper_index(&cols, self.cols()) + self.cols_begin,
            matrix: self.matrix
        }
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>, S: RangeBounds<usize>> Indexed<'b, (R, S)> for MatrixRefMut<'a, T> {
    type Output = MatrixRef<'b, T>;

    fn get(&'b self, (rows, cols): (R, S)) -> Self::Output {
        MatrixRef {
            rows_begin: get_lower_index(&rows, self.rows()) + self.rows.start,
            rows_end: get_upper_index(&rows, self.rows()) + self.rows.start,
            cols_begin: get_lower_index(&cols, self.cols()) + self.cols.start,
            cols_end: get_upper_index(&cols, self.cols()) + self.cols.end,
            matrix: self.matrix
        }
    }
}

impl<'a, 'b, T: 'b, R: RangeBounds<usize>, S: RangeBounds<usize>> IndexedMut<'b, (R, S)> for MatrixRefMut<'a, T> {
    type Output = MatrixRefMut<'b, T>;

    fn get_mut(&'b mut self, (rows, cols): (R, S)) -> Self::Output {
        MatrixRefMut {
            rows: (get_lower_index(&rows, self.rows()) + self.rows.start)..(get_upper_index(&rows, self.rows()) + self.rows.start),
            cols: (get_lower_index(&cols, self.cols()) + self.cols.start)..(get_upper_index(&cols, self.cols()) + self.cols.start),
            matrix: self.matrix
        }
    }
}

impl<'b, T: 'b> Indexed<'b, usize> for Matrix<T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, row: usize) -> Self::Output {
        VectorRef::create(&self[row])
    }
}

impl<'b, T: 'b> IndexedMut<'b, usize> for Matrix<T> {
    type Output = VectorRefMut<'b, T>;

    fn get_mut(&'b mut self, row: usize) -> Self::Output {
        VectorRefMut::create(&mut self[row])
    }
}

impl<'a, 'b, T: 'b> Indexed<'b, usize> for MatrixRef<'a, T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, row: usize) -> Self::Output {
        VectorRef::create(&self[row])
    }
}

impl<'a, 'b, T: 'b> Indexed<'b, usize> for MatrixRefMut<'a, T> {
    type Output = VectorRef<'b, T>;

    fn get(&'b self, row: usize) -> Self::Output {
        VectorRef::create(&self[row])
    }
}

impl<'a, 'b, T: 'b> IndexedMut<'b, usize> for MatrixRefMut<'a, T> {
    type Output = VectorRefMut<'b, T>;

    fn get_mut(&'b mut self, row: usize) -> Self::Output {
        VectorRefMut::create(&mut self[row])
    }
}

impl<T> Matrix<T>
where
    T: Add<T, Output = T> + Copy + Mul<T, Output = T>,
{
    ///
    /// Let T be the identity matrix (mxm where this matrix is mxn), in which the entries
    /// [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the values in transform.
    /// This function performs the multiplication A' := T * A, where A is this matrix
    ///
    pub fn transform_two_dims_left(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        self.as_mut().transform_two_dims_left(fst, snd, transform)
    }

    ///
    /// Let T be the identity matrix (nxn where this matrix is mxn), in which
    /// the entries [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the
    /// values in transform.
    /// This function performs the multiplication A' := A * T, where A is this matrix
    ///
    pub fn transform_two_dims_right(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        self.as_mut().transform_two_dims_right(fst, snd, transform)
    }
}

impl<'a, T> MatrixRefMut<'a, T>
where
    T: Add<T, Output = T> + Copy + Mul<T, Output = T>,
{
    ///
    /// Let T be the identity matrix (mxm where this matrix is mxn), in which the entries
    /// [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the values in transform.
    /// This function performs the multiplication A' := T * A, where A is this matrix
    ///
    pub fn transform_two_dims_left(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        assert!(fst < snd);
        self.assert_row_in_range(fst);
        self.assert_row_in_range(snd);
        for col in 0..self.cols() {
            let b = self[fst][col];
            self[fst][col] = self[fst][col] * transform[0] + self[snd][col] * transform[1];
            self[snd][col] = b * transform[2] + self[snd][col] * transform[3];
        }
    }

    ///
    /// Let T be the identity matrix (nxn where this matrix is mxn), in which
    /// the entries [fst,fst], [fst, snd], [snd, fst], [snd, snd] are replaced by the
    /// values in transform.
    /// This function performs the multiplication A' := A * T, where A is this matrix
    ///
    pub fn transform_two_dims_right(&mut self, fst: usize, snd: usize, transform: &[T; 4]) {
        assert!(fst < snd);
        self.assert_col_in_range(fst);
        self.assert_col_in_range(snd);
        for row in 0..self.rows() {
            let b = self[row][fst];
            self[row][fst] = self[row][fst] * transform[0] + self[row][snd] * transform[2];
            self[row][snd] = b * transform[1] + self[row][snd] * transform[3];
        }
    }
}

impl<T, M> Mul<M> for Matrix<T>
where
    T: Add<T, Output = T> + Clone + Mul<T, Output = T>,
    M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: M) -> Self::Output {
        self.as_ref().mul(rhs)
    }
}

impl<T, M> Add<M> for Matrix<T>
    where T: Add<T, Output = T> + Clone, M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn add(self, rhs: M) -> Self::Output {
        self += rhs;
        return self;
    }
}

impl<T, M> Sub<M> for Matrix<T>
    where T: Sub<T, Output = T> + Clone, M: MatrixView<T> 
{
    type Output = Matrix<T>;

    fn sub(self, rhs: M) -> Self::Output {
        self -= rhs;
        return self;
    }
}

impl<T> Neg for Matrix<T>
    where T: Neg + Clone
{
    type Output = Matrix<T::Output>;

    fn neg(self) -> Self::Output {
        self.as_ref().neg()
    }
}

impl<'a, T, M> Mul<M> for MatrixRef<'a, T>
where
    T: Add<T, Output = T> + Clone + Mul<T, Output = T>,
    M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: M) -> Self::Output {
        assert_eq!(self.cols(), rhs.rows());
        let cols = rhs.cols();
        let data: Vec<T> = (0..(cols * self.rows()))
            .map(|index: usize| {
                let row = index / cols;
                let col = index % cols;
                (1..self.cols())
                    .map(|k: usize| self[row][k].clone() * rhs.at(k, col).clone())
                    .fold(
                        self[row][0].clone() * rhs.at(0, col).clone(),
                        |acc: T, el: T| acc + el,
                    )
            })
            .collect();
        Matrix::new(data.into_boxed_slice(), self.rows())
    }
}

impl<'a, T, M> Add<M> for MatrixRef<'a, T>
    where T: Add<T, Output = T> + Clone, M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn add(self, rhs: M) -> Self::Output {
        Matrix::new((0..self.rows()).flat_map(|row| (0..self.cols()).map(move |col| *self.at(row, col) + *rhs.at(row, col))).collect::<Vec<T>>().into_boxed_slice(), self.rows())
    }
}

impl<'a, T, M> Sub<M> for MatrixRef<'a, T>
    where T: Sub<T, Output = T> + Clone, M: MatrixView<T> 
{
    type Output = Matrix<T>;

    fn sub(self, rhs: M) -> Self::Output {
        Matrix::new((0..self.rows()).flat_map(|row| (0..self.cols()).map(move |col| *self.at(row, col) - *rhs.at(row, col))).collect::<Vec<T>>().into_boxed_slice(), self.rows())
    }
}

impl<'a, T> Neg for MatrixRef<'a, T>
    where T: Neg + Clone
{
    type Output = Matrix<T::Output>;

    fn neg(self) -> Self::Output {
        Matrix::new((0..self.rows()).flat_map(|row| (0..self.cols()).map(move |col| -self.at(row, col).clone())).collect::<Vec<T::Output>>().into_boxed_slice(), self.rows())
    }
}

impl<'a, T, M> Mul<M> for MatrixRefMut<'a, T>
where
    T: Add<T, Output = T> + Clone + Mul<T, Output = T>,
    M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn mul(self, rhs: M) -> Self::Output {
        self.into_const().mul(rhs)
    }
}

impl<'a, T, M> Add<M> for MatrixRefMut<'a, T>
    where T: Add<T, Output = T> + Clone, M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn add(self, rhs: M) -> Self::Output {
        self.into_const().add(rhs)
    }
}

impl<'a, T, M> Sub<M> for MatrixRefMut<'a, T>
    where T: Sub<T, Output = T> + Clone, M: MatrixView<T> 
{
    type Output = Matrix<T>;

    fn sub(self, rhs: M) -> Self::Output {
        self.into_const().sub(rhs)
    }
}

impl<'a, T> Neg for MatrixRefMut<'a, T>
    where T: Neg + Clone
{
    type Output = Matrix<T::Output>;

    fn neg(self) -> Self::Output {
        self.into_const().neg()
    }
}

impl< T, M> AddAssign<M> for Matrix<T> 
    where T: AddAssign + Clone, M: MatrixView<T>
{
    fn add_assign(&mut self, rhs: M) {
        self.as_mut().add_assign(rhs);
    }
}

impl<T, M> SubAssign<M> for Matrix<T> 
    where T: SubAssign + Clone, M: MatrixView<T>
{
    fn sub_assign(&mut self, rhs: M) {
        self.as_mut().sub_assign(rhs);
    }
}


impl<'a, T, M> AddAssign<M> for MatrixRefMut<'a, T> 
    where T: AddAssign + Clone, M: MatrixView<T>
{
    fn add_assign(&mut self, rhs: M) {
        assert_sizes_match(self, &rhs);
        for row in 0..self.rows() {
            for col in 0..self.cols() {
                self[row][col].add_assign(rhs.at(row, col).clone());
            }
        }
    }
}

impl<'a, T, M> SubAssign<M> for MatrixRefMut<'a, T> 
    where T: SubAssign + Clone, M: MatrixView<T>
{
    fn sub_assign(&mut self, rhs: M) {
        assert_sizes_match(self, &rhs);
        for row in 0..self.rows() {
            for col in 0..self.cols() {
                self[row][col].sub_assign(rhs.at(row, col).clone());
            }
        }
    }
}

#[test]
fn test_matrix_get_rows() {
    #[rustfmt::skip]
	let mut m = Matrix::new(Box::new([1,  2,  3,
	                                  4,  5,  6,
								      7,  8,  9,
								      10, 11, 12]), 4);
    assert_eq!(3, m.cols());
    assert_eq!(4, m.rows());
    {
        let (fst_row, snd_row) = m.get_rows(0, 2);
        assert_eq!(&[1, 2, 3], fst_row.into_slice());
        assert_eq!(&[7, 8, 9], snd_row.into_slice());
    }
    {
        let (fst_row, snd_row) = m.get_rows(3, 2);
        assert_eq!(&[10, 11, 12], fst_row.into_slice());
        assert_eq!(&[7, 8, 9], snd_row.into_slice());
    }
}

#[test]
fn test_matrix_submatrix() {
    #[rustfmt::skip]
	let mut m = Matrix::new(Box::new([1,  2,  3,  7,
	                                  4,  5,  6,  11,
								      7,  8,  9,  2,
								      10, 11, 12, 4]), 4);
    assert_eq!(4, m.cols());
    assert_eq!(4, m.rows());

    let mut n = m.get_mut((1..3, 1..3));
    assert_eq!(5, n[0][0]);
    assert_eq!(9, n[1][1]);
    assert_eq!(2, n.rows());
    assert_eq!(2, n.cols());

    {
        let (mut r1, r2) = n.get_rows(1, 0);
        r1 += r2.as_const();
    }

    assert_eq!(7, m[2][0]);
    assert_eq!(13, m[2][1]);
    assert_eq!(15, m[2][2]);
    assert_eq!(2, m[2][3]);

    assert_eq!(2, m[0][1]);
    assert_eq!(5, m[1][1]);
    assert_eq!(11, m[3][1]);
}

#[test]
fn test_matrix_transform_two_dims_left() {
    #[rustfmt::skip]
	let mut m = Matrix::new(Box::new([1., 2., 3.,
	                                  4., 5., 6.,
								      7., 8., 9.]), 3);
    m.get_mut((.., ..))
        .transform_two_dims_left(0, 2, &[0., 1., 1., 0.]);

    #[rustfmt::skip]
	assert_eq!(&[7., 8., 9.,
	             4., 5., 6.,
				 1., 2., 3.], m.data());
    m.get_mut((.., ..))
        .transform_two_dims_left(1, 2, &[0.5, 0.5, 1.0, 0.5]);

    #[rustfmt::skip]
	assert_eq!(&[7.,  8.,  9.,
	             2.5, 3.5, 4.5,
				 4.5, 6.0, 7.5], m.data());
}

#[test]
fn test_matmul() {
    #[rustfmt::skip]
    let a = Matrix::new(Box::new([1, 2,
                                  3, 2]), 2);

    #[rustfmt::skip]
    let b = Matrix::new(Box::new([1, 2, 3,
                                  3, 4, 2]), 2);

    #[rustfmt::skip]
    assert_eq!(&[7, 10, 7, 
                 9, 14, 13], (a.get((.., ..)) * b.get((.., ..))).data());
}

#[test]
fn test_invert() {
    #[rustfmt::skip]
    let a = Matrix::new(Box::new([1., 2.,
                                  2., 0.]), 2);

    #[rustfmt::skip]
    assert_eq!(&[0.,  0.5,
                 0.5, -0.25], a.get((.., ..)).invert().unwrap().data());
}
