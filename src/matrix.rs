use super::alg::*;
use super::indexed::{Indexed, IndexedMut};
use std::mem::swap;
use std::ops::{
    Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, RangeBounds, Bound, Range,
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

impl<T> MatrixView<T> for Vector<T>
{
    fn rows(&self) -> usize {
        self.len()
    }

    fn cols(&self) -> usize {
        1
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        <Self as VectorView<T>>::at(self, row)
    }
}

impl<'a, T> MatrixView<T> for VectorRef<'a, T>
{
    fn rows(&self) -> usize {
        self.len()
    }

    fn cols(&self) -> usize {
        1
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        <Self as VectorView<T>>::at(self, row)
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

impl<T> MatrixViewMut<T> for Vector<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        <Self as VectorViewMut<T>>::at_mut(self, row)
    }
}

///
/// Represents a mxn matrix with elements of type T. Typical matrix operations
/// are not optimized much, so this type is only suitable for small matrices.
/// Instead, the focus lies on a convenient interface and a generic design.
///
#[derive(Debug, Clone)]
pub struct Matrix<T> {
    cols: usize,
    data: Box<[T]>,
}

#[derive(Debug)]
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

impl<'a, T> Clone for MatrixRef<'a, T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, T> Copy for MatrixRef<'a, T> {

}

impl<T> MatrixView<T> for Matrix<T> {
    
    fn rows(&self) -> usize {
        self.data.len() / self.cols
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        &self.data[row * self.cols + col]
    }
}

impl<T> MatrixViewMut<T> for Matrix<T> {
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        &mut self.data[row * self.cols + col]
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
        self.rows.end - self.rows.start
    }

    fn cols(&self) -> usize {
        self.cols.end - self.cols.start
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

    pub fn from_func<F>(rows: usize, cols: usize, mut f: F) -> Matrix<T> 
        where F: FnMut(usize, usize) -> T
    {
        let mut data = Vec::with_capacity(rows * cols);
        for row in 0..rows {
            data.extend((0..cols).map(|col| f(row, col)));
        }
        return Self::from_data(data.into_boxed_slice(), cols);
    }

    pub fn from_data(data: Box<[T]>, cols: usize) -> Matrix<T> {
        assert!(data.len() > 0, "Cannot create matrix with zero elements");
        assert!(
            data.len() % cols == 0,
            "Data length must be a multiple of column count, but got {} and {}",
            data.len(),
            cols
        );
        Matrix {
            data: data,
            cols: cols
        }
    }

    pub fn from_array<const R: usize, const C: usize>(array: [[T; C]; R]) -> Matrix<T> 
        where [T; C]: std::array::LengthAtMost32, [[T; C]; R]: std::array::LengthAtMost32
    {
        let data = std::array::IntoIter::new(array).flat_map(|row| std::array::IntoIter::new(row)).collect::<Vec<T>>().into_boxed_slice();
        Self::from_data(data, C)
    }

    pub fn from_nocopy<U>(value: Matrix<U>) -> Self
    where
        T: From<U>,
    {
        let cols = value.cols();
        let data: Vec<T> = value
            .data
            .into_vec()
            .into_iter()
            .map(|d| T::from(d))
            .collect();
        return Matrix {
            data: data.into_boxed_slice(),
            cols: cols,
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

    pub fn into_column_vector(self) -> Vector<T> {
        assert_eq!(1, self.cols(), "Matrix has {} columns, and therefore cannot be converted into a column vector", self.cols());
        Vector::new(self.data)
    }
}

impl<T> Matrix<T>
    where T: Clone 
{
    pub fn copy_of<M: MatrixView<T>>(value: M) -> Matrix<T> {
        Self::from_func(value.rows(), value.cols(), |row, col| value.at(row, col).clone())
    }
}

impl<T> Matrix<T>
    where T: Clone + AddAssign + Mul<Output = T>
{
    pub fn frobenius_square(&self) -> T {
        self.as_ref().frobenius_square()
    }
}

impl<T> Matrix<T>
    where T: Clone + MulAssign
{
    pub fn scal(&mut self, scal: T) {
        self.as_mut().scal(scal);
    }
}

impl<'a, T> MatrixRef<'a, T>
    where T: Clone + AddAssign + Mul<Output = T>
{
    pub fn frobenius_square(&self) -> T {
        let mut it = self.into_iter().map(|(x, _, _)| x.clone() * x.clone());
        let mut result = it.next().unwrap();
        for x in it {
            result += x;
        }
        return result;
    }
}

impl<'a, T> MatrixRefMut<'a, T>
    where T: Clone + AddAssign + Mul<Output = T>
{
    pub fn frobenius_square(&self) -> T {
        self.as_const().frobenius_square()
    }
}

impl<'a, T> MatrixRefMut<'a, T>
    where T: Clone + MulAssign
{
    pub fn scal(&mut self, scal: T) {
        for (x, _, _) in self.get_mut((.., ..)).into_iter() {
            *x *= scal.clone();
        }
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
    
    pub fn as_const(&self) -> MatrixRef<T> {
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
        Matrix::from_func(rows, cols, |_, _| T::zero())
    }
}

impl<'a, T, U: 'a> From<MatrixRef<'a, U>> for Matrix<T>
where
    T: From<&'a U>,
{
    fn from(value: MatrixRef<'a, U>) -> Self {
        let value_ref = &value;
        let data: Vec<T> = (value.rows_begin..value.rows_end)
            .flat_map(|row| {
                (value.cols_begin..value.cols_end)
                    .map(move |col| value_ref.matrix.at(row, col))
                    .map(|d| T::from(d))
            })
            .collect();
        return Matrix {
            data: data.into_boxed_slice(),
            cols: value.cols(),
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

impl<T> Matrix<T>
where
    T: Clone + Field
{
    ///
    /// Calculates the inverse of this matrix. Use only for small matrices, this
    /// is just simple gaussian elimination, and is neither very performant nor
    /// numerically stable!
    ///
    pub fn invert(&self) -> Result<Matrix<T>, ()> {
        self.as_ref().invert()
    }
}

impl<'a, T> MatrixRef<'a, T>
where
    T: Clone + Field
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
            // search for a non-null entry
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

impl<'a, T> MatrixRefMut<'a, T>
where
    T: Clone + Field
{
    ///
    /// Calculates the inverse of this matrix. Use only for small matrices, this
    /// is just simple gaussian elimination, and is neither very performant nor
    /// numerically stable!
    ///
    pub fn invert(&self) -> Result<Matrix<T>, ()> {
        self.as_const().invert()
    }
}

impl<'a, T> MatrixRefMut<'a, T> {
    
    pub fn get_rows<'b>(&'b mut self, fst: usize, snd: usize) -> (VectorRefMut<'b, T>, VectorRefMut<'b, T>) {
        self.assert_row_in_range(fst);
        self.assert_row_in_range(snd);
        let (fst_row, snd_row) = self.matrix.get_rows(fst + self.rows.start, snd + self.rows.start);
        return (fst_row.into_subrange(self.cols.start..self.cols.end), snd_row.into_subrange(self.cols.start..self.cols.end))
    }

    pub fn swap_rows(&mut self, fst: usize, snd: usize) {
        if fst != snd {
            self.assert_row_in_range(fst);
            self.assert_row_in_range(snd);
            let cols = self.cols();
            let (mut fst_row, mut snd_row) = self.get_rows(fst, snd);
            for col in 0..cols {
                swap(&mut fst_row[col], &mut snd_row[col]);
            }
        }
    }

    pub fn swap_cols(&mut self, fst: usize, snd: usize) {
        if fst != snd {
            self.assert_col_in_range(fst);
            self.assert_col_in_range(snd);
            for row in 0..self.rows() {
                self[row].swap(fst, snd);
            }
        }
    }

    pub fn assign(&mut self, rhs: Matrix<T>) {
        for (value, row, col) in rhs.into_iter() {
            self[row][col] = value;
        }
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        self.assert_row_in_range(index);
        &self.data[(index * self.cols)..(index * self.cols + self.cols)]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {

    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.assert_row_in_range(index);
        &mut self.data[(index * self.cols)..(index * self.cols + self.cols)]
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

fn get_lower_index<R: RangeBounds<usize>>(range: &R, _: usize) -> usize {
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
            cols_end: get_upper_index(&cols, self.cols()) + self.cols.start,
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
    where T: AddAssign<T> + Clone, M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn add(mut self, rhs: M) -> Self::Output {
        self += rhs;
        return self;
    }
}

impl<T, M> Sub<M> for Matrix<T>
    where T: SubAssign<T> + Clone, M: MatrixView<T> 
{
    type Output = Matrix<T>;

    fn sub(mut self, rhs: M) -> Self::Output {
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
        Matrix::from_func(self.rows(), rhs.cols(), |row, col|{
            (1..self.cols())
            .map(|k: usize| self[row][k].clone() * rhs.at(k, col).clone())
            .fold(
                self[row][0].clone() * rhs.at(0, col).clone(),
                |acc: T, el: T| acc + el,
            )
        })
    }
}

impl<'a, T, M> Add<M> for MatrixRef<'a, T>
    where T: Add<T, Output = T> + Clone, M: MatrixView<T>
{
    type Output = Matrix<T>;

    fn add(self, rhs: M) -> Self::Output {
        assert_sizes_match(&self, &rhs);
        Matrix::from_func(self.rows(), self.cols(), |row, col| self.at(row, col).clone() + rhs.at(row, col).clone())
    }
}

impl<'a, T, M> Sub<M> for MatrixRef<'a, T>
    where T: Sub<T, Output = T> + Clone, M: MatrixView<T> 
{
    type Output = Matrix<T>;

    fn sub(self, rhs: M) -> Self::Output {
        assert_sizes_match(&self, &rhs);
        Matrix::from_func(self.rows(), self.cols(), |row, col| self.at(row, col).clone() - rhs.at(row, col).clone())
    }
}

impl<'a, T> Neg for MatrixRef<'a, T>
    where T: Neg + Clone
{
    type Output = Matrix<T::Output>;

    fn neg(self) -> Self::Output {
        Matrix::from_func(self.rows(), self.cols(), |row, col| -self.at(row, col).clone())
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

pub struct MatrixIter<I> 
{
    data: I,
    row: usize,
    col: usize,
    cols: usize,
    skip: usize
}

impl<I> Iterator for MatrixIter<I> 
    where I: Iterator + ExactSizeIterator
{
    type Item = (I::Item, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(value) = self.data.next() {
            let result = (value, self.row, self.col);
            self.col += 1;
            if self.col >= self.cols {
                self.row += 1;
                self.cols = 0;
                for _ in 0..self.skip {
                    self.data.next();
                }
            }
            return Some(result);
        } else {
            return None;
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        debug_assert!((self.data.size_hint().0 - self.cols) % (self.cols + self.skip) == 0);
        let rows = (self.data.size_hint().0 - self.cols) / (self.cols + self.skip) + 1;
        let size = rows * self.cols - self.row * self.cols - self.col;
        return (size, Some(size));
    }
}

impl<I> ExactSizeIterator for MatrixIter<I> 
    where I: Iterator + ExactSizeIterator
{
}

impl<T> IntoIterator for Matrix<T> {
    type Item = (T, usize, usize);
    type IntoIter = MatrixIter<std::vec::IntoIter<T>>;

    fn into_iter(self) -> Self::IntoIter {
        let cols = self.cols();
        MatrixIter {
            data: self.data.into_vec().into_iter(),
            row: 0,
            col: 0,
            cols: cols,
            skip: 0
        }
    }
}

impl<'a, T> IntoIterator for MatrixRef<'a, T> {
    type Item = (&'a T, usize, usize);
    type IntoIter = MatrixIter<std::slice::Iter<'a, T>>;

    fn into_iter(self) -> Self::IntoIter {
        MatrixIter {
            data: self.matrix.data[(self.cols_begin + self.rows_begin * self.matrix.cols)..(self.cols_end + (self.rows_end - 1) * self.matrix.cols)].iter(),
            row: 0,
            col: 0,
            cols: self.cols(),
            skip: self.matrix.cols() - self.cols()
        }
    }
}

impl<'a, T> IntoIterator for MatrixRefMut<'a, T> {
    type Item = (&'a mut T, usize, usize);
    type IntoIter = MatrixIter<std::slice::IterMut<'a, T>>;

    fn into_iter(self) -> Self::IntoIter {
        let cols = self.cols();
        let parent_cols = self.matrix.cols();
        MatrixIter {
            data: self.matrix.data[(self.cols.start + self.rows.start * parent_cols)..(self.cols.end + (self.rows.end - 1) * parent_cols)].iter_mut(),
            row: 0,
            col: 0,
            cols: cols,
            skip: parent_cols - cols
        }
    }
}

impl<T, M> PartialEq<M> for Matrix<T>
    where T: PartialEq<T>, M: MatrixView<T>
{
    fn eq(&self, rhs: &M) -> bool {
        self.as_ref() == *rhs
    }
}

impl<'a, T, M> PartialEq<M> for MatrixRef<'a, T>
    where T: PartialEq<T>, M: MatrixView<T>
{
    fn eq(&self, rhs: &M) -> bool {
        assert_sizes_match(self, rhs);
        for row in 0..self.rows() {
            for col in 0..self.cols() {
                if self.at(row, col) != rhs.at(row, col) {
                    return false;
                }
            }
        }
        return true;
    }
}

impl<'a, T, M> PartialEq<M> for MatrixRefMut<'a, T>
    where T: PartialEq<T>, M: MatrixView<T>
{
    fn eq(&self, rhs: &M) -> bool {
        self.as_const() == *rhs
    }
}

impl<T, M> MulAssign<M> for Matrix<T> 
    where T: Add<T, Output = T> + Clone + Mul<T, Output = T>, M: MatrixView<T>
{
    fn mul_assign(&mut self, rhs: M) {
        self.as_mut().mul_assign(rhs);
    }
}

impl<'a, T, M> MulAssign<M> for MatrixRefMut<'a, T> 
    where T: Add<T, Output = T> + Clone + Mul<T, Output = T>, M: MatrixView<T>
{
    fn mul_assign(&mut self, rhs: M) {
        self.assign(self.as_const() * rhs);
    }
}

macro_rules! matlab {
    ($($($expr:expr),*);*) => {
        Matrix::from_array([$([$($expr),*]),*])
    };
    ($($($expr:expr)*);*) => {
        Matrix::from_array([$([$($expr),*]),*])
    };
}

#[test]
fn test_matlab_macro() {
    let m = matlab![1, 2, 3; 4, 5, 6];
    assert_eq!(1, m[0][0]);
    assert_eq!(5, m[1][1]);
    let m2 = matlab![1 2 3; 4 5 6];
    assert_eq!(1, m2[0][0]);
    assert_eq!(5, m2[1][1]);
}

#[test]
fn test_matrix_index() {
    #[rustfmt::skip]
    let mut m = Matrix::from_array([[1, 2, 3],
                                    [4, 5, 6]]);
    println!("{:?}", m);
    assert_eq!(5, m[1][1]);
    assert_eq!(5, m.as_ref()[1][1]);
    assert_eq!(5, m.as_mut()[1][1]);

    let mut sub = m.get_mut((.., 1..));
    assert_eq!(6, sub[1][1]);
    assert_eq!(6, sub.index_mut(1)[1]);
    assert_eq!(6, sub.as_const()[1][1]);

    let mut subsub = sub.get_mut((..=0, 1..2));
    assert_eq!(3, subsub[0][0]);
    assert_eq!(3, subsub.index_mut(0)[0]);
    assert_eq!(3, subsub.as_const()[0][0]);
}

#[test]
fn test_matrix_get_rows() {
    #[rustfmt::skip]
	let mut m = Matrix::from_array([[1,  2,  3],
	                                 [4,  5,  6],
								     [7,  8,  9],
								     [10, 11, 12]]);
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
	let mut m = Matrix::from_array([[1,  2,  3,  7],
	                                [4,  5,  6,  11],
								    [7,  8,  9,  2],
								    [10, 11, 12, 4]]);
    assert_eq!(4, m.cols());
    assert_eq!(4, m.rows());

    let mut n = m.get_mut((1..3, 1..3));
    assert_eq!(5, n[0][0]);
    assert_eq!(9, n[1][1]);
    assert_eq!(2, n.rows());
    assert_eq!(2, n.cols());

    {
        let (mut r1, r2) = n.get_rows(1, 0);
        println!("{:?}", r1);
        r1 += r2.as_const();
        println!("{:?}", r1);
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
	let mut m = Matrix::from_array([[1., 2., 3.],
	                                [4., 5., 6.],
								    [7., 8., 9.]]);
    m.get_mut((.., ..))
        .transform_two_dims_left(0, 2, &[0., 1., 1., 0.]);

    #[rustfmt::skip]
	assert_eq!(&Matrix::from_array([[7., 8., 9.],
	                                [4., 5., 6.],
				                    [1., 2., 3.]]), &m);
    m.get_mut((.., ..))
        .transform_two_dims_left(1, 2, &[0.5, 0.5, 1.0, 0.5]);

    #[rustfmt::skip]
	assert_eq!(&Matrix::from_array([[7.,  8.,  9. ],
	                                [2.5, 3.5, 4.5],
				                    [4.5, 6.0, 7.5]]), &m);
}

#[test]
fn test_matmul() {
    #[rustfmt::skip]
    let a = Matrix::from_array([[1, 2],
                                [3, 2]]);

    #[rustfmt::skip]
    let b = Matrix::from_array([[1, 2, 3],
                                [3, 4, 2]]);

    let prod = a * b;
    #[rustfmt::skip]
    assert_eq!(&Matrix::from_array([[7, 10, 7], 
                                    [9, 14, 13]]), &prod);
}

#[test]
fn test_invert() {
    #[rustfmt::skip]
    let a = Matrix::from_array([[1., 2.],
                                [2., 0.]]);

    let a_inv = a.invert().unwrap();
    #[rustfmt::skip]
    assert_eq!(&Matrix::from_array([[0.,  0.5],
                                    [0.5, -0.25]]), &a_inv);
}
