pub mod constant_value_matrix;
pub mod matrix_vector;
pub mod matrix_transpose;
pub mod submatrix;
pub mod diagonal_matrix;

pub trait MatrixView<T>: Sized {

    fn row_count(&self) -> usize;
    fn col_count(&self) -> usize;
    fn at(&self, row: usize, col: usize) -> &T;

    fn assert_row_in_range(&self, row: usize) {
        assert!(row < self.row_count(), "Row index {} out of range 0..{}", row, self.row_count());
    }

    fn assert_col_in_range(&self, col: usize) {
        assert!(col < self.col_count(), "Column index {} out of range 0..{}", col, self.col_count());
    }

    fn into_owned(self) -> MatrixOwned<T> 
        where T: Clone
    {
        MatrixOwned::from_fn(self.row_count(), self.col_count(), |i, j| self.at(i, j).clone())
    }
}

pub trait FromFnCreateable<T> : MatrixView<T> {

    fn from_fn<F>(rows: usize, cols: usize, f: F) -> Self
        where F: FnMut(usize, usize) -> T;
}

pub trait MatrixViewMut<T>: MatrixView<T> {

    fn at_mut(&mut self, row: usize, col: usize) -> &mut T;
    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize));
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixOwned<T> {
    col_count: usize,
    row_count: usize,
    data: Box<[T]>,
}

impl<T> MatrixView<T> for MatrixOwned<T> {
    
    fn row_count(&self) -> usize {
        self.row_count
    }

    fn col_count(&self) -> usize {
        self.col_count
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_col_in_range(col);
        self.assert_row_in_range(row);
        &self.data[row * self.col_count + col]
    }

    fn into_owned(self) -> MatrixOwned<T> {
        self
    }
}

impl<T> MatrixViewMut<T> for MatrixOwned<T> {
    
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_col_in_range(col);
        self.assert_row_in_range(row);
        &mut self.data[row * self.col_count() + col]
    }

    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        self.assert_row_in_range(fst.0);
        self.assert_row_in_range(snd.0);
        self.assert_col_in_range(fst.1);
        self.assert_col_in_range(snd.1);
        if fst == snd {
            return;
        }
        self.data.swap(fst.1 + fst.0 * self.col_count(), snd.1 + snd.0 * self.col_count());
    }
}


impl<T> MatrixOwned<T> {

    pub fn from_data(data: Box<[T]>, rows: usize, cols: usize) -> MatrixOwned<T> {
        if cols == 0 {
            assert!(
                data.len() == 0,
                "A zero-column matrix must have no data elements"
            );
        } else {
            assert!(
                data.len() % cols == 0,
                "Data length must be a multiple of column count, but got {} and {}",
                data.len(),
                cols
            );
        }
        assert!(data.len() == cols * rows);
        MatrixOwned {
            data: data,
            col_count: cols,
            row_count: rows
        }
    }

    pub fn from_array<const R: usize, const C: usize>(
        array: [[T; C]; R]
    ) -> MatrixOwned<T> 
    {
        let data = <[[T; C]; R] as std::iter::IntoIterator>::into_iter(array).flat_map(|row| 
            <[T; C] as std::iter::IntoIterator>::into_iter(row)
        ).collect::<Vec<T>>().into_boxed_slice();
        Self::from_data(data, R, C)
    }

    ///
    /// Converts the data stored in self into an iterator that
    /// yields all elements, row by row.
    /// 
    pub fn into_data_iter(self) -> std::vec::IntoIter<T> {
        self.data.into_vec().into_iter()
    }
}

impl<T> FromFnCreateable<T> for MatrixOwned<T> {

    fn from_fn<F>(rows: usize, cols: usize, mut f: F) -> Self
        where F: FnMut(usize, usize) -> T
    {
        let mut data = Vec::with_capacity(rows * cols);
        for row in 0..rows {
            data.extend((0..cols).map(|col| f(row, col)));
        }
        return Self::from_data(data.into_boxed_slice(), rows, cols);
    }
}

impl<'a, T, M> MatrixView<T> for &'a M
    where M: MatrixView<T>
{
    fn row_count(&self) -> usize {
        (*self).row_count()
    }

    fn col_count(&self) -> usize {
        (*self).col_count()
    }

    fn at(&self, row: usize, col: usize) -> &T {
        (*self).at(row, col)
    }
}

impl<'a, T, M> MatrixView<T> for &'a mut M
    where M: MatrixView<T>
{
    fn row_count(&self) -> usize {
        (**self).row_count()
    }

    fn col_count(&self) -> usize {
        (**self).col_count()
    }

    fn at(&self, row: usize, col: usize) -> &T {
        (**self).at(row, col)
    }
}

impl<'a, T, M> MatrixViewMut<T> for &'a mut M
    where M: MatrixViewMut<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        (*self).at_mut(row, col)
    }

    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        (*self).swap(fst, snd);
    }
}

#[test]
fn test_from_fn() {
    let a = MatrixOwned::from_array([[0, 1], [1, 2]]);
    let b = MatrixOwned::from_fn(2, 2, |i, j| i + j);
    assert_eq!(a, b);
}
