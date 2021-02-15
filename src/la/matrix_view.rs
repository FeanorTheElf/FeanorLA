use super::vector_view::*;

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
}

pub trait FromFnCreateable<T> : MatrixView<T> {

    fn from_fn<F>(rows: usize, cols: usize, f: F) -> Self
        where F: FnMut(usize, usize) -> T;
}

pub trait MatrixViewMut<T>: MatrixView<T> {

    fn at_mut(&mut self, row: usize, col: usize) -> &mut T;
    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize));
}

pub trait LifetimeMatrixMutRowIter<'a, T>: MatrixViewMut<T> {

    type RowRef: VectorViewMut<T>;
    type RowIter: Iterator<Item = Self::RowRef>;

    fn rows_mut(&'a mut self) -> Self::RowIter;

    fn get_row_mut(&'a mut self, index: usize) -> Self::RowRef {
        self.assert_row_in_range(index);
        self.rows_mut().nth(index).unwrap()
    }
}

pub trait MatrixMutRowIter<T>: for<'a> LifetimeMatrixMutRowIter<'a, T> {}
