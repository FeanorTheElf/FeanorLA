use super::vector_view::*;
use std::marker::PhantomData;

// Traits

pub struct MatrixRow<'a, T, M>
    where M: MatrixView<T> 
{
    row: usize,
    matrix: &'a M,
    item: PhantomData<T>
}

impl<'a, T, M> Copy for MatrixRow<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixRow<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, T, M> VectorView<T> for MatrixRow<'a, T, M> 
    where M: MatrixView<T>
{
    fn len(&self) -> usize {
        self.matrix.col_count()
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        self.matrix.at(self.row, i)
    }
}

pub struct MatrixRowIter<'a, T, M>
    where M: MatrixView<T> 
{
    current: MatrixRow<'a, T, M>
}

impl<'a, T, M> Iterator for MatrixRowIter<'a, T, M> 
    where M: MatrixView<T>
{
    type Item = MatrixRow<'a, T, M>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current;
        if result.row < result.matrix.row_count() {
            self.current.row += 1;
            return Some(result);
        } else {
            return None;
        }
    }
}

impl<'a, T, M> Copy for MatrixRowIter<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixRowIter<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}

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

    fn rows(&self) -> MatrixRowIter<T, Self> {
        MatrixRowIter {
            current: MatrixRow {
                row: 0,
                matrix: self,
                item: PhantomData
            }
        }
    }

    fn get_row(&self, row: usize) -> MatrixRow<T, Self> {
        self.assert_row_in_range(row);
        MatrixRow {
            row: row,
            matrix: self,
            item: PhantomData
        }
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
