use super::matrix_view::*;
use super::vector_view::*;
use super::vec::*;
use std::marker::PhantomData;

pub struct ColumnVector<V, T>
    where V: VectorView<T>
{
    base: Vector<V, T>,
    element: PhantomData<T>
}

impl<V, T> ColumnVector<V, T>
    where V: VectorView<T>
{
    pub fn new(vector: Vector<V, T>) -> Self {
        ColumnVector {
            base: vector,
            element: PhantomData
        }
    }
}

impl<V, T> MatrixView<T> for ColumnVector<V, T>
    where V: VectorView<T>
{
    fn row_count(&self) -> usize {
        self.base.len()
    }

    fn col_count(&self) -> usize {
        1
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.base.at(row)
    }
}

impl<V, T> MatrixViewMut<T> for ColumnVector<V, T>
    where V: VectorViewMut<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.base.at_mut(row)
    }

    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        self.assert_col_in_range(fst.1);
        self.assert_col_in_range(snd.1);
        unimplemented!();
    }
}