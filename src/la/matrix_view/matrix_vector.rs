use super::*;
use super::super::vector_view::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct ColumnVector<V, T>
    where V: VectorView<T>
{
    base: V,
    element: PhantomData<T>
}

impl<V, T> Copy for ColumnVector<V, T>
    where V: VectorView<T> + Copy
{}

impl<V, T> Clone for ColumnVector<V, T>
    where V: VectorView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<V, T> ColumnVector<V, T>
    where V: VectorView<T>
{
    pub fn new(vector: V) -> Self {
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
        self.base.swap(fst.0, snd.0);
    }
}

#[derive(Debug)]
pub struct RowVector<V, T>
    where V: VectorView<T>
{
    base: V,
    element: PhantomData<T>
}

impl<V, T> Copy for RowVector<V, T>
    where V: VectorView<T> + Copy
{}

impl<V, T> Clone for RowVector<V, T>
    where V: VectorView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<V, T> RowVector<V, T>
    where V: VectorView<T>
{
    pub fn new(vector: V) -> Self {
        RowVector {
            base: vector,
            element: PhantomData
        }
    }
}

impl<V, T> MatrixView<T> for RowVector<V, T>
    where V: VectorView<T>
{
    fn col_count(&self) -> usize {
        self.base.len()
    }

    fn row_count(&self) -> usize {
        1
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.base.at(col)
    }
}

impl<V, T> MatrixViewMut<T> for RowVector<V, T>
    where V: VectorViewMut<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        self.base.at_mut(row)
    }

    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        self.assert_row_in_range(fst.0);
        self.assert_row_in_range(snd.0);
        self.assert_col_in_range(fst.1);
        self.assert_col_in_range(snd.1);
        self.base.swap(fst.0, snd.0)
    }
}