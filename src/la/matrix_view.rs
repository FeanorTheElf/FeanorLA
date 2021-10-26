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
        let data = std::array::IntoIter::new(array).flat_map(|row| 
            std::array::IntoIter::new(row)
        ).collect::<Vec<T>>().into_boxed_slice();
        Self::from_data(data, R, C)
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

#[derive(Debug)]
pub struct OwnedMatrixRowMutRef<'b, T> {
    data: &'b mut [T]
}

impl<'b, T> VectorView<T> for OwnedMatrixRowMutRef<'b, T> {
    
    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<'b, T> VectorViewMut<T> for OwnedMatrixRowMutRef<'b, T> {
    
    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        self.data.swap(i, j);
    }
}

pub enum OwnedMatrixRowMutIter<'b, T> {
    RealRows(std::slice::ChunksExactMut<'b, T>),
    EmptyRows(usize)
}

trait EmptySlice: Sized {
    const EMPTY_SLICE: [Self; 0];
}

impl<T> EmptySlice for T {
    const EMPTY_SLICE: [T; 0] = [];
}

impl<'b, T> Iterator for OwnedMatrixRowMutIter<'b, T> {
    type Item = OwnedMatrixRowMutRef<'b, T>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            OwnedMatrixRowMutIter::RealRows(rows) => return rows.next().map(|r| OwnedMatrixRowMutRef {
                data: r
            }),
            OwnedMatrixRowMutIter::EmptyRows(count) => {
                if *count > 0 {
                    *count -= 1;

                    #[allow(const_item_mutation)]
                    return Some(OwnedMatrixRowMutRef {
                        data: &mut T::EMPTY_SLICE[..]
                    });
                } else {
                    return None;
                }
            }
        };
    }
}

impl<'b, T: 'b> LifetimeMatrixMutRowIter<'b, T> for MatrixOwned<T> {

    type RowRef = OwnedMatrixRowMutRef<'b, T>;
    type RowIter = OwnedMatrixRowMutIter<'b, T>;
    
    fn rows_mut(&'b mut self) -> OwnedMatrixRowMutIter<'b, T> {
        if self.col_count() != 0 {
            OwnedMatrixRowMutIter::RealRows(self.data.chunks_exact_mut(self.col_count()))
        } else {
            OwnedMatrixRowMutIter::EmptyRows(self.row_count)
        }
    }
}

impl<T: 'static> MatrixMutRowIter<T> for MatrixOwned<T> {}

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

impl<'a, 'b, T, M> LifetimeMatrixMutRowIter<'b, T> for &'a mut M
    where M: LifetimeMatrixMutRowIter<'b, T>
{
    type RowRef = M::RowRef;
    type RowIter = M::RowIter;

    fn rows_mut(&'b mut self) -> Self::RowIter {
        (*self).rows_mut()
    }
}

impl<'a, T, M> MatrixMutRowIter<T> for &'a mut M
    where M: MatrixMutRowIter<T>
{
}

#[test]
fn test_from_fn() {
    let a = MatrixOwned::from_array([[0, 1], [1, 2]]);
    let b = MatrixOwned::from_fn(2, 2, |i, j| i + j);
    assert_eq!(a, b);
}

#[test]
fn test_mut_row_iter() {
    let mut a = MatrixOwned::from_array([[0, 1], [3, 4]]);
    let mut it = a.rows_mut();
    let mut r1 = it.next().unwrap();
    let mut r2 = it.next().unwrap();
    assert!(it.next().is_none());
    assert_eq!(0, *r1.at_mut(0));
    assert_eq!(4, *r2.at_mut(1));
    *r1.at_mut(1) = 10;
    assert_eq!(10, *a.at(0, 1));
}