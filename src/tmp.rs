
pub struct OwnedMatrixRow<'b, T> {
    data: &'b [T]
}

impl<'b, T> Copy for OwnedMatrixRow<'b, T> {}

impl<'b, T> Clone for OwnedMatrixRow<'b, T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'b, T> VectorView<T> for OwnedMatrixRow<'b, T> {
    
    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

pub struct OwnedMatrixRowMut<'b, T> {
    data: &'b mut [T]
}

impl<'b, T> VectorView<T> for OwnedMatrixRowMut<'b, T> {
    
    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<'b, T> VectorViewMut<T> for OwnedMatrixRowMut<'b, T> {
    
    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }
}

impl<'b, T: 'b> LifetimeMatrixMutRowIter<'b, T> for MatrixOwned<T> {

    type RowRef = OwnedMatrixRowMut<'b, T>;
    
    fn rows(&'b mut self, i: usize) -> OwnedMatrixRowMut<'b, T> {
        OwnedMatrixRowMut {
            data: &mut self.data[i * self.cols..(i + 1) * self.cols]
        }
    }
}

impl<T: 'static> MutRowRefable<T> for MatrixOwned<T> {}

impl<'b, T: 'b> LifetimeMutRowsRefable<'b, T> for MatrixOwned<T> {
    
    fn get_rows(&'b mut self, fst: usize, snd: usize) -> (OwnedMatrixRowMut<'b, T>, OwnedMatrixRowMut<'b, T>) {
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
                OwnedMatrixRowMut { data: fst_row },
                OwnedMatrixRowMut {
                    data: &mut rest[snd_row_start..],
                },
            );
        } else {
            let part: &mut [T] = &mut self.data[(snd * cols)..((fst + 1) * cols)];
            let (snd_row, rest) = part.split_at_mut(cols);
            let fst_row_start = rest.len() - cols;
            return (
                OwnedMatrixRowMut {
                    data: &mut rest[fst_row_start..],
                },
                OwnedMatrixRowMut { data: snd_row },
            );
        }
    }
}

impl<T: 'static> MutRowsRefable<T> for MatrixOwned<T> {}











pub struct MatrixRow<'b, M, T> 
    where M: MatrixView<T>
{
    matrix: &'b M,
    row_index: usize,
    element: PhantomData<T>
}

pub struct MatrixRowMut<'b, M, T> 
    where M: MatrixViewMut<T>
{
    matrix: &'b mut M,
    row_index: usize,
    element: PhantomData<T>
}

impl<'a, M, T> Clone for MatrixRef<'a, M, T> 
    where M: MatrixView<T>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, M, T> Copy for MatrixRef<'a, M, T> 
    where M: MatrixView<T> {}

impl<'a, M, T> Clone for MatrixRow<'a, M, T> 
    where M: MatrixView<T>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, M, T> Copy for MatrixRow<'a, M, T> 
    where M: MatrixView<T> {}

impl<'a, M, T> VectorView<T> for MatrixRow<'a, M, T>
    where M: MatrixView<T> 
{
    fn len(&self) -> usize {
        self.matrix.cols()
    }

    fn at(&self, i: usize) -> &T {
        self.matrix.at(self.row_index, i)
    }
}

impl<'a, M, T> VectorView<T> for MatrixRowMut<'a, M, T>
    where M: MatrixViewMut<T> 
{
    fn len(&self) -> usize {
        self.matrix.col_count()
    }

    fn at(&self, i: usize) -> &T {
        self.matrix.at(self.row_index, i)
    }
}

impl<'a, M, T> VectorViewMut<T> for MatrixRowMut<'a, M, T>
    where M: MatrixViewMut<T> 
{
    fn at_mut(&mut self, i: usize) -> &mut T {
        self.matrix.at_mut(self.row_index, i)
    }
}