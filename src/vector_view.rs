
pub trait VectorView<T> {
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;

    fn assert_in_range(&self, index: usize) {
        assert!(index < self.len(), "Vector index {} out of range 0..{}", index, self.len());
    }
}

pub trait VectorViewMut<T>: VectorView<T> {
    fn at_mut(&mut self, index: usize) -> &mut T;
    fn swap(&mut self, i: usize, j: usize);
}
