
pub trait Indexed<'a, T> {
    type Output;

    fn get(&'a self, index: T) -> Self::Output;
}

pub trait IndexedMut<'a, T> {
    type Output;

    fn get_mut(&'a mut self, index: T) -> Self::Output;
}
