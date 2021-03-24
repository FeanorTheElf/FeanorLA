use super::alg::*;

use std::cell::RefCell;
use std::ops::{Add, Mul, Sub, Div, Rem, Neg};

pub struct RingReferencingEl<'a, R>
    where R: Ring
{
    ring: &'a R,
    storage: &'a RefCell<Vec<R::El>>,
    index: usize
}

impl<'a, R> Copy for RingReferencingEl<'a, R>
where R: Ring
{}

impl<'a, R> Clone for RingReferencingEl<'a, R>
where R: Ring
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, R> RingReferencingEl<'a, R>
where R: Ring
{
    pub fn create(ring: &'a R, storage: &'a RefCell<Vec<R::El>>, val: R::El) -> RingReferencingEl<'a, R> {
        let res_index = {
            let mut vec = storage.borrow_mut();
            let res_index = vec.len();
            vec.push(val);
            res_index
        };
        return RingReferencingEl {
            ring: ring,
            storage: storage,
            index: res_index
        };
    }

    pub fn value(&self) -> R::El {
        let vec = self.storage.borrow();
        return vec[self.index].clone();
    }
}

impl<'a, R> Add<RingReferencingEl<'a, R>> for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn add(self, rhs: Self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            self.ring.add_ref(vec[self.index].clone(), &vec[rhs.index])
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> Mul<RingReferencingEl<'a, R>> for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn mul(self, rhs: Self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            self.ring.mul_ref(&vec[self.index], &vec[rhs.index])
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> Sub<RingReferencingEl<'a, R>> for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn sub(self, rhs: Self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            self.ring.sub_ref_snd(vec[self.index].clone(), &vec[rhs.index])
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> Neg for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn neg(self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            self.ring.neg(vec[self.index].clone())
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> Div for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn div(self, rhs: Self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            if self.ring.is_euclidean() {
                self.ring.euclidean_div(vec[self.index].clone(), &vec[rhs.index])
            } else {
                self.ring.div(vec[self.index].clone(), &vec[rhs.index])
            }
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> Rem for RingReferencingEl<'a, R>
where R: Ring
{
    type Output = RingReferencingEl<'a, R>;

    fn rem(self, rhs: Self) -> Self::Output {
        let res_el = {
            let vec = self.storage.borrow();
            self.ring.euclidean_rem(vec[self.index].clone(), &vec[rhs.index])
        };
        return RingReferencingEl::create(self.ring, self.storage, res_el);
    }
}

impl<'a, R> std::fmt::Display for RingReferencingEl<'a, R>
where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.ring.format(&self.storage.borrow()[self.index], f, false)
    }
}

impl<'a, R> PartialEq for RingReferencingEl<'a, R>
where R: Ring
{
    fn eq(&self, rhs: &Self) -> bool {
        self.ring.eq(&self.storage.borrow()[self.index], &self.storage.borrow()[rhs.index])
    }
}

impl<'a, R> Eq for RingReferencingEl<'a, R>
where R: Ring
{}

pub trait PostprocessOutput<T> {

    type Output;

    fn postprocess(self, vec: &Vec<T>) -> Self::Output;
}

impl<'a, R> PostprocessOutput<R::El> for RingReferencingEl<'a, R> 
    where R: Ring
{
    type Output = R::El;

    fn postprocess(self, vec: &Vec<R::El>) -> Self::Output {
        return vec[self.index].clone()
    }
}

macro_rules! fixed_ring_env {
    ($ring:expr; $($var:ident),*; $content:tt) => {
        {
            let the_ring: &_ = $ring;
            let storage = std::cell::RefCell::from(Vec::new());
            let bind = |x| RingReferencingEl::create(the_ring, &storage, x);
            $(
                let $var = bind(($var).clone());
            )*
            let result = {
                $content
            };
            let storage_vec = storage.borrow();
            result.postprocess(&*storage_vec)
        }
    };
}

#[test]
fn test_exper() {
    let a = 4;
    let b = 10;
    let d = fixed_ring_env!{ &StaticRing::<i32>::RING; a, b; {
        let mut c = a + (b * a) + a + a;
        c = c / b;
        c = c * a;
        c
    }};
    assert_eq!(20, d);
}