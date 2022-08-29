pub mod bigint;

pub use bigint::*;
use super::ring::*;
use super::embedding::*;
use super::primitive::*;
use super::wrapper::*;

use std::cmp::Ordering;
use std::ops::Range;

pub trait IntegerRing: OrderedRing + EuclideanInfoRing + CanonicalIsomorphismInfo<StaticRing<i64>> + CanonicalIsomorphismInfo<BigIntRing> {

    fn to_float_approx(&self, el: &Self::El) -> f64;
    fn from_float_approx(&self, el: f64) -> Option<Self::El>;
    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self>;
    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self>;

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 {
        assert!(!self.is_zero(el));
        i64::RING.find_zero_floor(
            |power: &i64| if *power >= 0 && self.is_zero(&self.euclidean_div_pow_2(el.clone(), *power as u64)) {
                1
            } else {
                -1
            },
            0
        ) as u64
    }

    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool {
        !self.is_zero(&self.euclidean_rem(self.euclidean_div_pow_2(el.clone(), bit), &self.from_z(2)))
    }
    
    ///
    /// Given an increasing, continuous function f: R -> R that is negative for some x1 and 
    /// positive for some x2, finds the floor of some root of f (if f is strictly increasing, 
    /// this is unique).
    /// 
    /// # General case
    /// 
    /// This function also works in a slightly more general context. Assume that
    /// f(x) is negative for all sufficiently small x and positive for all suffiently 
    /// large x. Then this function will return the floor of some root of f. Note that
    /// this might not be a root of f, even if f has integral roots.
    /// 
    /// # Complexity
    /// 
    /// This function runs in O((T + log(d)) * log(d)) where d is the error made in 
    /// approx (i.e. the difference between the found root x and approx) and T is the
    /// time required for computing f on a value between x - d and x + d.
    /// 
    fn find_zero_floor<F>(&self, mut f: F, approx: Self::El) -> Self::El
        where F: FnMut(&Self::El) -> Self::El
    {
        let mut begin = approx.clone();
        let mut step = self.one();
        let two = self.from_z(2);
        while self.cmp(&f(&begin), &self.zero()) == Ordering::Greater {
            begin = self.sub_ref_snd(begin, &step);
            self.mul_assign(&mut step, two.clone());
        }
        let mut end = approx;
        step = self.one();
        while self.cmp(&f(&end), &self.zero()) == Ordering::Less {
            end = self.add_ref(end, &step);
            self.mul_assign(&mut step, two.clone());
        }
        return self.bisect(f, begin, end);
    }

    fn floor_div(&self, a: Self::El, b: &Self::El) -> Self::El {
        let (q, r) = self.euclidean_div_rem(a.clone(), b);
        if self.cmp(&r, &self.zero()) == Ordering::Less {
            self.sub(q, self.one())
        } else {
            q
        }
    }

    ///
    /// Given a continuous function f: R -> R that is negative on `begin` and 
    /// positive on `end`, finds the floor of some root of f. Note that even
    /// if f has integral roots, the returned value does not have to be a root
    /// of f.
    /// 
    /// # Complexity
    /// 
    /// This function runs in O((T + log(d)) * log(d)) where d is the difference between
    /// begin and end and T is the time required for computing f on a value between 
    /// begin and end. 
    /// 
    fn bisect<F>(&self, mut f: F, mut start: Self::El, mut end: Self::El) -> Self::El
        where F: FnMut(&Self::El) -> Self::El
    {
        assert!(self.cmp(&f(&start), &self.zero()) != Ordering::Greater);
        assert!(self.cmp(&f(&end), &self.zero()) != Ordering::Less);
        if self.is_zero(&f(&end)) {
            return end;
        }
        let two = self.from_z(2);
        loop {
            let mid = self.floor_div(self.add_ref(start.clone(), &end), &two);
            if self.is_eq(&mid, &start) {
                return start;
            }
            match self.cmp(&f(&mid), &self.zero()) {
                Ordering::Less => {
                    start = mid;
                },
                Ordering::Greater => {
                    end = mid;
                },
                _ => {
                    return mid;
                }
            }
        }
    }

    fn abs_cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        let zero = self.zero();
        match (self.cmp(&lhs, &zero), self.cmp(&rhs, &zero)) {
            (Ordering::Equal, Ordering::Equal) => Ordering::Equal,
            (Ordering::Less, Ordering::Equal) => Ordering::Greater,
            (Ordering::Greater, Ordering::Equal) => Ordering::Greater,
            (Ordering::Equal, Ordering::Less) => Ordering::Less,
            (Ordering::Equal, Ordering::Greater) => Ordering::Less,
            (Ordering::Less, Ordering::Less) => self.cmp(rhs, lhs),
            (Ordering::Greater, Ordering::Greater) => self.cmp(lhs, rhs),
            (Ordering::Less, Ordering::Greater) => {
                let neg_lhs = self.neg(lhs.clone());
                self.cmp(&neg_lhs, rhs)
            },
            (Ordering::Greater, Ordering::Less) => {
                let neg_rhs = self.neg(rhs.clone());
                self.cmp(lhs, &neg_rhs)
            }
        }
    }

    fn abs(&self, el: Self::El) -> Self::El
        where Self: OrderedRing
    {
        if self.cmp(&el, &self.zero()) == Ordering::Less {
            self.neg(el)
        } else {
            el
        }
    }

    ///
    /// Computes the n-th root of this number.
    /// 
    /// # Complexity
    /// 
    /// The asymptotic worst-case complexity is O(log(n)^2), however it
    /// will be quite fast on most inputs due to internal use of floating
    /// point approximations.
    /// 
    fn root_floor(&self, el: &Self::El, n: u64) -> Self::El {
        assert!(n > 0);
        let root_approx = self.to_float_approx(el).powf(1. / n as f64);
        if n % 2 == 0 {
            return self.find_zero_floor(
                |x| self.sub_ref_snd(self.mul(self.abs(x.clone()), self.pow(x, (n - 1) as u32)), el), 
                self.from_float_approx(root_approx).unwrap_or(self.zero())
            );
        } else {
            return self.find_zero_floor(
                |x| self.sub_ref_snd(self.pow(x, n as u32), el), 
                self.from_float_approx(root_approx).unwrap_or(self.zero())
            );
        }
    }
}

pub fn range_iter<I: IntegerRing>(range: Range<El<I>>, ring: I) -> impl Iterator<Item = El<I>> {
    let end_minus_2 = ring.sub(range.end, ring.from_z(2));
    std::iter::repeat(()).scan(ring.sub(range.start, ring.one()), move |state, ()| {
        if ring.cmp(&*state, &end_minus_2) != Ordering::Greater {
            ring.add_assign(state, ring.one());
            return Some(state.clone());
        } else {
            return None;
        }
    })
}

impl<'a, R: IntegerRing> IntegerRing for &'a R {

    fn to_float_approx(&self, el: &Self::El) -> f64 { (**self).to_float_approx(el) }
    fn from_float_approx(&self, el: f64) -> Option<Self::El>  { (**self).from_float_approx(el) }
    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { (**self).mul_pow_2(el, power) }
    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { (**self).euclidean_div_pow_2(el, power) }
}

impl IntegerRing for StaticRing<i64> {

    fn to_float_approx(&self, el: &Self::El) -> f64 {
        *el as f64
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El>  {
        Some(el as i64)
    }

    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { 
        el << power
    }

    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> {
        el / (1 << power)
    }

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 {
        assert!(!self.is_zero(el));
        (i64::BITS - el.abs().leading_zeros() - 1) as u64
    }

    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool {
        (el.abs() >> bit) & 1 == 1
    }
}

impl<R: IntegerRing> RingElWrapper<R> 
{
    pub fn to_float_approx(&self) -> f64 {
        self.parent_ring().to_float_approx(self.val())
    }

    pub fn root_floor(&self, n: u64) -> Self {
        let result = self.parent_ring().root_floor(self.val(), n);
        return RingElWrapper::new(result, self.parent_ring().clone());
    }

    pub fn floor_div(self, rhs: &RingElWrapper<R>) -> Self {
        let (el, ring) = self.destruct();
        let result = ring.floor_div(el, rhs.val());
        return RingElWrapper::new(result, ring);
    }

    pub fn euclidean_div_pow_2(self, exponent: u64) -> Self {
        let (el, ring) = self.destruct();
        let result = ring.euclidean_div_pow_2(el, exponent);
        return RingElWrapper::new(result, ring);
    }
}

#[test]
fn test_find_zero_floor() {
    let f = |x: &BigInt| BigInt::RING.mul_ref(x, x) - 234867;
    assert_eq!(BigInt::from(484), BigInt::RING.find_zero_floor(f, BigInt::ZERO));

    let f = |x: &BigInt| x.clone();
    assert_eq!(BigInt::ZERO, BigInt::RING.find_zero_floor(f, BigInt::ZERO));
}

#[test]
fn test_find_zero_floor_i64() {
    let f = |x: &i64| x.pow(3) - *x + 478;
    assert_eq!(-8, i64::RING.find_zero_floor(f, 0));
}

#[test]
fn test_root_floor() {
    let n = BigInt::from(7681).pow(32);
    assert_eq!(BigInt::from(7681), BigInt::RING.root_floor(&n, 32));
}

#[test]
fn test_range_iter() {
    let i = |x| BigInt::from(x);
    assert_eq!(
        vec![i(1), i(2)],
        range_iter(i(1)..i(3), BigInt::RING).collect::<Vec<_>>()
    );
    assert_eq!(
        Vec::<BigInt>::new(),
        range_iter(i(1)..i(-3), BigInt::RING).collect::<Vec<_>>()
    );
}