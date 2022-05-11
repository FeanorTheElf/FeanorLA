use super::super::ring::*;
use super::super::embedding::*;
use super::super::primitive::*;
use super::super::bigint::*;
use super::super::wrapper::*;
use std::cmp::Ordering;

pub trait IntegerRing: Ring + CanonicalIsomorphismInfo<StaticRing<i64>> + CanonicalIsomorphismInfo<BigIntRing> {

    fn to_float_approx(&self, el: &Self::El) -> f64;
    fn from_float_approx(&self, el: f64) -> Option<Self::El>;
    
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
        where F: FnMut(&Self::El) -> Self::El, Self: OrderedRing + EuclideanInfoRing
    {
        let mut begin = approx.clone();
        let mut step = self.one();
        while self.cmp(&f(&begin), &self.zero()) == Ordering::Greater {
            begin = self.sub_ref_snd(begin, &step);
            take_mut::take(&mut step, |x| self.add(x.clone(), x));
        }
        let mut end = approx;
        step = self.one();
        while self.cmp(&f(&end), &self.zero()) == Ordering::Less {
            end = self.add_ref(end, &step);
            take_mut::take(&mut step, |x| self.add(x.clone(), x));
        }
        return self.bisect(f, begin, end);
    }

    fn floor_div(&self, a: Self::El, b: &Self::El) -> Self::El 
        where Self: OrderedRing + EuclideanInfoRing
    {
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
        where F: FnMut(&Self::El) -> Self::El,
            Self: OrderedRing + EuclideanInfoRing
    {
        assert!(self.cmp(&f(&start), &self.zero()) != Ordering::Greater);
        assert!(self.cmp(&f(&end), &self.zero()) != Ordering::Less);
        if self.is_zero(&f(&end)) {
            return end;
        }
        let two = self.from_z(2);
        loop {
            let mid = self.floor_div(self.add_ref(start.clone(), &end), &two);
            if self.eq(&mid, &start) {
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

    fn abs_cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering
        where Self: OrderedRing
    {
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
    fn root_floor(&self, el: &Self::El, n: usize) -> Self::El
        where Self: OrderedRing + EuclideanInfoRing
    {
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

impl<'a, R: IntegerRing> IntegerRing for &'a R {

    fn to_float_approx(&self, el: &Self::El) -> f64 { (**self).to_float_approx(el) }
    fn from_float_approx(&self, el: f64) -> Option<Self::El>  { (**self).from_float_approx(el) }

}

impl IntegerRing for StaticRing<i64> {

    fn to_float_approx(&self, el: &Self::El) -> f64 {
        *el as f64
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El>  {
        Some(el as i64)
    }
}

impl IntegerRing for BigIntRing {

    fn to_float_approx(&self, el: &Self::El) -> f64 {
        el.to_float_approx()
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El>  {
        BigInt::from_float_approx(el)
    }
}

impl<R: IntegerRing> IntegerRing for WrappingRing<R> {
    
    fn to_float_approx(&self, el: &Self::El) -> f64 {
        self.wrapped_ring().to_float_approx(el.val())
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El>  {
        self.wrapped_ring().from_float_approx(el).map(|x| (self.wrapped_ring().clone()).bind_by_value(x))
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
fn test_root_floor() {
    let n = BigInt::from(7681).pow(32);
    assert_eq!(BigInt::from(7681), BigInt::RING.root_floor(&n, 32));
}