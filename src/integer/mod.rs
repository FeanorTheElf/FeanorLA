pub mod bigint;
pub mod primes;
pub mod roots;
pub mod bigint_soo;

pub use bigint::*;
use super::prelude::*;
use super::wrapper::*;

use std::cmp::Ordering;
use std::ops::Range;

pub trait IntegerRing: OrderedRing + EuclideanInfoRing + CanonicalIsomorphismInfo<StaticRing<i64>> + CanonicalIsomorphismInfo<BigIntRing> {

    fn to_float_approx(&self, el: &Self::El) -> f64;
    fn from_float_approx(&self, el: f64) -> Option<Self::El>;
    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self>;
    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self>;

    fn is_odd(&self, el: &Self::El) -> bool {
        !self.is_eq(&self.mul_pow_2(self.euclidean_div_pow_2(el.clone(), 1), 1), el)
    }

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 {
        assert!(!self.is_zero(el));
        roots::find_zero_floor(
            &i64::RING,
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

    fn floor_div(&self, a: Self::El, b: &Self::El) -> Self::El {
        let (q, r) = self.euclidean_div_rem(a.clone(), b);
        if self.cmp(&r, &self.zero()) == Ordering::Less {
            self.sub(q, self.one())
        } else {
            q
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

    fn abs(&self, el: Self::El) -> Self::El {
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
            return roots::find_zero_floor(
                self,
                |x| self.sub_ref_snd(self.mul(self.abs(x.clone()), self.pow(x, (n - 1) as u32)), el), 
                self.from_float_approx(root_approx).unwrap_or(self.zero())
            );
        } else {
            return roots::find_zero_floor(
                self,
                |x| self.sub_ref_snd(self.pow(x, n as u32), el), 
                self.from_float_approx(root_approx).unwrap_or(self.zero())
            );
        }
    }

    ///
    /// Generates a uniformly random number from the range 0 to end_exclusive, using
    /// entropy from the given rng.
    /// 
    /// Uses rejection sampling.
    /// 
    fn get_uniformly_random_oorandom(
        &self,
        rng: &mut oorandom::Rand32,
        end_exclusive: &El<Self>
    ) -> El<Self> {
        self.get_uniformly_random(|| rng.rand_u32(), end_exclusive)
    }

    ///
    /// Generates a uniformly random number from the range 0 to end_exclusive, using
    /// entropy from the given rng.
    /// 
    /// Uses rejection sampling.
    /// 
    fn get_uniformly_random<G>(
        &self,
        rng: G, 
        end_exclusive: &El<Self>
    ) -> El<Self> 
        where G: FnMut() -> u32
    {
        integer_ring_get_uniformly_random(self, rng, end_exclusive)
    }

    fn highest_dividing_power_of_two(&self, el: &El<Self>) -> usize {
        roots::bisect(
            &i64::RING,
            |k| {
                if self.is_eq(el, &self.mul_pow_2(self.euclidean_div_pow_2(el.clone(), *k as u64), *k as u64)) { -1 } else { 1 }
            },
            0,
            self.abs_log2_floor(el) as i64 + 1
        ) as usize
    }
}

fn integer_ring_get_uniformly_random<G, I>(
    ring: &I,
    mut rng: G, 
    end_exclusive: &El<I>
) -> El<I> 
    where I: IntegerRing, G: FnMut() -> u32
{
    assert!(ring.cmp(end_exclusive, &ring.zero()) == std::cmp::Ordering::Greater);
    let block_size = u32::BITS as u64;
    let k = ring.abs_log2_floor(&end_exclusive) + 1;
    loop {
        let mut i = 0;
        let mut current = ring.zero();
        while i + block_size < k {
            i += block_size;
            current = ring.add(ring.mul_pow_2(current, block_size), ring.from_z(rng() as i64));
        }
        let random_most_significant_bits = (rng() & ((1 << k % block_size) - 1)) as i64;
        current = ring.add(current, ring.mul_pow_2(ring.from_z(random_most_significant_bits), k - (k % block_size)));
        debug_assert_eq!(i + k % block_size, k);
        if ring.cmp(&current, end_exclusive) == std::cmp::Ordering::Less {
            return current;
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

impl<'a, R> IntegerRing for R 
    where R: RingDecorator, R::DecoratedRing: IntegerRing
{

    fn to_float_approx(&self, el: &Self::El) -> f64 { self.decorated_ring().to_float_approx(el) }
    fn from_float_approx(&self, el: f64) -> Option<Self::El>  { self.decorated_ring().from_float_approx(el) }
    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { self.decorated_ring().mul_pow_2(el, power) }
    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { self.decorated_ring().euclidean_div_pow_2(el, power) }
    fn abs_log2_floor(&self, el: &El<Self>) -> u64 { self.decorated_ring().abs_log2_floor(el) }
    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool { self.decorated_ring().abs_is_bit_set(el, bit) }
    fn floor_div(&self, a: Self::El, b: &Self::El) -> Self::El { self.decorated_ring().floor_div(a, b) }
    fn abs_cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering { self.decorated_ring().abs_cmp(lhs, rhs) }
    fn root_floor(&self, el: &Self::El, n: u64) -> Self::El { self.decorated_ring().root_floor(el, n) }
    fn get_uniformly_random_oorandom(&self, rng: &mut oorandom::Rand32, end_exclusive: &El<Self>) -> El<Self> { self.decorated_ring().get_uniformly_random_oorandom(rng, end_exclusive) }
    fn highest_dividing_power_of_two(&self, el: &El<Self>) -> usize { self.decorated_ring().highest_dividing_power_of_two(el) }
    fn abs(&self, el: Self::El) -> Self::El { self.decorated_ring().abs(el) }
    fn is_odd(&self, el: &Self::El) -> bool { self.decorated_ring().is_odd(el) }

    fn get_uniformly_random<G>(
        &self,
        rng: G, 
        end_exclusive: &El<Self>
    ) -> El<Self> 
        where G: FnMut() -> u32
    { self.decorated_ring().get_uniformly_random(rng, end_exclusive) }
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

impl<R: IntegerRing> IntegerRing for WrappingRing<R> {

    fn to_float_approx(&self, el: &Self::El) -> f64 { 
        self.wrapped_ring().to_float_approx(el.val())
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El>  { 
        self.wrapped_ring().from_float_approx(el).map(|x| self.from(x))
    }

    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { 
        self.from(self.wrapped_ring().mul_pow_2(el.into_val(), power))
    }

    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { 
        self.from(self.wrapped_ring().euclidean_div_pow_2(el.into_val(), power))
    }

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 { 
        self.wrapped_ring().abs_log2_floor(el.val())
    }

    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool { 
        self.wrapped_ring().abs_is_bit_set(el.val(), bit)
    }

    fn abs(&self, el: Self::El) -> Self::El { 
        self.from(self.wrapped_ring().abs(el.into_val()))
    }

    fn floor_div(&self, a: Self::El, b: &Self::El) -> Self::El {
        self.from(self.wrapped_ring().floor_div(a.into_val(), b.val()))
    }

    fn get_uniformly_random<G>(
        &self,
        rng: G, 
        end_exclusive: &El<Self>
    ) -> El<Self> 
        where G: FnMut() -> u32
    { self.from(self.wrapped_ring().get_uniformly_random(rng, end_exclusive.val())) }

    fn abs_cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        self.wrapped_ring().abs_cmp(lhs.val(), rhs.val())
    }

    fn root_floor(&self, el: &Self::El, n: u64) -> Self::El {
        self.from(self.wrapped_ring().root_floor(el.val(), n))
    }

    fn get_uniformly_random_oorandom(&self, rng: &mut oorandom::Rand32, end_exclusive: &El<Self>) -> El<Self> { 
        self.from(self.wrapped_ring().get_uniformly_random_oorandom(rng, end_exclusive.val()))
    }

    fn highest_dividing_power_of_two(&self, el: &El<Self>) -> usize { 
        self.wrapped_ring().highest_dividing_power_of_two(el.val())
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

    pub fn sqrt_floor(&self) -> Self {
        self.root_floor(2)
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
fn test_root_floor() {
    let n = BigInt::RING.pow(&BigInt::from(7681), 32);
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

#[cfg(test)]
use oorandom::Rand32;

#[test]
fn test_integer_ring_get_uniformly_random() {
    let ring = BigInt::RING;
    let end_exclusive = BigInt::from(18897856102); // more than 34 bits
    let mut rng = Rand32::new(0);
    let data: Vec<BigInt> = (0..1000).map(|_| integer_ring_get_uniformly_random(&ring, || rng.rand_u32(), &end_exclusive)).collect();
    let limit = BigInt::RING.mul_pow_2(BigInt::RING.one(), 34);
    assert!(data.iter().any(|x| BigInt::RING.is_geq(x, &limit))); // failure probability is less than 10^-37
    let bit_average = data.iter().map(|x| if ring.abs_is_bit_set(x, 33) { 0. } else { 1. }).sum::<f64>() / 1000.;
    assert!(0.3 < bit_average && bit_average < 0.7);  // failure probability is less than 10^-34
}

#[test]
fn test_i64_euclidean_div_pow_2() {
    assert_eq!(0, i64::RING.euclidean_div_pow_2(-1, 1))
}