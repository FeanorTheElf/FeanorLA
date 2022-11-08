use std::marker::PhantomData;
use std::ops::{Add, Sub, BitXor, RangeBounds, Bound};
use std::cmp::PartialOrd;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Bitset64 {
    bits: u64
}

fn smallest_k_elements_between(superset: Bitset64, k: usize, lower: usize, upper: usize) -> Bitset64 {
    debug_assert!(k <= superset.len());
    debug_assert!(superset.is_subset(Bitset64::from(lower..upper)));
    if k == 0 {
        return Bitset64::new();
    } else if k == 1 {
        return Bitset64::singleton(superset.smallest_element().unwrap());
    }
    let middle = (lower + upper) / 2;
    let lower_half = superset.intersect_range(..middle);
    let upper_half = superset.intersect_range(middle..);
    let size = lower_half.len();
    if size < k {
        return smallest_k_elements_between(upper_half, k - size, middle, upper) + lower_half;
    } else {
        return smallest_k_elements_between(lower_half, k, lower, middle);
    }
}

impl Bitset64 {

    pub fn new() -> Self {
        Bitset64 {
            bits: 0
        }
    }

    pub fn singleton(a: usize) -> Self {
        let mut result = Bitset64::new();
        result.insert(a);
        return result;
    }

    pub fn is_subset(&self, rhs: Bitset64) -> bool {
        (*self - rhs).empty()
    }

    pub fn intersect(&self, rhs: Bitset64) -> Bitset64 {
        Bitset64 {
            bits: self.bits & rhs.bits
        }
    }

    pub fn intersect_range<R: RangeBounds<usize>>(&self, r: R) -> Bitset64 {
        self.intersect(Bitset64::from(r))
    }

    pub fn union(&self, rhs: Bitset64) -> Bitset64 {
        Bitset64 {
            bits: self.bits | rhs.bits
        }
    }

    pub fn any_element(&self) -> Option<usize> {
        self.smallest_element()
    }

    pub fn symm_diff(&self, rhs: Bitset64) -> Bitset64 {
        Bitset64 {
            bits: self.bits ^ rhs.bits
        }
    }

    pub fn diff(&self, rhs: Bitset64) -> Bitset64 {
        Bitset64 {
            bits: self.bits & !rhs.bits
        }
    }

    pub fn contains(&self, v: usize) -> bool {
        assert!(v < 64);
        ((self.bits >> v) & 1) != 0
    }

    pub fn insert(&mut self, v: usize) {
        assert!(v < 64);
        self.bits |= 1 << v;
    }

    pub fn delete(&mut self, v: usize) {
        assert!(v < 64);
        self.bits &= !(1 << v);
    }

    pub fn toggle(&mut self, v: usize) {
        assert!(v < 64);
        self.bits ^= 1 << v;
    }

    pub fn empty(&self) -> bool {
        self.bits == 0
    }

    pub fn replace(&mut self, old: usize, new: usize) {
        self.delete(old);
        self.insert(new);
    }

    pub fn greatest_element(&self) -> Option<usize> {
        63u32.checked_sub(self.bits.leading_zeros()).map(|x| x as usize)
    }

    pub fn smallest_element(&self) -> Option<usize> {
        let v = self.bits.trailing_zeros() as usize;
        if v < 64 { Some(v) } else { None }
    }

    pub fn smallest_element_geq(&self, bound: usize) -> Option<usize> {
        let v = if bound >= 64 { 64 } else { (self.bits >> bound).trailing_zeros() as usize + bound };
        if v < 64 { Some(v) } else { None }
    }

    pub fn len(&self) -> usize {
        self.bits.count_ones() as usize
    }

    pub fn iter(&self) -> Bitset64Iter {
        Bitset64Iter {
            _phantom: PhantomData,
            bitset: *self
        }
    }

    pub fn smallest_elements(&self, k: usize) -> Self {
        smallest_k_elements_between(*self, k, 0, 64)
    }
}

impl Add<Bitset64> for Bitset64 {
    type Output = Bitset64;

    fn add(self, rhs: Bitset64) -> Self::Output {
        self.union(rhs)
    }
}

impl Sub<Bitset64> for Bitset64 {
    type Output = Bitset64;

    fn sub(self, rhs: Bitset64) -> Self::Output {
        self.diff(rhs)
    }
}

impl BitXor<Bitset64> for Bitset64 {
    type Output = Bitset64;

    fn bitxor(self, rhs: Bitset64) -> Self::Output {
        self.symm_diff(rhs)
    }
}

impl PartialOrd<Bitset64> for Bitset64 {

    fn partial_cmp(&self, rhs: &Bitset64) -> Option<std::cmp::Ordering> {
        if self == rhs {
            Some(std::cmp::Ordering::Equal)
        } else if self.is_subset(*rhs) {
            Some(std::cmp::Ordering::Less)
        } else if rhs.is_subset(*self) {
            Some(std::cmp::Ordering::Greater)
        } else {
            None
        }
    }
}

impl std::fmt::Debug for Bitset64 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "bitset {{")?;
        for el in self.iter().take(1) {
            write!(f, "{}", el)?;

        }
        for el in self.iter().skip(1) {
            write!(f, ", {}", el)?;
        }
        return write!(f, "}}");
    }
}

#[derive(Debug, Clone)]
pub struct Bitset64Iter<'a> {
    bitset: Bitset64,
    _phantom: PhantomData<&'a Bitset64>
}

impl<'a> Iterator for Bitset64Iter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.bitset.smallest_element();
        if let Some(v) = result {
            self.bitset.delete(v);
        }
        return result;
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.bitset.len(), Some(self.bitset.len()))
    }
}

impl<'a> DoubleEndedIterator for Bitset64Iter<'a> {
    
    fn next_back(&mut self) -> Option<Self::Item> {
        let result = self.bitset.greatest_element();
        if let Some(v) = result {
            self.bitset.delete(v);
        }
        return result;
    }
}

impl<'a> ExactSizeIterator for Bitset64Iter<'a> {}
impl<'a> std::iter::FusedIterator for Bitset64Iter<'a> {}

impl<R: RangeBounds<usize>> From<R> for Bitset64 {

    fn from(r: R) -> Self {
        let lower = match r.start_bound() {
            Bound::Included(a) => *a,
            Bound::Excluded(a) => *a + 1,
            Bound::Unbounded => 0
        };
        assert!(lower <= 63);
        let upper = match r.end_bound() {
            Bound::Included(a) => *a + 1,
            Bound::Excluded(a) => *a,
            Bound::Unbounded => 64
        };
        assert!(upper <= 64);
        if lower >= upper {
            Bitset64::new()
        } else {
            Bitset64 {
                bits: (u64::MAX >> (64 + lower - upper)) << lower
            }
        }
    }
}

impl std::iter::FromIterator<usize> for Bitset64 {
    
    fn from_iter<I: IntoIterator<Item = usize>>(iter: I) -> Self {
        let mut result = Bitset64::new();
        for c in iter {
            result.insert(c);
        }
        return result;
    }
}

impl IntoIterator for Bitset64 {
    type Item = usize;
    type IntoIter = Bitset64Iter<'static>;

    fn into_iter(self) -> Self::IntoIter {
        Bitset64Iter {
            _phantom: PhantomData,
            bitset: self
        }
    }
}

#[derive(Debug, Clone)]
pub struct RangePowerset {
    end: u64,
    current: Bitset64
}

impl RangePowerset {

    pub fn peek(&self) -> Option<Bitset64> {
        if self.current.bits == self.end {
            return None;
        } else {
            return Some(self.current);
        }
    }
}

impl Iterator for RangePowerset {
    type Item = Bitset64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current.bits == self.end {
            return None;
        } else {
            let result = self.current;
            self.current.bits += 1;
            return Some(result);
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = (self.end - self.current.bits) as usize;
        (size, Some(size))
    }
}
impl ExactSizeIterator for RangePowerset {}
impl std::iter::FusedIterator for RangePowerset {}

pub fn range_powerset(end_exclusive: usize) -> RangePowerset {
    assert!(end_exclusive < 64);
    RangePowerset {
        end: 1 << end_exclusive,
        current: Bitset64::new()
    }
}

#[allow(unused)]
macro_rules! bitset {
    ($($element:expr),*) => {
        {
            #[allow(unused_mut)]
            let mut result_set = Bitset64::new();
            $(result_set.insert($element);)*
            result_set
        }
    };
}

#[test]
fn test_from_range() {
    assert_eq!(bitset!{2, 3, 4}, Bitset64::from(2..5));
    assert_eq!(bitset!{3, 4, 5, 6}, Bitset64::from(3..=6));
    assert_eq!(bitset!{0, 1, 2}, Bitset64::from(..3));
    assert_eq!(bitset!{63}, Bitset64::from(63..));
}

#[test]
fn test_operations() {
    let a = bitset!{1, 2, 3, 4};
    let b = bitset!{1, 3, 5, 6};
    assert_eq!(bitset!{1, 2, 3, 4, 5, 6}, a + b);
    assert_eq!(bitset!{1, 3}, a.intersect(b));
    assert_eq!(bitset!{2, 4, 5, 6}, a ^ b);
    assert_eq!(bitset!{2, 4}, a - b);
}

#[test]
fn test_subset() {
    let a = Bitset64::from(0..5);
    let b = Bitset64::from(1..4);
    assert!(b.is_subset(a));
    assert!(b < a);
    assert!(!(a < a));
    assert!(a <= a);
    assert!(!(b >= a));
}

#[test]
fn test_smallest_elements() {
    let a = bitset!{0, 2, 3, 4, 6, 9, 15};
    assert_eq!(Some(0), a.smallest_element());
    assert_eq!(Some(6), a.smallest_element_geq(5));
    assert_eq!(bitset!{0, 2, 3}, a.smallest_elements(3));
    assert_eq!(bitset!{0, 2, 3, 4, 6, 9}, a.smallest_elements(6));
    assert_eq!(bitset!{}, a.smallest_elements(0));
}

#[test]
fn test_range_powerset() {
    let mut iter = range_powerset(5);
    assert_eq!(32, iter.clone().count());
    assert_eq!(Some(bitset!{}), iter.next());
    assert_eq!(Some(bitset!{0}), iter.next());
    assert_eq!(30, iter.count());
}