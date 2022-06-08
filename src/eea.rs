use super::ring::*;

use std::mem::swap;
use std::cmp::Ordering;

///
/// For a, b computes s, t, d such that `s*a + t*b == d` is a greatest 
/// common divisor of a and b. d is only unique up to units, and s, t 
/// are not unique at all. No guarantees are given on which
/// of these solutions is returned. For integers, see signed_eea 
/// which gives more guarantees.
/// 
/// The given ring must be euclidean
/// 
pub fn eea<R>(ring: &R, fst: R::El, snd: R::El) -> (R::El, R::El, R::El) 
    where R::El: Clone, R: EuclideanInfoRing
{
    assert!(ring.is_euclidean().can_use());

    let (mut a, mut b) = (fst, snd);
    let (mut sa, mut ta) = (ring.one(), ring.zero());
    let (mut sb, mut tb) = (ring.zero(), ring.one());

    // invariant: sa * a + ta * b = fst, sb * a + tb * b = snd
    while !ring.is_eq(&b, &ring.zero()) {
        let (quot, rem) = ring.euclidean_div_rem(a, &b);
        ta = ring.sub(ta, ring.mul_ref(&quot, &tb));
        sa = ring.sub(sa, ring.mul_ref(&quot, &sb));
        a = rem;
        swap(&mut a, &mut b);
        swap(&mut sa, &mut sb);
        swap(&mut ta, &mut tb);
    }
    return (sa, ta, a);
}

/// 
/// For integers a, b finds the smallest integers s, t so that 
/// `s*a + t*b == gcd(a, b)` is the greatest common divisor of a, b.
/// 
/// Details: s and t are not unique, this function will return 
/// the smallest tuple (s, t) (ordered by the total ordering
/// given by `(s, t) ≤ (u, v) :<=> |s| ≤ |u| and |t| ≤ |v|`). 
/// In the case |a| = |b|, there are two minimal elements, in this case, it is 
/// unspecified whether this function returns (±1, 0, a) or (0, ±1, a). 
/// We define the greatest common divisor gcd(a, b) as the minimal
/// element of the set of integers dividing a and b (ordered by divisibility), 
/// whose sign matches the sign of a.
/// 
/// In particular, have 
/// ```
/// # use feanor_la::eea::signed_gcd;
/// # use feanor_la::primitive::RingEl;
/// assert_eq!(2, signed_gcd(6, 8, &i64::RING));
/// assert_eq!(0, signed_gcd(0, 0, &i64::RING)); 
/// assert_eq!(5, signed_gcd(0, -5, &i64::RING));
/// assert_eq!(-5, signed_gcd(-5, 0, &i64::RING)); 
/// assert_eq!(-1, signed_gcd(-1, 1, &i64::RING));
/// assert_eq!(1, signed_gcd(1, -1, &i64::RING));
/// ```
/// and therefore `signed_eea(6, 8) == (-1, 1, 2)`, 
/// `signed_eea(-6, 8) == (-1, -1, -2)`, 
/// `signed_eea(8, -6) == (1, 1, 2)`, 
/// `signed_eea(0, 0) == (0, 0, 0)`
/// 
pub fn signed_eea<R>(fst: R::El, snd: R::El, ring: &R) -> (R::El, R::El, R::El)
    where R: OrderedRing + EuclideanInfoRing
{
    if ring.is_zero(&fst) {
        return match ring.cmp(&snd, &ring.zero()) {
            Ordering::Equal => (ring.zero(), ring.zero(), ring.zero()),
            Ordering::Less => (ring.zero(), ring.neg(ring.one()), ring.neg(snd)),
            Ordering::Greater => (ring.zero(), ring.one(), snd)
        };
    }
    let fst_negative = ring.cmp(&fst, &ring.zero());

    let (s, t, d) = eea(ring, fst, snd);
    
    // the sign is not consistent (potentially toggled each iteration), 
    // so normalize here
    if ring.cmp(&d, &ring.zero()) == fst_negative {
        return (s, t, d);
    } else {
        return (ring.neg(s), ring.neg(t), ring.neg(d));
    }
}

#[cfg(test)]
use super::primitive::*;

#[test]
fn test_eea_sign() {
    assert_eq!((2, -1, 1), signed_eea(3, 5, &i64::RING));
    assert_eq!((-1, 2, 1), signed_eea(5, 3, &i64::RING));
    assert_eq!((2, 1, -1), signed_eea(-3, 5, &i64::RING));
    assert_eq!((-1, -2, 1), signed_eea(5, -3, &i64::RING));
    assert_eq!((2, 1, 1), signed_eea(3, -5, &i64::RING));
    assert_eq!((-1, -2, -1), signed_eea(-5, 3, &i64::RING));
    assert_eq!((2, -1, -1), signed_eea(-3, -5, &i64::RING));
    assert_eq!((-1, 2, -1), signed_eea(-5, -3, &i64::RING));
    assert_eq!((0, 0, 0), signed_eea(0, 0, &i64::RING));
    assert_eq!((1, 0, 4), signed_eea(4, 0, &i64::RING));
    assert_eq!((0, 1, 4), signed_eea(0, 4, &i64::RING));
    assert_eq!((1, 0, -4), signed_eea(-4, 0, &i64::RING));
    assert_eq!((0, -1, 4), signed_eea(0, -4, &i64::RING));
}

#[test]
fn test_signed_eea() {
    assert_eq!((-1, 1, 2), signed_eea(6, 8, &i64::RING));
    assert_eq!((2, -1, 5), signed_eea(15, 25, &i64::RING));
    assert_eq!((4, -7, 2), signed_eea(32, 18, &i64::RING));
}

/// 
/// Finds a greatest common divisor of a and b.
/// 
/// The gcd of two elements in a euclidean ring is the (w.r.t divisibility) greatest
/// element that divides both elements. It is unique up to multiplication with units. 
/// This function makes no guarantees on which of these will be returned.
/// 
/// If this is required, see also signed_gcd that gives precise statement on the
/// sign of the gcd in case of two integers.
/// 
/// The given ring must be euclidean
/// 
pub fn gcd<R>(ring: &R, a: R::El, b: R::El) -> R::El 
    where R: EuclideanInfoRing, R::El: Clone
{
    assert!(ring.is_euclidean().can_use());
    
    let (_, _, d) = eea(ring, a, b);
    return d;
}

/// 
/// Finds the greatest common divisor of a and b.
/// 
/// The gcd is only unique up to multiplication by units, so in this case up to sign.
/// However, this function guarantees the following behavior w.r.t different signs:
/// 
/// a < 0 => gcd(a, b) < 0
/// a > 0 => gcd(a, b) > 0
/// sign of b is irrelevant
/// gcd(0, 0) = 0
/// 
pub fn signed_gcd<R>(a: R::El, b: R::El, ring: &R) -> R::El
    where R: EuclideanInfoRing + OrderedRing
{
    let (_, _, d) = signed_eea(a, b, ring);
    return d;
}

pub fn lcm<R>(ring: &R, fst: R::El, snd: R::El) -> R::El
    where R: EuclideanInfoRing
{
    ring.euclidean_div(ring.mul_ref(&fst, &snd), &gcd(ring, fst, snd))
}

#[test]
fn test_gcd() {
    assert_eq!(3, signed_gcd(15, 6, &i64::RING));
    assert_eq!(3, signed_gcd(6, 15, &i64::RING));

    assert_eq!(7, signed_gcd(0, 7, &i64::RING));
    assert_eq!(7, signed_gcd(7, 0, &i64::RING));
    assert_eq!(0, signed_gcd(0, 0, &i64::RING));

    assert_eq!(1, signed_gcd(9, 1, &i64::RING));
    assert_eq!(1, signed_gcd(1, 9, &i64::RING));

    assert_eq!(1, signed_gcd(13, 300, &i64::RING));
    assert_eq!(1, signed_gcd(300, 13, &i64::RING));

    assert_eq!(-3, signed_gcd(-15, 6, &i64::RING));
    assert_eq!(3, signed_gcd(6, -15, &i64::RING));
    assert_eq!(-3, signed_gcd(-6, -15, &i64::RING));
}
