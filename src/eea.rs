use super::alg::*;

use std::mem::swap;

///
/// For a, b computes s, t, d such that s*a + t*b = d is a greatest common divisor of a and b.
/// d is only unique up to units, and s, t are not unique at all. No guarantees are given on which
/// of these solutions is returned. For integers, see signed_eea which gives more guarantees.
/// 
pub fn eea<R: EuclideanRing>(ring: &R, fst: R::El, snd: R::El) -> (R::El, R::El, R::El) 
    where R::El: Clone
{
    let (mut a, mut b) = (fst, snd);
    let (mut sa, mut ta) = (ring.one(), ring.zero());
    let (mut sb, mut tb) = (ring.zero(), ring.one());

    // invariant: sa * a + ta * b = fst, sb * a + tb * b = snd
    while !ring.eq(&b, &ring.zero()) {
        let (quot, rem) = ring.div_rem(a, b.clone());
        ta = ring.sub(ta, ring.mul(quot.clone(), tb.clone()));
        sa = ring.sub(sa, ring.mul(quot.clone(), sb.clone()));
        a = rem;
        swap(&mut a, &mut b);
        swap(&mut sa, &mut sb);
        swap(&mut ta, &mut tb);
    }
    return (sa, ta, a);
}

/// 
/// For integers a, b finds the smallest integers s, t so that `s*a + t*b = gcd(a, b)` is the greatest common divisor of a, b.
/// 
/// Details: s and t are not unique, this function will return the smallest tuple (s, t) (ordered by the total ordering
/// given by `(s, t) ≤ (u, v) :<=> |s| ≤ |u| and |t| ≤ |v|`). In the case |a| = |b|, there are two minimal elements, in this case, it is 
/// unspecified whether this function returns (±1, 0, a) or (0, ±1, a). We define the greatest common divisor gcd(a, b) as the minimal
/// element of the set of integers dividing a and b (ordered by divisibility), whose sign matches the sign of a.
/// 
/// In particular, have `gcd(6, 8) == 2`, `gcd(0, 0) == 0`, `gcd(0, x) == |x|`, `gcd(x, 0) == x`, `gcd(-1, 1) == -1`, `gcd(1, -1) == 1`
/// and therefore `signed_eea(6, 8) == (-1, 1, 2)`, `signed_eea(-6, 8) == (-1, -1, -2)`, `signed_eea(8, -6) == (1, 1, 2)`, `signed_eea(0, 0) == (0, 0, 0)`
/// 
pub fn signed_eea<Int: Integer + Copy>(fst: Int, snd: Int) -> (Int, Int, Int) {
    if fst == Int::zero() {
        if snd == Int::zero() {
            return (Int::zero(), Int::zero(), Int::zero());
        } else if snd < Int::zero() {
            return (Int::zero(), -Int::one(), -snd);
        } else {
            return (Int::zero(), Int::one(), snd);
        }
    }

    let (s, t, d) = eea(&StaticRing::<Int>::RING, fst, snd);
    
    // the sign is not consistent (potentially toggled each iteration), so normalize here
    if (d < Int::zero()) == (fst < Int::zero()) {
        return (s, t, d);
    } else {
        return (-s, -t, -d);
    }
}

#[test]
fn test_eea_sign() {
    assert_eq!((2, -1, 1), signed_eea(3, 5));
    assert_eq!((-1, 2, 1), signed_eea(5, 3));
    assert_eq!((2, 1, -1), signed_eea(-3, 5));
    assert_eq!((-1, -2, 1), signed_eea(5, -3));
    assert_eq!((2, 1, 1), signed_eea(3, -5));
    assert_eq!((-1, -2, -1), signed_eea(-5, 3));
    assert_eq!((2, -1, -1), signed_eea(-3, -5));
    assert_eq!((-1, 2, -1), signed_eea(-5, -3));

    assert_eq!((0, 0, 0), signed_eea(0, 0));

    assert_eq!((1, 0, 4), signed_eea(4, 0));
    assert_eq!((0, 1, 4), signed_eea(0, 4));
    assert_eq!((1, 0, -4), signed_eea(-4, 0));
    assert_eq!((0, -1, 4), signed_eea(0, -4));
}

#[test]
fn test_signed_eea() {
    assert_eq!((-1, 1, 2), signed_eea(6, 8));
    assert_eq!((2, -1, 5), signed_eea(15, 25));
    assert_eq!((4, -7, 2), signed_eea(32, 18));
}

/// Finds the greatest common divisor of a and b
/// 
///  a < 0 => ggT(a, b) < 0
/// a > 0 => ggT(a, b) > 0
/// sign of b is irrelevant
/// ggT(0, 0) = 0
/// 
/// for details, see eea
pub fn gcd<Int: Integer + Copy>(a: Int, b: Int) -> Int {
    let (_, _, d) = signed_eea(a, b);
    return d;
}

#[test]
fn test_gcd() {
    assert_eq!(3, gcd(15, 6));
    assert_eq!(3, gcd(6, 15));

    assert_eq!(7, gcd(0, 7));
    assert_eq!(7, gcd(7, 0));
    assert_eq!(0, gcd(0, 0));

    assert_eq!(1, gcd(9, 1));
    assert_eq!(1, gcd(1, 9));

    assert_eq!(1, gcd(13, 300));
    assert_eq!(1, gcd(300, 13));

    assert_eq!(-3, gcd(-15, 6));
    assert_eq!(3, gcd(6, -15));
    assert_eq!(-3, gcd(-6, -15));
}