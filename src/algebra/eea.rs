use super::super::alg::*;

use std::mem::swap;

///
/// For a, b computes s, t, d such that `s*a + t*b == d` is a greatest 
/// common divisor of a and b. d is only unique up to units, and s, t 
/// are not unique at all. No guarantees are given on which
/// of these solutions is returned. For integers, see signed_eea 
/// which gives more guarantees.
/// 
/// The given ring must be euclidean
/// 
pub fn eea<R: Ring>(ring: &R, fst: R::El, snd: R::El) -> (R::El, R::El, R::El) 
    where R::El: Clone
{
    assert!(ring.is_euclidean());

    let (mut a, mut b) = (fst, snd);
    let (mut sa, mut ta) = (ring.one(), ring.zero());
    let (mut sb, mut tb) = (ring.zero(), ring.one());

    // invariant: sa * a + ta * b = fst, sb * a + tb * b = snd
    while !ring.eq(&b, &ring.zero()) {
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
/// In particular, have `signed_gcd(6, 8) == 2`, `signed_gcd(0, 0) == 0`, 
/// `signed_gcd(0, x) == |x|`, `signed_gcd(x, 0) == x`, 
/// `signed_gcd(-1, 1) == -1`, `signed_gcd(1, -1) == 1`
/// and therefore `signed_eea(6, 8) == (-1, 1, 2)`, 
/// `signed_eea(-6, 8) == (-1, -1, -2)`, 
/// `signed_eea(8, -6) == (1, 1, 2)`, 
/// `signed_eea(0, 0) == (0, 0, 0)`
/// 
pub fn signed_eea<Int: Integer>(fst: Int, snd: Int) -> (Int, Int, Int) {
    if fst == Int::zero() {
        if snd == Int::zero() {
            return (Int::zero(), Int::zero(), Int::zero());
        } else if snd < Int::zero() {
            return (Int::zero(), -Int::one(), -snd);
        } else {
            return (Int::zero(), Int::one(), snd);
        }
    }
    let fst_negative = fst < Int::zero();

    let (s, t, d) = eea(&Int::RING, fst, snd);
    
    // the sign is not consistent (potentially toggled each iteration), 
    // so normalize here
    if (d < Int::zero()) == fst_negative {
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
pub fn gcd<R: Ring>(ring: &R, a: R::El, b: R::El) -> R::El 
    where R::El: Clone
{
    assert!(ring.is_euclidean());
    
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
pub fn signed_gcd<Int: Integer>(a: Int, b: Int) -> Int {
    let (_, _, d) = signed_eea(a, b);
    return d;
}

#[test]
fn test_gcd() {
    assert_eq!(3, signed_gcd(15, 6));
    assert_eq!(3, signed_gcd(6, 15));

    assert_eq!(7, signed_gcd(0, 7));
    assert_eq!(7, signed_gcd(7, 0));
    assert_eq!(0, signed_gcd(0, 0));

    assert_eq!(1, signed_gcd(9, 1));
    assert_eq!(1, signed_gcd(1, 9));

    assert_eq!(1, signed_gcd(13, 300));
    assert_eq!(1, signed_gcd(300, 13));

    assert_eq!(-3, signed_gcd(-15, 6));
    assert_eq!(3, signed_gcd(6, -15));
    assert_eq!(-3, signed_gcd(-6, -15));
}

#[test]
fn exp() {
    use super::fractions::*;
    use super::poly::*;
    use super::rat::*;
    use super::super::alg_env::*;

    let mut var_ring = MultivariatePolyRing::new(r64::RING);
    let A = var_ring.adjoint("A");
    let B = var_ring.adjoint("B");
    let base_ring = FieldOfFractions::new(var_ring);
    let A = base_ring.from(A);
    let B = base_ring.from(B);
    let coord_ring_x = PolyRing::adjoint(base_ring, "x");
    let A = coord_ring_x.from(A);
    let B = coord_ring_x.from(B);
    let x = coord_ring_x.unknown();
    let A = coord_ring_x.bind::<RingAxiomsEuclideanRing>(A);
    let B = coord_ring_x.bind(B);
    let x = coord_ring_x.bind(x);
    
    let f = x.clone().pow(3) + A.clone() * x.clone() + B.clone();
    let g = x.clone().pow(4) - A.clone() * x.clone().pow(2) * 2 - B * x.clone() * 8 + A.pow(2);
    println!("{}", g);
    println!("{}", f);
    println!("{}", g.clone() % f.clone());
    println!("{}", f.clone() % (g.clone() % f.clone()));
    println!("{}", (g.clone() % f.clone()) % (f.clone() % (g.clone() % f.clone())));
    println!("{}", f.ring().display(&gcd(&f.ring(), f, g)));
    assert!(false);
}