use super::alg::Integer;

use std::mem::swap;

pub fn eea<Int: Integer>(fst: Int, snd: Int) -> (Int, Int) {
    let (mut a, mut b) = if fst > snd { (fst, snd) } else { (snd, fst) };
    let (mut sa, mut ta) = (1, 0);
    let (mut sb, mut tb) = (0, 1);
    while b != 0 {
        ta -= a / b * tb;
        sa -= a / b * sb;
        a = a % b;
        swap(&mut a, &mut b);
        swap(&mut sa, &mut sb);
        swap(&mut ta, &mut tb);
    }
    return if fst > snd { (sa, ta) } else { (ta, sa) };
}

pub fn gcd<Int: Integer>(a: Int, b: Int) -> Int {
    let (s, t) = eea(a, b);
    return s * a + t * b;
}