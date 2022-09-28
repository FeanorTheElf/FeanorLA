use super::super::prelude::*;
use super::*;

pub struct IntegralCubic<'a, R: IntegerRing + OrderedRing> {
    p: &'a R::El,
    q: &'a R::El,
    ring: R
}

impl<'a, R: IntegerRing + OrderedRing + EuclideanInfoRing> IntegralCubic<'a, R> {

    pub fn new(p: &'a R::El, q: &'a R::El, ring: R) -> Self {
        IntegralCubic {
            p, q, ring
        }
    }

    pub fn eval(&self, x: &R::El) -> R::El {
        self.ring.add_ref(self.ring.add(self.ring.pow(x, 3), self.ring.mul_ref(x, &self.p)), &self.q)
    }

    ///
    /// Calls the given closure on all integral roots of the cubic.
    /// The closure is only called once, even if the root is a multiple roots.
    /// 
    /// The algorithm relies on bisection, and works reliably with very large numbers.
    /// 
    fn calc_integral_roots<F>(&self, mut f: F)
        where F: FnMut(R::El)
    {
        let Z = &self.ring;
        let mut process_potential_root = |x: R::El| if Z.is_zero(&self.eval(&x)) { f(x) };
        let diff_disc = Z.mul(Z.from_z(-12), self.p.clone());
        if Z.cmp(&diff_disc, &Z.zero()) != Ordering::Greater {
            // zero or one maximum/minimum, so there can be at most one root
            process_potential_root(find_zero_floor(&Z, |x| self.eval(x), Z.zero()));
        } else {
            let root_size_bound = if Z.abs_cmp(&self.p, &self.q) == std::cmp::Ordering::Less {
                Z.abs(self.q.clone())
            } else {
                Z.abs(self.p.clone())
            };
            let extremum_floor = Z.floor_div(Z.root_floor(&diff_disc, 2), &Z.from_z(6));
            // on the intervals [a0, a1], [b0, b1], [c0, c1], the function is monotonous, 
            // hence has at most one root
            let a0 = Z.neg(root_size_bound.clone());
            let a1 = Z.neg(extremum_floor.clone());
            let b0 = Z.sub_ref_fst(&a1, Z.one());
            let b1 = extremum_floor;
            let c0 = Z.add_ref(Z.one(), &b1);
            let c1 = root_size_bound;
            let a0_val = self.eval(&a0);
            let a1_val = self.eval(&a1);
            let b0_val = self.eval(&b0);
            let b1_val = self.eval(&b1);
            let c0_val = self.eval(&c0);
            let c1_val = self.eval(&c1);
            let zero = Z.zero();
            if Z.cmp(&a0_val, &zero) != Ordering::Greater && Z.cmp(&a1_val, &zero) != Ordering::Less {
                process_potential_root(bisect(&Z, |x| self.eval(&x), a0, a1));
            }
            if Z.cmp(&b0_val, &zero) != Ordering::Less && Z.cmp(&b1_val, &zero) != Ordering::Greater {
                process_potential_root(bisect(&Z, |x| Z.neg(self.eval(&x)), b0, b1));
            }
            if Z.cmp(&c0_val, &zero) != Ordering::Greater && Z.cmp(&c1_val, &zero) != Ordering::Less {
                process_potential_root(bisect(&Z, |x| self.eval(&x), c0, c1));
            }
        }
    }

    pub fn find_integral_roots(&self) -> impl Iterator<Item = R::El> {
        let mut roots = [None, None, None];
        let mut i = 0;
        self.calc_integral_roots(|root| {
            roots[i] = Some(root);
            i += 1;
        });
        return <[Option<_>; 3] as std::iter::IntoIterator>::into_iter(roots).filter_map(|x| x);
    }
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
pub fn find_zero_floor<R, F>(ring: &R, mut f: F, approx: El<R>) ->  El<R>
    where R: IntegerRing, F: FnMut(& El<R>) ->  El<R>
{
    let mut begin = approx.clone();
    let mut step = ring.one();
    let two = ring.from_z(2);
    while ring.cmp(&f(&begin), &ring.zero()) == Ordering::Greater {
        begin = ring.sub_ref_snd(begin, &step);
        ring.mul_assign(&mut step, &two.clone());
    }
    let mut end = approx;
    step = ring.one();
    while ring.cmp(&f(&end), &ring.zero()) == Ordering::Less {
        end = ring.add_ref(end, &step);
        ring.mul_assign(&mut step, &two.clone());
    }
    return bisect(ring, f, begin, end);
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
pub fn bisect<R, F>(ring: &R, mut f: F, mut start: El<R>, mut end: El<R>) -> El<R>
    where R: IntegerRing, F: FnMut(&El<R>) -> El<R>
{
    assert!(!ring.is_pos(&f(&start)));
    assert!(!ring.is_neg(&f(&end)));
    if ring.is_zero(&f(&end)) {
        return end;
    }
    let two = ring.from_z(2);
    loop {
        let mid = ring.floor_div(ring.add_ref(start.clone(), &end), &two);
        if ring.is_eq(&mid, &start) {
            return start;
        }
        match ring.cmp(&f(&mid), &ring.zero()) {
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

#[test]
fn test_find_zero_floor() {
    let f = |x: &BigInt| BigInt::RING.sub(BigInt::RING.mul_ref(x, x), BigInt::from(234867));
    assert_eq!(BigInt::from(484), find_zero_floor(&BigInt::RING, f, BigInt::ZERO));

    let f = |x: &BigInt| x.clone();
    assert_eq!(BigInt::ZERO, find_zero_floor(&BigInt::RING, f, BigInt::ZERO));
}

#[test]
fn test_find_zero_floor_i64() {
    let f = |x: &i64| x.pow(3) - *x + 478;
    assert_eq!(-8, find_zero_floor(&i64::RING, f, 0));
}

#[test]
fn test_integral_cubic_eval() {
    let i = z_hom(&BigInt::RING);
    assert_eq!(
        i(1),
        IntegralCubic { p: &i(1), q: &i(1), ring: &BigInt::RING }.eval(&BigInt::ZERO)
    );
    assert_eq!(
        i(-3),
        IntegralCubic { p: &i(-2), q: &i(1), ring: &BigInt::RING }.eval(&i(-2))
    )
}

#[test]
fn test_find_integral_roots() {
    let i = z_hom(&BigInt::RING);
    assert_eq!(
        Vec::<BigInt>::new(),
        IntegralCubic { p: &i(1), q: &i(1), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
    assert_eq!(
        vec![i(1)],
        IntegralCubic { p: &i(-2), q: &i(1), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
    assert_eq!(
        vec![i(-3), i(1), i(2)],
        IntegralCubic { p: &i(-7), q: &i(6), ring: &BigInt::RING }.find_integral_roots().collect::<Vec<_>>()
    );
}
