use crate::la::vec::*;

use super::prelude::*;

pub mod karatsuba;

fn naive_convoluted_mul<R, U, V, W>(
    mut out: Vector<W, R::El>,
    lhs: Vector<U, R::El>,
    rhs: Vector<V, R::El>,
    ring: &R
)
    where R: Ring, U: VectorView<R::El>, V: VectorView<R::El>, W: VectorViewMut<R::El>
{
    assert!(out.len() >= lhs.len() + rhs.len() - 1);
    for i in 0..lhs.len() {
        for j in 0..rhs.len() {
            ring.add_assign(out.at_mut(i + j), ring.mul_ref(lhs.at(i), rhs.at(j)));
        }
    }
}

#[cfg(test)]
pub fn bench_convoluted_mul_impl_i128<F>(bencher: &mut test::Bencher, mut f: F, n: usize) 
    where F: FnMut(&Vec<i128>, &Vec<i128>, &mut Vec<i128>)
{
    let a = (0..n).map(|i| i as i128).collect::<Vec<_>>();
    let b = (0..n).map(|i| 2 * i as i128).collect::<Vec<_>>();
    let mut c = Vec::new();
    c.resize(2 * n, 0);
    let expected = (n * n * (n + 1) - n * (n + 1) * (2 * n + 1) / 3) as i128;
    bencher.iter(|| {
        c.clear();
        c.resize(2 * n, 0);
        f(&a, &b, &mut c);
        assert_eq!(expected, *c.at(n));
    });
}

#[bench]
fn bench_karatsuba_i128_100(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| karatsuba::karatsuba::<_, 4>(&mut out[..], &l[..], &r[..], &i128::RING), 
        100
    );
}

#[bench]
fn bench_naive_i128_100(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| naive_convoluted_mul(Vector::new(out), Vector::new(l), Vector::new(r), &i128::RING), 
        100
    );
}

#[bench]
fn bench_karatsuba_i128_1000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| karatsuba::karatsuba::<_, 4>(&mut out[..], &l[..], &r[..], &i128::RING), 
        1000
    );
}

#[bench]
fn bench_naive_i128_1000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| naive_convoluted_mul(Vector::new(out), Vector::new(l), Vector::new(r), &i128::RING), 
        1000
    
    );
}