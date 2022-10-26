use crate::la::vec::*;

use super::prelude::*;

pub mod karatsuba;

///
/// Trait for rings to provide information on which convoluted multiplication algorithms
/// will be most efficient.
/// 
pub trait MulAlgorithmLimit : Ring {

    ///
    /// Gives an approximate threshold at which naive convoluted nxn multiplication becomes
    /// slower than karatsuba nxn multiplication.
    /// 
    fn karatsuba_threshold(&self) -> usize;
}

fn naive_convoluted_mul<R, U, V, W>(
    lhs: Vector<U, R::El>,
    rhs: Vector<V, R::El>,
    mut out: Vector<W, R::El>,
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
    where F: FnMut(Vector<&Vec<i128>, i128>, Vector<&Vec<i128>, i128>, Vector<&mut Vec<i128>, i128>)
{
    let a = Vector::from_fn(n, |i| i as i128);
    let b = Vector::from_fn(n, |i| 2 * i as i128);
    let mut c = Vector::zero(2 * n).into_owned();
    let expected = (n * n * (n + 1) - n * (n + 1) * (2 * n + 1) / 3) as i128;
    bencher.iter(|| {
        c.assign(Vector::zero(2 * n));
        f(a.as_ref(), b.as_ref(), c.as_mut());
        assert_eq!(expected, *c.at(n));
    });
}

#[bench]
fn bench_karatsuba_i128_100(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| karatsuba::karatsuba_mul(l, r, out, &i128::RING), 
        100
    );
}

#[bench]
fn bench_naive_i128_100(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| naive_convoluted_mul(l, r, out, &i128::RING), 
        100
    );
}

#[bench]
fn bench_karatsuba_i128_1000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| karatsuba::karatsuba_mul(l, r, out, &i128::RING), 
        1000
    );
}

#[bench]
fn bench_naive_i128_1000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| naive_convoluted_mul(l, r, out, &i128::RING), 
        1000
    
    );
}
#[bench]
fn bench_karatsuba_i128_10000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| karatsuba::karatsuba_mul(l, r, out, &i128::RING), 
        10000
    );
}

#[bench]
fn bench_naive_i128_10000(bencher: &mut test::Bencher) {
    bench_convoluted_mul_impl_i128(
        bencher, 
        |l, r, out| naive_convoluted_mul(l, r, out, &i128::RING), 
        10000
    );
}