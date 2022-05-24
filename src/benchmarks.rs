use super::wrapper::*;
use super::fq::fq_small::*;
use super::poly::*;

use test::Bencher;

#[bench]
fn bench_poly_multiplication(b: &mut Bencher) {
    let poly_ring = PolyRing::adjoint(F1369, "x");
    let a = poly_ring.bind(poly_ring.from(poly_ring.base_ring().generator()));
    let x = poly_ring.bind(poly_ring.unknown());
    let f = (0..=100).map(|i| x.pow(i) * i as i64).sum::<RingElWrapper<&_>>();
    let g = (0..=100).map(|i| x.pow(i) * (100 - i) as i64 * &a).sum::<RingElWrapper<&_>>();
    let h = (0..=100).map(|n| x.pow(n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>()
         + (0..100).map(|n| x.pow(200 - n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>();
    b.iter(|| {
        assert_eq!(h, &f * &g);
    });
}

#[bench]
fn bench_poly_multiplication_huge(b: &mut Bencher) {
    let poly_ring = PolyRing::adjoint(F1369, "x");
    let a = poly_ring.bind(poly_ring.from(poly_ring.base_ring().generator()));
    let x = poly_ring.bind(poly_ring.unknown());
    let f = (0..=1000).map(|i| x.pow(i) * i as i64).sum::<RingElWrapper<&_>>();
    let g = (0..=1000).map(|i| x.pow(i) * (1000 - i) as i64 * &a).sum::<RingElWrapper<&_>>();
    let h = (0..=1000).map(|n| x.pow(n) * &a * ((1000 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>()
         + (0..1000).map(|n| x.pow(2000 - n) * &a * ((1000 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>();
    b.iter(|| {
        assert_eq!(h, &f * &g);
    });
}
