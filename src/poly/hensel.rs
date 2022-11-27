use super::super::prelude::*;
use super::super::fq::*;
use super::super::eea::*;
use super::*;

pub fn map_poly<P, S, F>(base_ring: &P, f: &El<P>, target_ring: &S, mut map_fn: F) -> El<S>
    where P: PolyRing, S: PolyRing, F: FnMut(El<P::BaseRing>) -> El<S::BaseRing> 
{
    if base_ring.is_zero(f) {
        return target_ring.zero();
    }
    return target_ring.from_coefficients(
        (0..=base_ring.deg(f).unwrap()).map(|i| map_fn(base_ring.coefficient_at(f, i)))
    );
}

fn reduce_bezout<P>(ring: &P, f: &El<P>, g: &El<P>, mut s: El<P>, mut t: El<P>) -> (El<P>, El<P>) 
    where P: PolyRing + EuclideanInfoRing
{
    let (quo, rem) = ring.euclidean_div_rem(s, g);
    s = rem;
    let mut gf_coeff = quo;
    let (quo, rem) = ring.euclidean_div_rem(t, f);
    t = rem;
    gf_coeff = ring.add(gf_coeff, quo);
    s = ring.add(s, ring.mul_ref(&gf_coeff, &g));
    return (s, t);
}

///
/// Performs hensel lifting, i.e. for a factorization `h = fg mod p^e` computes a factorization
/// `h = f' g' mod p^f` for `f > e` such that `f' = f, g' = g mod p^e`. Note that this can fail
/// if `f` is not coprime to `g` in `Fp[X]`, in which case `Err(())` is returned.
/// 
/// This function performs operations also in the polynomial ring over the prime field, `Fp[X]`.
/// To allow controlling this behavior, it takes this ring as a final parameter. Note that in the
/// common use case that we want to lift a factorization of `h` in `Fp[X]` to some ring `Z/p^eZ[X]`,
/// this can just be the same as the base ring.
/// 
/// # Uniqueness
/// 
/// If `h`, `f` and `g` are monic, and `f'`, `g'` are required to be monic, this factorization is
/// unique. In this case, the function will also produce monic `f'`, `g'`. However, if this is
/// not the case, there might be multiple lifts with above properties, and `f'` resp. `g'` can
/// even have larger degree than `f` and `g`. In such a situation, the output is arbitrary, but
/// at least we ensure that `deg(f') + deg(g') == deg(h)`.
/// 
/// # Examples
/// ```
/// # use feanor_la::prelude::*;
/// # use feanor_la::poly::*;
/// # use feanor_la::poly::hensel::*;
/// # use feanor_la::fq::zn_small::*;
/// # use feanor_la::wrapper::*;
/// // we define the ring hierarchy
/// let Fp = WrappingRing::new(ZnEl::<7>::RING);
/// let S = WrappingRing::new(ZnEl::<{7 * 7 * 7}>::RING);
/// let FpX = Fp.std_poly_ring("X");
/// let SX = S.std_poly_ring("X");
/// 
/// // now we define the polynomials
/// let x = FpX.unknown();
/// let f = x.pow(2) + 1;
/// let g = x.pow(2) + 3;
/// let x = SX.unknown();
/// let h = x.pow(4) + x.pow(3) * 49 + x.pow(2) * 11 + &x * 147 + 10;
/// 
/// // we have h = fg mod 7, so can do hensel lifting
/// let (lifted_f, lifted_g) = hensel_lift(FpX.wrapped_ring(), SX.wrapped_ring(), f.into_val(), g.into_val(), h.val(), FpX.wrapped_ring()).unwrap();
/// 
/// // and the result is correct
/// let lifted_f = SX.wrap(lifted_f);
/// let lifted_g = SX.wrap(lifted_g);
/// assert_eq!(h, lifted_f * lifted_g);
/// ```
/// 
#[allow(non_snake_case)]
pub fn hensel_lift<P, S, T>(RX: &P, SX: &S, f: El<P>, g: El<P>, h: &El<S>, FpX: &T) -> Result<(El<S>, El<S>), ()>
    where P: PolyRing, P::BaseRing: IntegerQuotientRing, S: PolyRing + DivisibilityInfoRing, S::BaseRing: IntegerQuotientRing, T: PolyRing + EuclideanInfoRing, T::BaseRing: IntegerQuotientRing
{
    let Fp = FpX.base_ring();
    let R = RX.base_ring();
    let S = SX.base_ring();
    let (p_base, e_base) = R.characteristic().is_prime_power().unwrap();
    let (p_ext, e_ext) = S.characteristic().is_prime_power().unwrap();
    let (p_fp, e_fp) = FpX.characteristic().is_prime_power().unwrap();
    assert_eq!(p_base, p_ext);
    assert_eq!(p_base, p_fp);
    assert_eq!(1, e_fp);
    assert!(e_base < e_ext);
    assert!(e_base >= 1);
    let p = p_base.into_val().to_i64().unwrap();

    // define lots of lifting functions
    let fp_lifting_ring = Fp.lifting_ring();
    let rbase_lifting_ring = R.lifting_ring();
    let rext_lifting_ring = S.lifting_ring();
    let R_Fp = |x: El<P::BaseRing>| Fp.from_z_gen(R.lift(&x, &rbase_lifting_ring), &rbase_lifting_ring);
    let S_R = |x: El<S::BaseRing>| R.from_z_gen(S.lift(&x, &rext_lifting_ring), &rext_lifting_ring);
    let S_Fp = |x: El<S::BaseRing>| Fp.from_z_gen(S.lift(&x, &rext_lifting_ring), &rext_lifting_ring);
    let R_S = |x: El<P::BaseRing>| S.from_z_gen(R.lift(&x, &rbase_lifting_ring), &rbase_lifting_ring);
    let Fp_S = |x: El<T::BaseRing>| S.from_z_gen(Fp.lift(&x, &fp_lifting_ring), &fp_lifting_ring);
    let RX_FpX = |x: &El<P>| map_poly(RX, x, FpX, R_Fp);
    let SX_RX = |x: &El<S>| map_poly(SX, x, RX, S_R);
    let SX_FpX = |x: &El<S>| map_poly(SX, x, FpX, S_Fp);
    let RX_SX = |x: &El<P>| map_poly(RX, x, SX, R_S);
    let FpX_SX = |x: &El<T>| map_poly(FpX, x, SX, Fp_S);

    assert!(RX.is_eq(&RX.mul_ref(&f, &g), &SX_RX(&h)));

    let reduced_f = RX_FpX(&f);
    let reduced_g = RX_FpX(&g);
    let (mut s, mut t, d) = eea(&FpX, reduced_f.clone(), reduced_g.clone());
    if !FpX.is_unit(&d) {
        return Err(());
    }
    FpX.scale(&mut s, &Fp.div(Fp.one(), &FpX.coefficient_at(&d, 0)));
    FpX.scale(&mut t, &Fp.div(Fp.one(), &FpX.coefficient_at(&d, 0)));
    
    let mut f = RX_SX(&f);
    let mut g = RX_SX(&g);
    for e in e_base..e_ext {
        let pe = S.pow(&S.from_z(p), e as u32);
        let delta_h = SX_FpX(&SX.quotient(&SX.sub_ref_fst(&h, SX.mul_ref(&f, &g)), &SX.from(pe.clone())).unwrap());
        let mut f1 = FpX.mul_ref(&t, &delta_h);
        let mut g1 = FpX.mul_ref(&s, &delta_h);
        (f1, g1) = reduce_bezout(FpX, &reduced_g, &reduced_f, f1, g1);
        f = SX.add(f, SX.scaled(FpX_SX(&f1), &pe));
        g = SX.add(g, SX.scaled(FpX_SX(&g1), &pe));
    }
    return Ok((f, g));
}

#[cfg(test)]
use super::super::fq::zn_small::*;

#[test]
fn test_hensel_lift() {
    let base_ring = WrappingRing::new(PolyRingImpl::adjoint(ZnEl::<5>::RING, "X"));
    let target_ring = WrappingRing::new(PolyRingImpl::adjoint(ZnEl::<25>::RING, "X"));
    let reduce = |x: &El<WrappingRing<PolyRingImpl<StaticRing<ZnEl<25>>>>>| map_poly(target_ring.wrapped_ring(), x.val(), base_ring.wrapped_ring(), 
        |y| base_ring.wrapped_ring().base_ring().from_z(target_ring.wrapped_ring().base_ring().lift(&y, &i64::RING))
    );
    let x = target_ring.unknown();
    let f = x.pow(3) + x.pow(2) * 100 + &x * 7 + 5;
    let g = x.pow(2) + &x * 5 + 6;
    let h = &f * &g;
    let (f_actual, g_actual) = hensel_lift(base_ring.wrapped_ring(), target_ring.wrapped_ring(), reduce(&f), reduce(&g), h.val(), base_ring.wrapped_ring()).unwrap();
    let f_actual = target_ring.wrap(f_actual);
    let g_actual = target_ring.wrap(g_actual);
    assert_eq!(f, f_actual);
    assert_eq!(g, g_actual);
}