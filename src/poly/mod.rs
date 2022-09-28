use super::prelude::*;
use super::wrapper::*;
use monomial_order::*;

use std::borrow::Cow;

pub mod ops;
pub mod factoring;
pub mod uni_var;
pub mod multi_var;
pub mod sumation;
pub mod monomial_order;
pub mod buchberger;

use vector_map::VecMap;
pub use uni_var::*;
pub use multi_var::*;

pub fn normalized_monomial<I: DoubleEndedIterator<Item = usize>>(exponents: I) -> Vec<usize> {
    let mut result: Vec<usize> = exponents.rev().skip_while(|x| *x == 0).collect();
    result.reverse();
    return result;
}

pub trait MultiPolyRing: Ring + CanonicalIsomorphismInfo<MultivariatePolyRing<Self::BaseRing>> + RingExtension {

    fn unknown_count(&self) -> usize;
    fn get_name(&self, var: usize) -> &'static str;
    fn get_var(&self, name: &str) -> usize;
    fn as_poly(&self, var: usize) -> El<Self>;
    fn coefficient_at(&self, x: &El<Self>, monomial: &Vec<usize>) -> El<Self::BaseRing>;
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: El<Self>, values: Vec<&El<S>>, ring: &S) -> El<S>;
    fn deg(&self, var: usize, x: &El<Self>) -> Option<usize>;
    fn total_deg(&self, x: &El<Self>) -> Option<usize>;
    fn for_monomials<'a, F: FnMut(&Vec<usize>, &'a El<Self::BaseRing>)>(&'a self, x: &'a El<Self>, f: F);

    ///
    /// Computes an expression of the expression `lhs - rhs * m * a` where m is
    /// a monomial and a a scalar. We have a dedicated function for this to allow ring
    /// to provide very efficient implementations (in particular avoiding
    /// intermediate objects). This function is used very frequently during
    /// algorithms using elimination (Groebner basis), and so its performance
    /// is somewhat critical.
    /// 
    fn sub_scaled(&self, lhs: El<Self>, rhs: &El<Self>, monomial: &Vec<usize>, scalar: &El<Self::BaseRing>) -> El<Self> {
        self.sub(lhs, self.scale(rhs.clone(), monomial, scalar))
    }

    fn scale(&self, x: El<Self>, monomial: &Vec<usize>, scalar: &El<Self::BaseRing>) -> El<Self> {
        self.mul(
            self.mul(x, self.product(
                (0..monomial.len()).map(|i| self.pow(&self.as_poly(i), monomial[i] as u32)
            ))),
            self.embedding()(scalar.clone())
        )
    }

    fn lt<'a, M: MonomialOrder>(&'a self, x: &'a El<Self>, order: M) -> Option<(Cow<'a, Vec<usize>>, &'a El<Self::BaseRing>)> {
        let mut result: Option<(Cow<'a, Vec<usize>>, _)> = None;
        self.for_monomials(x, |m, c| if let Some((current_m, _)) = &result {
            if order.cmp(&m, current_m) == std::cmp::Ordering::Greater {
                result = Some((Cow::Owned(m.clone()), c));
            }
        } else {
            result = Some((Cow::Owned(m.clone()), c));
        });
        return result;
    }

    fn lm<'a, M: MonomialOrder>(&'a self, x: &'a El<Self>, order: M) -> Option<Cow<'a, Vec<usize>>> {
        self.lt(x, order).map(|x| x.0)
    }

    ///
    /// # Attention
    /// 
    /// The ring underlying the univariate result polynomial ring should be self. This
    /// is already strongly suggested by the lifetime contract, and allows cheaper elevating/
    /// de-elevating since no no ring objects have to be created. After elevating all variables,
    /// use then `as_constant` to access the value as element of the base ring.
    /// 
    fn elevate_var_ring<'a>(&'a self, var: usize) -> PolyRingImpl<&'a Self>;
    fn elevate_var<'a>(&'a self, var: usize, x: El<Self>) -> El<PolyRingImpl<&'a Self>>;
    fn de_elevate_var(&self, var: usize, x: El<PolyRingImpl<&Self>>) -> El<Self>;
}

impl<R, P> CanonicalEmbeddingInfo<MultivariatePolyRing<R>> for P
    where R: Ring, P: RingDecorator, P::DecoratedRing: CanonicalEmbeddingInfo<MultivariatePolyRing<R>>
{
    fn has_embedding(&self, from: &MultivariatePolyRing<R>) -> RingPropValue { self.decorated_ring().has_embedding(from) }
    fn embed(&self, from: &MultivariatePolyRing<R>, el: El<MultivariatePolyRing<R>>) -> Self::El  { self.decorated_ring().embed(from, el) }
}

impl<R, P> CanonicalIsomorphismInfo<MultivariatePolyRing<R>> for P
    where R: Ring, P: RingDecorator, P::DecoratedRing: CanonicalIsomorphismInfo<MultivariatePolyRing<R>>
{
    fn has_isomorphism(&self, from: &MultivariatePolyRing<R>) -> RingPropValue { self.decorated_ring().has_isomorphism(from) }
    fn preimage(&self, from: &MultivariatePolyRing<R>, el: El<Self>) -> El<MultivariatePolyRing<R>>  { self.decorated_ring().preimage(from, el) }
}

impl<R> MultiPolyRing for R
    where R: RingDecorator, R::DecoratedRing: MultiPolyRing
{
    fn unknown_count(&self) -> usize { self.decorated_ring().unknown_count() }
    fn get_name(&self, var: usize) -> &'static str { self.decorated_ring().get_name(var) }
    fn get_var(&self, name: &str) -> usize { self.decorated_ring().get_var(name) }
    fn as_poly(&self, var: usize) -> El<Self> { self.decorated_ring().as_poly(var) }
    fn coefficient_at(&self, x: &El<Self>, monomial: &Vec<usize>) -> El<Self::BaseRing> { self.decorated_ring().coefficient_at(x, monomial) }
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: El<Self>, values: Vec<&El<S>>, ring: &S) -> El<S> { self.decorated_ring().evaluate_at(f, values, ring) }
    fn deg(&self, var: usize, x: &El<Self>) -> Option<usize> { self.decorated_ring().deg(var, x) }
    fn total_deg(&self, x: &El<Self>) -> Option<usize> { self.decorated_ring().total_deg(x) }
    fn elevate_var_ring<'a>(&'a self, var: usize) -> PolyRingImpl<&'a Self> { PolyRingImpl::adjoint(self, self.get_name(var)) }
    fn for_monomials<'a, F: FnMut(&Vec<usize>, &'a El<Self::BaseRing>)>(&'a self, x: &'a El<Self>, f: F) { self.decorated_ring().for_monomials(x, f) }
    fn lt<'a, M: MonomialOrder>(&'a self, x: &'a El<Self>, order: M) -> Option<(Cow<'a, Vec<usize>>, &'a El<Self::BaseRing>)> { self.decorated_ring().lt(x, order) }
    fn sub_scaled(&self, lhs: El<Self>, rhs: &El<Self>, monomial: &Vec<usize>, scalar: &El<Self::BaseRing>) -> El<Self> { self.decorated_ring().sub_scaled(lhs, rhs, monomial, scalar) }
    fn scale(&self, x: El<Self>, monomial: &Vec<usize>, scalar: &El<Self::BaseRing>) -> El<Self> { self.decorated_ring().scale(x, monomial, scalar) }
    
    fn elevate_var<'a>(&'a self, var: usize, x: El<Self>) -> El<PolyRingImpl<&'a Self>> { 
        let result = self.decorated_ring().elevate_var(var, x);
        let map = self.elevate_var_ring(var).lift_hom::<_, &R::DecoratedRing>(|x| x);
        return map(result);
    }

    fn de_elevate_var(&self, var: usize, x: El<PolyRingImpl<&Self>>) -> El<Self> {
        let map = self.elevate_var_ring(var).lift_hom::<_, &R::DecoratedRing>(|x| x);
        self.decorated_ring().de_elevate_var(var, map(x))
    }
}

impl<R> WrappingRing<R>
    where R: MultiPolyRing
{
    pub fn unknown_count(&self) -> usize { self.wrapped_ring().unknown_count() }
    pub fn get_unknown(&self, i: usize) -> &'static str { self.wrapped_ring().get_name(i) }
    pub fn as_poly(&self, var: &str) -> El<Self> { self.from(self.wrapped_ring().as_poly(self.wrapped_ring().get_var(var))) }

    pub fn elevate_var_ring<'a>(&'a self, var: usize) -> WrappingRing<PolyRingImpl<&'a R>> {
        WrappingRing::new(self.wrapped_ring().elevate_var_ring(var))
    }
}

impl<R> RingElWrapper<R>
    where R: MultiPolyRing
{
    pub fn monomial_coefficient<'a>(&'a self, monomial: Vec<usize>) -> RingElWrapper<&'a R::BaseRing> {
        let coeff = self.parent_ring().coefficient_at(self.val(), &monomial);
        RingElWrapper::new(coeff, self.parent_ring().base_ring())
    }

    pub fn evaluate_at<'a, S>(self, values: &'a Vec<RingElWrapper<S>>) -> RingElWrapper<&'a S> 
        where S: Ring + CanonicalEmbeddingInfo<R::BaseRing>
    {
        let (el, ring) = self.destruct();
        assert_eq!(values.len(), ring.unknown_count());
        let result = ring.evaluate_at(el, values.iter().map(|x| x.val()).collect(), values[0].parent_ring());
        return RingElWrapper::new(result, values[0].parent_ring());
    }

    pub fn var_deg(&self, var: usize) -> Option<usize> { self.parent_ring().deg(var, self.val()) }
    pub fn total_deg(&self) -> Option<usize> { self.parent_ring().total_deg(self.val()) }

    pub fn elevate_var<'a>(&'a self, var: usize) -> RingElWrapper<PolyRingImpl<&'a R>>{
        RingElWrapper::new(self.parent_ring().elevate_var(var, self.val().clone()), self.parent_ring().elevate_var_ring(var))
    }
}

impl<'a, R> RingElWrapper<PolyRingImpl<&'a R>>
    where R: MultiPolyRing
{
    pub fn de_elevate_var<'b>(&'b self) -> RingElWrapper<&'b R> {
        let ring = self.parent_ring().base_ring().decorated_ring();
        RingElWrapper::new(ring.de_elevate_var(ring.get_var(self.parent_ring().unknwon_name()), self.val().clone()), ring)
    }
}

pub trait PolyRing: Ring + CanonicalIsomorphismInfo<PolyRingImpl<Self::BaseRing>> + RingExtension {

    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>>;
    fn deg(&self, x: &El<Self>) -> Option<usize>;
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>;
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>);
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: &El<Self>, x: El<S>, ring: &S) -> El<S>;
    fn derive(&self, el: El<Self>) -> El<Self>;
    fn unknown(&self) -> El<Self>;
    fn unknwon_name(&self) -> &str;

    fn normalize(&self, mut f: El<Self>) -> (El<Self>, El<Self::BaseRing>) {
        assert!(self.base_ring().is_field().can_use());
        let lc = self.lc(&f).unwrap().clone();
        let lc_inv = self.base_ring().div(self.base_ring().one(), &lc);
        self.scale(&mut f, &lc_inv);
        return (f, lc);
    }
}

impl<R, P> CanonicalEmbeddingInfo<PolyRingImpl<R>> for P
    where R: Ring, P: RingDecorator, P::DecoratedRing: CanonicalEmbeddingInfo<PolyRingImpl<R>>
{
    fn has_embedding(&self, from: &PolyRingImpl<R>) -> RingPropValue { self.decorated_ring().has_embedding(from) }
    fn embed(&self, from: &PolyRingImpl<R>, el: El<PolyRingImpl<R>>) -> Self::El  { self.decorated_ring().embed(from, el) }
}

impl<R, P> CanonicalIsomorphismInfo<PolyRingImpl<R>> for P
    where R: Ring, P: RingDecorator, P::DecoratedRing: CanonicalIsomorphismInfo<PolyRingImpl<R>>
{
    fn has_isomorphism(&self, from: &PolyRingImpl<R>) -> RingPropValue { self.decorated_ring().has_isomorphism(from) }
    fn preimage(&self, from: &PolyRingImpl<R>, el: El<Self>) -> El<PolyRingImpl<R>>  { self.decorated_ring().preimage(from, el) }
}

impl<P> PolyRing for P
    where P: RingDecorator, P::DecoratedRing: PolyRing
{
    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>> { self.decorated_ring().lc(x) }
    fn unknwon_name(&self) -> &str { self.decorated_ring().unknwon_name() }
    fn deg(&self, x: &El<Self>) -> Option<usize>  { self.decorated_ring().deg(x) }
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>  { self.decorated_ring().coefficient_at(x, i) }
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>) { self.decorated_ring().scale(x, coeff) }
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: &El<Self>, x: El<S>, ring: &S) -> El<S> { self.decorated_ring().evaluate_at(f, x, ring) }
    fn derive(&self, el: El<Self>) -> El<Self> { self.decorated_ring().derive(el) }
    fn normalize(&self, f: El<Self>) -> (El<Self>, El<Self::BaseRing>) { self.decorated_ring().normalize(f) }
    fn unknown(&self) -> El<Self> { self.decorated_ring().unknown() }
}

impl<S, P> FnOnce<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    type Output = RingElWrapper<S>;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<S, P> FnMut<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<S, P> Fn<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        let (x, ring) = x.destruct();
        RingElWrapper::new(self.parent_ring().evaluate_at(self.val(), x, &ring), ring)
    }
}

impl<P> RingElWrapper<P>
    where P: PolyRing
{
    pub fn lc<'a>(&'a self) -> Option<RingElWrapper<&'a P::BaseRing>> {
        self.parent_ring().lc(self.val()).map(|x| RingElWrapper::new(x, self.parent_ring().base_ring()))
    }

    pub fn deg(&self) -> Option<usize> {
        self.parent_ring().deg(self.val())
    }

    pub fn coefficient_at<'a>(&'a self, i: usize) -> RingElWrapper<&'a P::BaseRing> {
        RingElWrapper::new(self.parent_ring().coefficient_at(self.val(), i), self.parent_ring().base_ring())
    }

    pub fn roots(self) -> VecMap<RingElWrapper<P::BaseRing>, usize>
        where P: UfdInfoRing
    {
        self.factor().into_iter().filter_map(|(f, e)| {
            if f.deg() == Some(1) {
                Some(((-f.coefficient_at(0) / f.coefficient_at(1)).clone_ring(), e))
            } else {
                None
            }
        }).collect()
    }

    pub fn normalize(&mut self) {
        let coeff = self.lc().unwrap().inv().clone_ring();
        self.parent_ring().clone().scale(&mut self.val_mut(), coeff.val());
    }

    pub fn scale(&mut self, coeff: &RingElWrapper<P::BaseRing>) {
        self.parent_ring().clone().scale(&mut self.val_mut(), coeff.val());
    }

    pub fn scaled(mut self, coeff: &RingElWrapper<P::BaseRing>) -> RingElWrapper<P> {
        self.scale(coeff);
        return self;
    }
}

impl<P> WrappingRing<P>
    where P: PolyRing
{
    pub fn unknown(&self) -> El<Self> {
        self.from(self.wrapped_ring().unknown())
    }
}