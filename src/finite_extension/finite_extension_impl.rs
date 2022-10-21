use super::super::prelude::*;
use super::super::poly::*;
use super::super::la::inversion::*;
use super::super::poly::ops::poly_format;
use super::super::fq::*;
use super::super::combinatorics::iters::*;
use super::*;
use super::super::integer::*;
use super::super::la::vec::*;
use super::super::la::mat::*;

use std::marker::PhantomData;
use std::iter::FromIterator;

#[derive(Clone)]
pub struct FiniteExtensionImpl<R, V = VectorOwned<El<R>>, W = VectorOwned<El<R>>>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El>
{
    base_ring: R,
    mipo_values: Vector<V, R::El>,
    element_vector: PhantomData<W>,
    is_integral_cache: RingPropValueCache,
    gen_name: &'static str
}

impl<R> FiniteExtensionImpl<R, VectorOwned<R::El>, VectorOwned<R::El>>
    where R: Ring
{
    pub fn adjoin_element<P>(base_ring: R, mipo: El<P>, poly_ring: &P, gen_name: &'static str) -> Self
        where P: PolyRing<BaseRing = R>
    {
        let ring = PolyRingImpl::adjoint(base_ring.clone(), "X");
        let mipo = poly_ring.preimage(&ring, mipo);
        assert!(base_ring.is_one(&ring.lc(&mipo).unwrap()), "Only adjoining monic polynomials is supported");
        let degree = mipo.iter().enumerate().filter(|(_, x)| !base_ring.is_zero(x)).map(|(i, _)| i).max().unwrap();
        let mipo_values = Vector::new(
            mipo.raw_data()
                .into_iter()
                .map(|x| base_ring.neg(x))
                .take(degree)
                .collect::<Vec<_>>()
        );
        return Self::new(base_ring, mipo_values, gen_name);
    }
}

impl<R, V, W> FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    ///
    /// Creates the ring R[X]/(f) where 
    /// ```text
    /// f = mipo_values[0] + mipo_values[1] * X + ... + mipo_values[n - 1] * X^(n - 1) - X^n
    /// ```
    /// 
    pub const fn new(base_ring: R, mipo_values: Vector<V, R::El>, gen_name: &'static str) -> Self {
        FiniteExtensionImpl {
            base_ring: base_ring, 
            mipo_values: mipo_values,
            element_vector: PhantomData,
            is_integral_cache: RingPropValueCache::new(),
            gen_name: gen_name
        }
    }

    fn assert_valid_element(&self, el: &El<Self>) {
        assert_eq!(el.len(), self.rank_as_module());
    }

    pub fn polynomial_repr(&self, poly_ring: &PolyRingImpl<&R>, el: El<Self>) -> El<PolyRingImpl<R>> {
        el.into_owned().raw_data().into_iter().enumerate()
            .map(|(i, x)| poly_ring.mul(poly_ring.pow(&poly_ring.unknown(), i as u32), poly_ring.from(x)))
            .fold(poly_ring.zero(), |a, b| poly_ring.add(a, b))
    }

    fn create_multiplication_matrix(&self, el: El<Self>) -> Matrix<MatrixOwned<R::El>, R::El> {
        let d = self.rank_as_module();
        let mut matrix = Matrix::zero_ring(d, d, &self.base_ring).into_owned();
        for (j, x) in el.into_owned().raw_data().to_vec().into_iter().enumerate() {
            *matrix.at_mut(j, 0) = x;
        }
        for i in 1..d {
            for j in 1..d {
                *matrix.at_mut(j, i) = matrix.at(j - 1, i - 1).clone();
            }
            let last_el = matrix.at(d - 1, i - 1).clone();
            for j in 0..d {
                self.base_ring().add_assign(matrix.at_mut(j, i), self.base_ring().mul_ref(&self.mipo_values[j], &last_el));
            }
        }
        return matrix;
    }

    ///
    /// Returns the (monic) polynomial f such that this ring is isomorphic to
    /// `base_ring[X] / (f)`.
    /// 
    pub fn generating_polynomial(&self, poly_ring: &PolyRingImpl<&R>) -> El<PolyRingImpl<&R>> {
        self.mipo_values.iter()
            .cloned()
            .chain(std::iter::once(self.base_ring.neg(self.base_ring.one())))
            .scan(poly_ring.one(), |state, coeff| {
                let result = poly_ring.mul_ref(state, &poly_ring.from(coeff));
                poly_ring.mul_assign(state, &poly_ring.unknown());
                return Some(result);
            }).fold(poly_ring.zero(), |a, b| poly_ring.add(a, b))
    }

    pub fn in_base_ring(&self, el: &El<Self>) -> Option<El<R>> {
        if el.iter().skip(1).all(|x| self.base_ring().is_zero(x)) {
            Some(el.at(0).clone())
        } else {
            None
        }
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<R> for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring.has_embedding(from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(from, el))
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<&R> for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, from: &&R) -> RingPropValue {
        self.base_ring().has_embedding(*from)
    }

    fn embed(&self, from: &&R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(*from, el))
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<R> for FiniteExtensionImpl<&R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring().has_embedding(&from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(&from, el))
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<FiniteExtensionImpl<R, V, W>> for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, from: &Self) -> RingPropValue {
        if !self.base_ring().has_embedding(from.base_ring()).can_use() {
            return self.base_ring().has_embedding(from.base_ring());
        } else if self.rank_as_module() == from.rank_as_module() && self.mipo_values.as_ref().eq(from.mipo_values.as_ref(), self.base_ring()) {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    fn embed(&self, from: &Self, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        Vector::new(el.iter().cloned().map(|x| self.base_ring().embed(from.base_ring(), x)).collect())
    }
}

impl<R, V, W> CanonicalIsomorphismInfo<FiniteExtensionImpl<R, V, W>> for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_isomorphism(&self, from: &Self) -> RingPropValue {
        if !self.base_ring().has_isomorphism(from.base_ring()).can_use() {
            return self.base_ring().has_isomorphism(from.base_ring());
        } else if self.rank_as_module() == from.rank_as_module() && self.mipo_values.as_ref().eq(from.mipo_values.as_ref(), self.base_ring()) {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    fn preimage(&self, from: &Self, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        Vector::new(el.iter().cloned().map(|x| self.base_ring().preimage(from.base_ring(), x)).collect())
    }
}

impl<R, V, W> RingBase for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    type El = Vector<W, R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(&lhs);
        self.assert_valid_element(rhs);

        lhs.add_assign(rhs.as_ref(), &self.base_ring);
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(lhs);
        self.assert_valid_element(rhs);
        // I have somewhat optimized the following code: In the "simple" case that the
        // degree is a compile-time constant, W is array-base and the base ring multiplication
        // is simple as well, the generated assembly vnow looks pretty close to optimal 
        // - no memory allocations, no loops, just some register arithmetic

        // when we start, result contains the upper half of the product lhs * rhs as
        // polynomials; then we reduce that piece by piece
        let mut result: Vector<W, R::El> = Vector::new((self.rank_as_module()..=(self.rank_as_module() * 2 - 2)).map(|i|
            self.base_ring().sum(((i - self.rank_as_module() + 1)..self.rank_as_module()).map(|j| 
                self.base_ring().mul_ref(lhs.at(j), rhs.at(i - j))
            ))
        ).chain(std::iter::once(self.base_ring().zero())).collect());

        for i in (0..=(self.rank_as_module() - 2)).rev() {
            let value = std::mem::replace(result.at_mut(i), self.base_ring().zero());
            for j in 0..self.rank_as_module() {
                let index = if i + j < self.rank_as_module() {
                    i + j
                } else {
                    i + j - self.rank_as_module()
                };
                self.base_ring().add_assign(result.at_mut(index), self.base_ring().mul_ref(self.mipo_values.at(j), &value));
            }
        }

        for i in 0..self.rank_as_module() {
            for j in 0..=i {
                self.base_ring().add_assign(result.at_mut(i), self.base_ring().mul_ref(lhs.at(j), rhs.at(i - j)));
            }
        }

        return result;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        self.assert_valid_element(&val);
        val.scale(&self.base_ring.neg(self.base_ring.one()), &self.base_ring);
        return val;
    }

    fn zero(&self) -> Self::El {
        self.from(self.base_ring.zero())
    }

    fn one(&self) -> Self::El {
        self.from(self.base_ring.one())
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        lhs.as_ref().eq(rhs.as_ref(), &self.base_ring)
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        val.as_ref().eq(Vector::zero_ring(self.rank_as_module(), &self.base_ring), &self.base_ring)
    }

    fn is_one(&self, val: &Self::El) -> bool {
        self.base_ring.is_one(val.at(0)) && (
            self.rank_as_module() == 1 || 
            val.as_ref().subvector(1..).eq(Vector::zero_ring(self.rank_as_module() - 1, &self.base_ring), &self.base_ring)
        )
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        self.base_ring.is_neg_one(val.at(0)) && (
            self.rank_as_module() == 1 || 
            val.as_ref().subvector(1..).eq(Vector::zero_ring(self.rank_as_module() - 1, &self.base_ring), &self.base_ring)
        )
    }

    fn characteristic(&self) -> BigInt {
        self.base_ring().characteristic()
    }

    default fn is_integral(&self) -> RingPropValue {
        return RingPropValue::Unknown;
    }

    fn is_noetherian(&self) -> bool {
        self.base_ring.is_noetherian()
    }

    default fn is_field(&self) -> RingPropValue {
        return RingPropValue::Unknown;
    }

    fn div(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(&lhs);
        self.assert_valid_element(rhs);
        assert!(!self.is_zero(rhs));
        if self.base_ring.is_field().can_use() {
            let multiplication_matrix = self.create_multiplication_matrix(rhs.clone());
            <R as MatrixSolve>::solve_linear_equation(&self.base_ring, multiplication_matrix, &mut Matrix::col_vec(lhs.as_mut())).unwrap();
            return lhs;
        } else {
            unimplemented!()
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        let poly_ring = PolyRingImpl::adjoint(&self.base_ring, self.gen_name);
        poly_ring.format(&self.polynomial_repr(&poly_ring, el.clone()), f, in_prod)
    }
}

impl<R, V, W> RingExtension for FiniteExtensionImpl<R, V, W> 
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    type BaseRing = R;
    type Embedding = StandardEmbedding<R, Self>;

    fn is_extension(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn base_ring(&self) -> &R {
        &self.base_ring
    }

    fn from(&self, el: El<R>) -> El<Self> {
        Vector::new(
            std::iter::once(el)
            .chain(std::iter::repeat(self.base_ring.zero()).take(self.rank_as_module() - 1))
            .collect()
        )
    }

    fn embedding(&self) -> Self::Embedding {
        embedding(self.base_ring().clone(), self.clone())
    }
}

impl<R, V, W> FiniteExtension for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    type ModuleVectorView = W;

    fn generator(&self) -> El<Self> {
        Vector::new(
            std::iter::once(self.base_ring.zero())
            .chain(std::iter::once(self.base_ring.one()))
            .chain(std::iter::repeat(self.base_ring.zero()).take(self.rank_as_module() - 2))
            .collect()
        )
    }

    fn rank_as_module(&self) -> usize {
        self.mipo_values.len()
    }

    fn as_module_el(&self, x: El<Self>) -> Vector<Self::ModuleVectorView, El<Self::BaseRing>> {
        x
    }
}

impl<R, V, W> RingBase for FiniteExtensionImpl<R, V, W>
    where R: Ring, for<'a> PolyRingImpl<&'a R>: UfdInfoRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn is_integral(&self) -> RingPropValue {
        if !self.is_integral_cache.is_computed() {
            let poly_ring = PolyRingImpl::adjoint(&self.base_ring, "X");
            let can_compute = self.base_ring.is_integral() & poly_ring.is_ufd();
            if can_compute.can_use() {
                self.is_integral_cache.set(RingPropValue::True & poly_ring.is_prime(&self.generating_polynomial(&poly_ring)))
            } else {
                self.is_integral_cache.set(can_compute)
            }
        }
        self.is_integral_cache.get()
    }

    fn is_field(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() {
            self.is_integral()
        } else {
            self.base_ring.is_field()
        }
    }
}

impl<R, V, W> DivisibilityInfoRing for FiniteExtensionImpl<R, V, W>
    where R: Ring, for<'a> PolyRingImpl<&'a R>: UfdInfoRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        // currently only implemented in case we have a field
        self.is_field()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if self.is_zero(rhs) && self.is_zero(lhs) {
            Some(self.one())
        } else if self.is_zero(rhs) {
            None
        } else {
            Some(self.div(lhs.clone(), rhs))
        }
    }
}

impl<R, V, W> PartialEq for FiniteExtensionImpl<R, V, W>
    where R: Ring + PartialEq, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn eq(&self, rhs: &Self) -> bool {
        self.base_ring() == rhs.base_ring() && self.rank_as_module() == rhs.rank_as_module() && self.mipo_values.as_ref().eq(rhs.mipo_values.as_ref(), self.base_ring())
    }
}

impl<R, V, W> std::fmt::Debug for FiniteExtensionImpl<R, V, W>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Ring extension of {:?} generated by equation X^{} = ", &self.base_ring, self.rank_as_module())?;
        poly_format(&self.base_ring, self.mipo_values.as_ref(), f, "X")?;
        return Ok(());
    }
}

#[derive(Clone)]
pub struct FiniteRingExtensionIterFn<R, V, W> 
    where  R: FiniteRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    base_iter: MultiProduct<FiniteRingElementIter<R>, fn(&[El<R>]) -> Vector<W, El<R>>, El<FiniteExtensionImpl<R, V, W>>>,
    vector_type: PhantomData<V>
}

impl<R, V, W> FiniteRingIterFn<FiniteExtensionImpl<R, V, W>> for FiniteRingExtensionIterFn<R, V, W> 
    where  R: FiniteRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn next(&mut self, _: &FiniteExtensionImpl<R, V, W>) -> Option<El<FiniteExtensionImpl<R, V, W>>> {
        <_ as Iterator>::next(&mut self.base_iter)
    }
}

impl<R, V, W> FiniteRing for FiniteExtensionImpl<R, V, W>
    where  R: FiniteRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    type IterFn = FiniteRingExtensionIterFn<R, V, W>;

    fn size(&self) -> BigInt {
        BigInt::RING.pow(&self.base_ring().size(), self.rank_as_module() as u32)
    }

    fn iter_fn(&self) -> Self::IterFn {
        fn convert<R, W>(data: &[El<R>]) -> Vector<W, El<R>> 
            where R: Ring, W: VectorViewMut<R::El> + Clone + FromIterator<El<R>> + std::fmt::Debug
        {
            Vector::new(data.iter().cloned().collect())
        }
        FiniteRingExtensionIterFn{
            base_iter: multi_cartesian_product(
                std::iter::repeat(self.base_ring.clone()).map(|ring| finite_field_elements(ring)), 
                convert::<R, W>
            ),
            vector_type: PhantomData
        }
    }
    
    fn random_element<G>(&self, mut rng: G) -> El<Self> 
        where G: FnMut() -> u32
    {
        Vector::new((0..self.rank_as_module()).map(|_| self.base_ring().random_element(&mut rng)).collect())
    }
}

impl<R, V, W> HashableElRing for FiniteExtensionImpl<R, V, W>
    where  R: HashableElRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El) {
        for x in el.iter() {
            self.base_ring().hash(h, x);
        }
    }
}


#[cfg(test)]
use super::super::fq::zn_small::*;

#[test]
fn test_format() {
    let ring: FiniteExtensionImpl<_, _> = FiniteExtensionImpl::new(i64::RING, Vector::from_array([-1, 0]), "α");
    let i = ring.generator();
    assert_eq!("α", format!("{}", ring.display(&i)));
    assert_eq!("-1", format!("{}", ring.display(&ring.mul_ref(&i, &i))));
}

#[test]
fn test_division() {
    let base = ZnEl::<7>::RING;
    let field: FiniteExtensionImpl<_, _> = FiniteExtensionImpl::new(base, Vector::from_array([ZnEl::project(-1), ZnEl::project(0)]), "α");
    assert!(field.is_field().can_use());
    assert!(field.is_eq(&field.from(ZnEl::project(1) / ZnEl::project(2)), &field.div(field.from_z(1), &field.from_z(2))));

    let a = field.add(field.generator(), field.from_z(3));
    let b = field.add(field.mul(field.generator(), field.from_z(2)), field.from_z(4));
    let x = field.div(a.clone(), &b);
    assert!(field.is_eq(&a, &field.mul(b, x)));
}

#[test]
fn test_mul() {
    let base = i64::RING;
    let extension: FiniteExtensionImpl<_, _> = FiniteExtensionImpl::new(base, Vector::from_array([-1, 0]), "i");
    let a = extension.sub(extension.one(), extension.generator());
    let b = extension.add(extension.one(), extension.generator());
    assert!(extension.is_eq(&extension.from_z(2), &extension.mul(a, b)));
    assert!(extension.is_eq(&extension.from_z(-1), &extension.pow(&extension.generator(), 2)));
}

#[test]
fn test_adjoin_element() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i64::RING, "x"));
    let x = ring.unknown();
    let f = x.pow(3) + x.pow(2) + 1;
    let ring = WrappingRing::new(FiniteExtensionImpl::adjoin_element(i64::RING, f.val().clone(), f.parent_ring(), "x"));
    let x = ring.generator();
    assert_eq!(x.ring().zero(), f(x));
}
