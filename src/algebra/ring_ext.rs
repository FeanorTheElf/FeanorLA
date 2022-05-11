use super::super::prelude::*;
use super::poly::*;
use super::super::la::mat::*;
use super::super::la::algorithms::*;
use super::poly::ops::poly_format;

use std::marker::PhantomData;
use std::iter::FromIterator;

#[derive(Clone)]
pub struct SimpleRingExtension<R, V, W = VectorOwned<<R as RingBase>::El>>
    where R: Ring, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El>
{
    base_ring: R,
    mipo_values: Vector<V, R::El>,
    element_vector: PhantomData<W>,
    is_integral_cache: RingPropValueCache
}

impl<R, W> SimpleRingExtension<R, VectorOwned<R::El>, W>
    where R: CanonicalIsomorphismInfo<R>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    pub fn adjoin_element<F>(base_ring: R, mipo: F) -> Self
        where F: FnOnce(&PolyRing<&R>) -> <PolyRing<R> as RingBase>::El
    {
        let ring = PolyRing::adjoint(&base_ring, "X");
        let mut mipo = mipo(&ring);
        let scaling_factor = base_ring.neg(base_ring.div(base_ring.one(), ring.lc(&mipo).expect("Minimal polynomial must not be constant when creating a new ring")));
        mipo = ring.mul(mipo, ring.from(scaling_factor));
        let degree = mipo.iter().enumerate().filter(|(_, x)| !base_ring.is_zero(x)).map(|(i, _)| i).max().unwrap();
        let mipo_values = Vector::new(mipo.raw_data().into_vec().into_iter().take(degree).collect::<Vec<_>>());
        return Self::new(base_ring, mipo_values);
    }
}

impl<R, V, W> SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    ///
    /// Creates the ring R[X]/(f) where 
    /// ```text
    /// f = mipo_values[0] + mipo_values[1] * X + ... + mipo_values[n - 1] * X^(n - 1) - X^n
    /// ```
    /// 
    pub const fn new(base_ring: R, mipo_values: Vector<V, R::El>) -> Self {
        SimpleRingExtension {
            base_ring: base_ring, 
            mipo_values: mipo_values,
            element_vector: PhantomData,
            is_integral_cache: RingPropValueCache::new()
        }
    }

    pub fn degree(&self) -> usize {
        self.mipo_values.len()
    }

    fn assert_valid_element(&self, el: &<Self as RingBase>::El) {
        assert_eq!(el.len(), self.degree());
    }

    pub fn polynomial_repr(&self, poly_ring: &PolyRing<&R>, el: <Self as RingBase>::El) -> El<PolyRing<R>> {
        el.into_owned().raw_data().into_vec().into_iter().enumerate()
            .map(|(i, x)| poly_ring.mul(poly_ring.pow(&poly_ring.unknown(), i as u32), poly_ring.from(x)))
            .fold(poly_ring.zero(), |a, b| poly_ring.add(a, b))
    }

    pub fn generator(&self) -> El<Self> {
        Vector::new(
            std::iter::once(self.base_ring.zero())
            .chain(std::iter::once(self.base_ring.one()))
            .chain(std::iter::repeat(self.base_ring.zero()).take(self.degree() - 2))
            .collect()
        )
    }

    pub fn from(&self, el: R::El) -> El<Self> {
        Vector::new(
            std::iter::once(el)
            .chain(std::iter::repeat(self.base_ring.zero()).take(self.degree() - 1))
            .collect()
        )
    }

    fn create_multiplication_matrix(&self, el: El<Self>) -> Matrix<MatrixOwned<R::El>, R::El> {
        let d = self.degree();
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
                *matrix.at_mut(j, i) =self.base_ring.mul_ref(&self.mipo_values[j], &last_el);
            }
        }
        return matrix;
    }

    ///
    /// Returns the (monic) polynomial f such that this ring is isomorphic to
    /// `base_ring[X] / (f)`.
    /// 
    pub fn generating_polynomial(&self, poly_ring: &PolyRing<&R>) -> El<PolyRing<&R>> {
        self.mipo_values.iter()
            .cloned()
            .chain(std::iter::once(self.base_ring.neg(self.base_ring.one())))
            .scan(poly_ring.one(), |state, coeff| {
                let result = poly_ring.mul_ref(state, &poly_ring.from(coeff));
                take_mut::take_or_recover(state, || poly_ring.unspecified_element(), 
                    |x| poly_ring.mul(x, poly_ring.unknown())
                );
                return Some(result);
            }).fold(poly_ring.zero(), |a, b| poly_ring.add(a, b))
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<R> for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone + FromIterator<R::El>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, _from: &R) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &R, el: R::El) -> Self::El {
        self.from(el)
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<&R> for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone + FromIterator<R::El>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, _from: &&R) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &&R, el: R::El) -> Self::El {
        self.from(el)
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<R> for SimpleRingExtension<&R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone + FromIterator<R::El>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, _from: &R) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &R, el: R::El) -> Self::El {
        self.from(el)
    }
}

impl<R, V, W> CanonicalEmbeddingInfo<SimpleRingExtension<R, V, W>> for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone + FromIterator<R::El>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_embedding(&self, _from: &Self) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &Self, el: Self::El) -> Self::El {
        el
    }
}

impl<R, V, W> CanonicalIsomorphismInfo<SimpleRingExtension<R, V, W>> for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone + FromIterator<R::El>, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn has_isomorphism(&self, _from: &Self) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &Self, el: Self::El) -> Self::El {
        el
    }
}

impl<R, V, W> RingBase for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    type El = Vector<W, R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(&lhs);
        self.assert_valid_element(rhs);

        lhs.add_assign(rhs.as_ref(), &self.base_ring);
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(&lhs);
        self.assert_valid_element(rhs);

        let mut result = (0..(2 * self.degree())).map(|_| self.base_ring.zero()).collect::<Vec<_>>();
        for i in 0..lhs.len() {
            for j in 0..rhs.len() {
                take_mut::take_or_recover(&mut result[i + j], || self.base_ring.unspecified_element(), |y| self.base_ring.add(y, 
                    self.base_ring.mul_ref(&lhs[i], &rhs[j])
                ));
            }
        }
        for i in (self.degree()..result.len()).rev() {
            let x = result.pop().unwrap();
            for j in (i - self.degree())..i {
                take_mut::take_or_recover(&mut result[j], || self.base_ring.unspecified_element(), |y| self.base_ring.add(y, 
                    self.base_ring.mul_ref(&self.mipo_values[j + self.degree() - i], &x)
                ));
            }
        }
        return Vector::new(result.into_iter().collect());
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

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        lhs.as_ref().eq(rhs.as_ref(), &self.base_ring)
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        val.as_ref().eq(Vector::zero_ring(self.degree(), &self.base_ring), &self.base_ring)
    }

    fn is_one(&self, val: &Self::El) -> bool {
        self.base_ring.is_one(val.at(0)) && (
            self.degree() == 1 || 
            val.as_ref().subvector(1..).eq(Vector::zero_ring(self.degree() - 1, &self.base_ring), &self.base_ring)
        )
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        self.base_ring.is_neg_one(val.at(0)) && (
            self.degree() == 1 || 
            val.as_ref().subvector(1..).eq(Vector::zero_ring(self.degree() - 1, &self.base_ring), &self.base_ring)
        )
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
        let poly_ring = PolyRing::adjoint(&self.base_ring, "α");
        poly_ring.format(&self.polynomial_repr(&poly_ring, el.clone()), f, in_prod)
    }
}

impl<R, V, W> RingBase for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, for<'a> PolyRing<&'a R>: UfdInfoRing, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn is_integral(&self) -> RingPropValue {
        if !self.is_integral_cache.is_computed() {
            let poly_ring = PolyRing::adjoint(&self.base_ring, "X");
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

impl<R, V, W> std::fmt::Debug for SimpleRingExtension<R, V, W>
    where R: CanonicalIsomorphismInfo<R>, V: VectorView<R::El> + Clone, W: VectorViewMut<R::El> + Clone + FromIterator<R::El> + std::fmt::Debug
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Ring extension of {:?} generated by equation X^{} = ", &self.base_ring, self.degree())?;
        poly_format(&self.base_ring, self.mipo_values.as_ref(), f, "X")?;
        return Ok(());
    }
}

#[test]
fn test_format() {
    let ring: SimpleRingExtension<_, _> = SimpleRingExtension::new(i64::RING, Vector::from_array([-1, 0]));
    let i = ring.generator();
    assert_eq!("α", format!("{}", ring.display(&i)));
    assert_eq!("-1", format!("{}", ring.display(&ring.mul_ref(&i, &i))));
}
