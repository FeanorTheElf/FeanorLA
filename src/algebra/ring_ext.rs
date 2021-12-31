use super::super::alg::*;
use super::poly::*;
use super::super::la::vec::*;
use super::super::la::vector_view::*;

#[derive(Debug, Clone)]
pub struct SimpleRingExtension<R>
    where R: Ring
{
    base_ring: R,
    mipo_values: Box<[R::El]>
}

impl<R> SimpleRingExtension<R>
    where R: Ring
{
    ///
    /// Creates the ring R[X]/(f) where 
    /// ```
    /// f = mipo_values[0] + mipo_values[1] * X + ... + mipo_values[n - 1] * X^(n - 1) - X^n
    /// ```
    /// 
    pub fn new(base_ring: R, mipo_values: Box<[R::El]>) -> Self {
        SimpleRingExtension {
            base_ring, mipo_values
        }
    }

    pub fn adjoin_element<F>(base_ring: R, mipo: F) -> Self
        where F: FnOnce(&PolyRing<&R>) -> <PolyRing<R> as Ring>::El
    {
        let ring = PolyRing::adjoint(&base_ring, "X");
        let mut mipo = mipo(&ring);
        let scaling_factor = base_ring.neg(base_ring.div(base_ring.one(), ring.lc(&mipo).expect("Minimal polynomial must not be constant when creating a new ring")));
        mipo = ring.mul(mipo, ring.from(scaling_factor));
        let degree = mipo.iter().enumerate().filter(|(_, x)| !base_ring.is_zero(x)).map(|(i, _)| i).max().unwrap();
        let mipo_values = mipo.into_iter().take(degree).collect::<Vec<_>>().into_boxed_slice();
        return Self::new(base_ring, mipo_values);
    }

    pub fn degree(&self) -> usize {
        self.mipo_values.len()
    }

    fn assert_valid_element(&self, el: &<Self as Ring>::El) {
        assert_eq!(el.len(), self.degree());
    }

    pub fn polynomial_repr(&self, poly_ring: &PolyRing<&R>, el: <Self as Ring>::El) -> <PolyRing<R> as Ring>::El {
        el.raw_data().into_vec().into_iter().enumerate()
            .map(|(i, x)| poly_ring.mul(poly_ring.pow(&poly_ring.unknown(), i as u32), poly_ring.from(x)))
            .fold(poly_ring.zero(), |a, b| poly_ring.add(a, b))
    }

    pub fn generator(&self) -> <Self as Ring>::El {
        Vector::unit_vector_ring(1, self.degree(), &self.base_ring)
    }
}

impl<R> Ring for SimpleRingExtension<R>
    where R: Ring
{
    type El = Vector<VectorOwned<R::El>, R::El>;

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
                result[i + j] = self.base_ring.mul_ref(&lhs[i], &rhs[j]);
            }
        }
        for i in (self.degree()..result.len()).rev() {
            let x = result.pop().unwrap();
            for j in (i - self.degree())..i {
                take_mut::take(&mut result[j], |y| self.base_ring.add(y, 
                    self.base_ring.mul_ref(&self.mipo_values[j + self.degree() - i], &x)
                ));
            }
        }
        return Vector::new(VectorOwned::new(result.into_boxed_slice()));
    }

    fn neg(&self, val: Self::El) -> Self::El {
        self.assert_valid_element(&val);
        val.neg(&self.base_ring)
    }

    fn zero(&self) -> Self::El {
        Vector::zero_ring(self.degree(), &self.base_ring).into_owned()
    }

    fn one(&self) -> Self::El {
        Vector::unit_vector_ring(0, self.degree(), &self.base_ring)
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

    default fn is_integral(&self) -> bool {
        // Cannot check whether a simple ring extension is integral if the corresponding polynomial ring does not support primality information
        return false;
    }

    default fn is_euclidean(&self) -> bool {
        // Cannot check whether a simple ring extension is euclidean if the corresponding polynomial ring does not support primality information
        return false;
    }

    default fn is_field(&self) -> bool {
        // Cannot check whether a simple ring extension is a field if the corresponding polynomial ring does not support primality information
        return false
    }

    default fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        panic!("Not a euclidean domain")
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid_element(&lhs);
        self.assert_valid_element(rhs);
        assert!(!self.is_zero(rhs));
        unimplemented!()
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        let poly_ring = PolyRing::adjoint(&self.base_ring, "α");
        poly_ring.format(&self.polynomial_repr(&poly_ring, el.clone()), f, in_prod)
    }
}

#[test]
fn test_format() {
    let ring = SimpleRingExtension::new(i64::RING, vec![-1, 0].into_boxed_slice());
    let i = ring.generator();
    assert_eq!("α", format!("{}", ring.display(&i)));
    assert_eq!("-1", format!("{}", ring.display(&ring.mul_ref(&i, &i))));
}