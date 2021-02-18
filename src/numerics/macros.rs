use super::super::la::mat::*;

pub trait ApproxEq {
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool;
}

impl ApproxEq for f64 {
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool {
        (self - rhs).abs() < delta
    }
}

impl<M, T> ApproxEq for Matrix<M, T>
    where M: MatrixView<T>, T: ApproxEq
{
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool {
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                if !self.at(row, col).approx_eq(rhs.at(row, col), delta) {
                    return false;
                }
            }
        }
        return true;
    }
}

impl<V, T> ApproxEq for Vector<V, T>
    where V: VectorView<T>, T: ApproxEq
{
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool {
        for i in 0..self.len() {
            if !self.at(i).approx_eq(rhs.at(i), delta) {
                return false;
            }
        }
        return true;
    }
}

macro_rules! assert_approx_eq {
    ($left:expr, $right:expr, $delta:expr) => {
        if !($left).approx_eq(($right), ($delta))
        {
            panic!(
                r#"assertion failed: `(left == right +-{})`
  left: `{:?}`,
 right: `{:?}`"#,
                $delta, $left, $right
            );
        }
    };
}
