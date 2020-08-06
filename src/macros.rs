use super::matrix::MatrixView;

#[cfg(test)]
pub trait ApproxEq {
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool;
}

#[cfg(test)]
impl ApproxEq for f64 {
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool {
        (self - rhs).abs() < delta
    }
}

#[cfg(test)]
impl<M> ApproxEq for M
    where M: MatrixView<f64>
{
    fn approx_eq(&self, rhs: &Self, delta: f64) -> bool {
        for row in 0..self.rows() {
            for col in 0..self.cols() {
                if !self.at(row, col).approx_eq(rhs.at(row, col), delta) {
                    return false;
                }
            }
        }
        return true;
    }
}

#[cfg(test)]
macro_rules! assert_approx_eq {
    ($left:expr, $right:expr, $delta:expr) => {
        if !($left).approx_eq(($right), ($delta))
        {
            panic!(
                r#"assertion failed: `(left == right +-{:?})`
  left: `{:?}`,
 right: `{:?}`"#,
                $delta, $left, $right
            );
        }
    };
}
