pub use crate::float::ApproxEq;

macro_rules! assert_approx_eq {
    ($left:expr, $right:expr, $delta:expr) => {
        if !($left).approx_eq(($right), &($delta))
        {
            panic!(
                r#"assertion failed: `(left == right, allowed error: {})`
  left: `{:?}`,
 right: `{:?}`"#,
                $delta, $left, $right
            );
        }
    };
}
