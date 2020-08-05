
#[cfg(test)]
fn approx_eq(lhs: &[f64], rhs: &[f64], delta: f64) -> bool {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(a, b)| a - b)
        .map(f64::abs)
        .all(|a| a < delta)
}

#[cfg(test)]
macro_rules! assert_approx_eq {
    ($left:expr, $right:expr, $delta:expr) => {
        if !($left)
            .iter()
            .zip(($right).iter())
            .map(|(a, b)| a - b)
            .map(f64::abs)
            .all(|a| a < $delta)
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
