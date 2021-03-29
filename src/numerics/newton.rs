use super::super::float::*;
use super::super::la::mat::*;

type Point<'a, T> = Vector<VectorRef<'a, VectorOwned<T>, T>, T>;

pub struct NotConvergent;

///
/// The parameter `jacobian_solve` should on input x calculate a 
/// (potentially least-squares) solution of the linear equation
/// system J_f(x) y = f(x) where J_f(x) is the Jacobian matrix
/// of f at the point x.
/// 
pub fn newton_base<T, S>(mut jacobian_solve: S, starting_point: Point<T>, error: T, max_iterations: usize) -> Result<Point<T>, NotConvergent>
    where T: Float,
        S: FnMut(Point<T>) -> Point<T>
{
    let mut current = starting_point.to_owned();
    for _ in 0..max_iterations {
        current = current.as_ref() - jacobian_solve(current.as_ref());
    }
    return Err(NotConvergent);
}