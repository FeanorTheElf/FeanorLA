use super::super::float::*;
use super::super::la::mat::*;

type Point<T> = Vector<VectorOwned<T>, T>;
type PointRef<'a, T> = Vector<&'a VectorOwned<T>, T>;

#[derive(Debug)]
pub struct NotConvergent<T>(pub Point<T>);

///
/// Applies the Newton-method to find a root of the given function
/// that maps from R^k to R^k.
/// 
/// The parameter `jacobian_solve` should on input x and y calculate a 
/// (potentially least-squares) solution of the linear equation
/// system J_f(x) z = y where J_f(x) is the Jacobian matrix of f at 
/// the point x.
/// 
pub fn newton_multidim<T, F, S>(
    mut function: F,
    mut jacobian_solve: S, 
    starting_point: Point<T>, 
    error: T, 
    max_iterations: usize
) -> Result<Point<T>, NotConvergent<T>>
    where T: Float,
        F: FnMut(PointRef<T>) -> Point<T>,
        S: FnMut(PointRef<T>, Point<T>) -> Point<T>
{
    let mut current = starting_point.to_owned();
    for _ in 0..max_iterations {
        let value = function(current.as_ref());
        if value.approx_eq(&Vector::zero(value.len()), &error) {
            return Ok(current.to_owned());
        }
        current = current.as_ref() - jacobian_solve(current.as_ref(), value);
    }
    return Err(NotConvergent(current));
}

#[cfg(test)]
use super::super::alg::*;
#[cfg(test)]
use super::super::algebra::poly::*;
#[cfg(test)]
use super::super::la::algorithms::*;

#[test]
fn test_newton_multidim() {
    // the function should be f(x, y) = (xy + 2, x + 2 * y)
    let mut ring = MultivariatePolyRing::new(&f64::RING);
    let x = ring.adjoint("x");
    let y = ring.adjoint("y");
    let f1 = ring.add(ring.mul_ref(&x, &y), ring.from(2.));
    let f2 = ring.add(x.clone(), ring.mul_ref(&y, &ring.from(2.)));

    #[rustfmt::skip]
    let jacobian = Matrix::from_array([[ring.derive(&f1, ring.get_var("x")), ring.derive(&f1, ring.get_var("y"))],
                                       [ring.derive(&f2, ring.get_var("x")), ring.derive(&f2, ring.get_var("y"))]]);
                                       
    let function = Vector::from_array([f1, f2]);
    println!("Function: {}", Matrix::col_vec(function.as_ref()).display(&ring));
    println!("Jacobian: {}", jacobian.display(&ring));
    
    let function = |point: PointRef<f64>| {
        ring.evaluate_vector_at(function.as_ref(), &point)
    };
    let jacobian_solve = |point: PointRef<f64>, mut y: Point<f64>| {
        f64::RING.solve_linear_equation(
            ring.evaluate_matrix_at(jacobian.as_ref(), &point), 
            &mut Matrix::col_vec(y.as_mut())
        ).unwrap();
        return y;
    };
    let root = newton_multidim(function, jacobian_solve, Vector::from_array([10., 10.]), 0.0001, 30).unwrap();
    assert_approx_eq!(Vector::from_array([2., -1.]), &root, 0.001);
}