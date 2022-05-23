use super::prelude::*;

fn add_assign_intersect<R, U, V>(mut dst: Vector<U, R::El>, val: Vector<V, R::El>, ring: &R)
    where U: VectorViewMut<R::El>, V: VectorView<R::El>, R: Ring
{
    let l = std::cmp::min(dst.len(), val.len());
    dst.subvector_mut(..l).add_assign(val.subvector(..l), ring);
}

fn sub_assign_intersect<R, U, V>(mut dst: Vector<U, R::El>, val: Vector<V, R::El>, ring: &R)
    where U: VectorViewMut<R::El>, V: VectorView<R::El>, R: Ring
{
    let l = std::cmp::min(dst.len(), val.len());
    dst.subvector_mut(..l).sub_assign(val.subvector(..l), ring);
}

fn karatsuba_impl<R, U, V, W>(
    lhs: Vector<Subvector<U, R::El>, R::El>,
    rhs: Vector<Subvector<V, R::El>, R::El>,
    mut out: Vector<Subvector<W, R::El>, R::El>,
    tmp: &mut [R::El],
    ring: &R,
    block_size: usize
)
    where R: Ring, U: VectorView<R::El> + Copy, V: VectorView<R::El> + Copy, W: VectorViewMut<R::El>
{
    assert!(tmp.len() == 4 * block_size - 4);
    assert!(lhs.len() <= block_size);
    assert!(rhs.len() <= block_size);

    if lhs.len() == 0 || rhs.len() == 0 {
        return;
    }

    if block_size == 1 {
        ring.add_assign(out.at_mut(0), ring.mul_ref(lhs.at(0), rhs.at(0)));
        return;
    }

    let block_half = block_size/2;

    let (lhs_combined, rest) = tmp.split_at_mut(block_half);
    let (rhs_combined, rest) = rest.split_at_mut(block_half);
    let (mul_result, rest) = rest.split_at_mut(block_size);
    let mut lhs_combined = Vector::new(lhs_combined);
    let mut rhs_combined = Vector::new(rhs_combined);
    let mut mul_result = Vector::new(mul_result);

    mul_result.assign(Vector::zero_ring(block_size, ring));

    karatsuba_impl(
        lhs.into_subvector_intersect(..block_half), 
        rhs.into_subvector_intersect(..block_half), 
        mul_result.subvector_mut(..), 
        rest, 
        ring, 
        block_half
    );

    add_assign_intersect(out.subvector_mut_intersect(..block_size), mul_result.as_ref(), ring);
    sub_assign_intersect(out.subvector_mut_intersect(block_half..(block_size + block_half)), mul_result.as_ref(), ring);

    mul_result.assign(Vector::zero_ring(block_size, ring));

    karatsuba_impl(
        lhs.into_subvector_intersect(block_half..), 
        rhs.into_subvector_intersect(block_half..), 
        mul_result.subvector_mut(..), 
        rest, 
        ring, 
        block_half
    );

    add_assign_intersect(out.subvector_mut_intersect(block_size..), mul_result.as_ref(), ring);
    sub_assign_intersect(out.subvector_mut_intersect(block_half..(block_size + block_half)), mul_result.as_ref(), ring);

    lhs_combined.assign(lhs.into_subvector_intersect(..block_half));
    add_assign_intersect(lhs_combined.as_mut(), lhs.into_subvector_intersect(block_half..), ring);

    rhs_combined.assign(rhs.into_subvector_intersect(..block_half));
    add_assign_intersect(rhs_combined.as_mut(), rhs.into_subvector_intersect(block_half..), ring);

    karatsuba_impl(
        lhs_combined.subvector(..), 
        rhs_combined.subvector(..),
        out.into_subvector_mut_intersect(block_half..), 
        rest, 
        ring, 
        block_half
    );

}

#[test]
fn test_karatsuba_impl() {
    let a = Vector::from_array([1, 2, 3]);
    let b = Vector::from_array([3, 4, 5]);
    let mut c = Vector::from_array([0, 0, 0, 0, 0]);
    let mut tmp = (0..12).collect::<Vec<_>>();
    karatsuba_impl(a.subvector(..), b.subvector(..), c.subvector_mut(..), &mut tmp[..], &i64::RING, 4);
    assert_eq!(Vector::from_array([3, 10, 22, 22, 15]), c);
}