
pub mod define_fq {
    pub use crate::ring_extension::simple_extension::*;
    pub use crate::fq::zn_small::*;
    pub use crate::la::vector_view::compile_time_vector::*;
    pub use crate::la::vec::*;
    pub use crate::prelude::*;
}

mod internal_definitions {
    use super::define_fq::*;

    // because of a (probable) compiler bug, `ZnEl<2>` does not work here
    type F2El = ZnElImpl<2, true>;
    type F7El = ZnElImpl<7, true>;
    type F37El = ZnElImpl<37, true>;

    gen_const_vector!(ConstVector2F2; F2El; V0, V1);
    pub type F4MipoType = ConstVector2F2<{F2El::project(-1)}, {F2El::project(-1)}>;
    pub const F4_MIPO: Vector<F4MipoType, F2El> = Vector::new(F4MipoType {});

    gen_const_vector!(ConstVector2F7; F7El; V0, V1);
    pub type F49MipoType = ConstVector2F7<{F7El::project(-1)}, {F7El::project(0)}>;
    pub const F49_MIPO: Vector<F49MipoType, F7El> = Vector::new(F49MipoType {});

    gen_const_vector!(ConstVector2F37; F37El; V0, V1);
    pub type F1369MipoType = ConstVector2F37<{F37El::project(2)}, {F37El::project(33)}>;
    pub const F1369_MIPO: Vector<F1369MipoType, F37El> = Vector::new(F1369MipoType {});
}

use internal_definitions::*;
use define_fq::*;

pub type F2Type = StaticRing::<ZnEl<2>>;
pub type F3Type = StaticRing::<ZnEl<3>>;
pub type F5Type = StaticRing::<ZnEl<5>>;
pub type F7Type = StaticRing::<ZnEl<7>>;
pub const F2: F2Type = ZnEl::<2>::RING;
pub const F3: F3Type = ZnEl::<3>::RING;
pub const F5: F5Type = ZnEl::<5>::RING;
pub const F7: F7Type = ZnEl::<7>::RING;

pub type F4Type = SimpleRingExtension<StaticRing::<ZnEl<2>>, F4MipoType, VectorArray<ZnEl<2>, 2>>;
pub type F49Type = SimpleRingExtension<StaticRing::<ZnEl<7>>, F49MipoType, VectorArray<ZnEl<7>, 2>>;
pub type F1369Type = SimpleRingExtension<StaticRing<ZnEl<37>>, F1369MipoType, VectorArray<ZnEl<37>, 2>>;
pub const F4: F4Type = SimpleRingExtension::new(F2, F4_MIPO);
pub const F49: F49Type = SimpleRingExtension::new(F7, F49_MIPO);
pub const F1369: F1369Type = F1369Type::new(ZnEl::<37>::RING, F1369_MIPO);

#[cfg(test)]
use super::FiniteRing;

#[test]
fn test_arithmetic() {
    let a = F4.generator();
    assert!(!F4.eq(&a, &F4.one()));
    assert!(F4.eq(&F4.pow(&a, 3), &F4.one()));
    assert!(F4.eq(&F4.pow(&F4.add(a, F4.one()), 3), &F4.one()));
}

#[test]
fn test_division() {
    let a = F4.generator();
    let b = F4.add_ref(F4.one(), &a);
    let c = F4.div(F4.one(), &a);
    assert!(F4.eq(&b, &c));
}

#[test]
fn test_finite_odd_field_is_field() {
    assert!(F49.is_field().can_use());
}

#[test]
fn test_is_ring() {
    fn as_ring<R: Ring>(_: R) {
        assert!(true);
    }
    as_ring(F49);
}

#[test]
fn test_is_finite_ring() {
    fn as_ring<R: FiniteRing>(_: R) {
        assert!(true);
    }
    as_ring(F49);
}

#[test]
#[ignore] // currently, we cannot factor polynomials over Z/2Z
fn test_finite_even_field_is_field() {
    assert_eq!(RingPropValue::True, F4.is_field());
}