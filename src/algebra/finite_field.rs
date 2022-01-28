
mod internal_definitions{
    pub use super::super::ring_ext::*;
    pub use super::super::zn::*;
    pub use super::super::super::la::vec::*;
    pub use super::super::super::la::const_vector::*;
    pub use super::super::super::alg::*;

    // because of a (probable) compiler bug, `ZnEl<2>` does not work here
    type F2El = ZnElImpl<2, true>;
    
    gen_const_vector!(ConstVector2F2; F2El; V0, V1);
    pub type F4MipoType = ConstVector2F2<{F2El::project(-1)}, {F2El::project(-1)}>;
    pub const F4_MIPO: Vector<F4MipoType, F2El> = Vector::new(F4MipoType {});
}

use internal_definitions::*;

pub const F2: StaticRing::<ZnEl<2>> = ZnEl::<2>::RING;
pub const F3: StaticRing::<ZnEl<3>> = ZnEl::<3>::RING;
pub const F5: StaticRing::<ZnEl<5>> = ZnEl::<5>::RING;
pub const F7: StaticRing::<ZnEl<7>> = ZnEl::<7>::RING;
pub const F4: SimpleRingExtension<StaticRing::<ZnEl<2>>, F4MipoType> = SimpleRingExtension::new(F2, F4_MIPO);

#[test]
fn test_zero_size() {
    assert_eq!(0, std::mem::size_of_val(&F4));
}

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
    println!("Found gen");
    let b = F4.add_ref(F4.one(), &a);
    println!("Found b");
    let c = F4.div(F4.one(), &a);
    println!("Found c");
    assert!(F4.eq(&b, &c));
}
