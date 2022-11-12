
use super::fq_small;
use super::zn_small::*;
use super::super::finite_extension::finite_extension_impl::*;
use super::super::la::vec::*;
use super::super::prelude::*;

mod internal {

    use super::super::fq_small::define_fq::*;

    type F2El = ZnElImpl<2, true>;

    gen_const_vector!(ConstVector3F2; F2El; V0, V1, V2);
    pub type F8MipoType = ConstVector3F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}>;
    pub const F8_MIPO: Vector<F8MipoType, F2El> = Vector::new(F8MipoType {});

    gen_const_vector!(ConstVector4F2; F2El; V0, V1, V2, V3);
    pub type F16MipoType = ConstVector4F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F16_MIPO: Vector<F16MipoType, F2El> = Vector::new(F16MipoType {});

    gen_const_vector!(ConstVector5F2; F2El; V0, V1, V2, V3, V4);
    pub type F32MipoType = ConstVector5F2<{F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F32_MIPO: Vector<F32MipoType, F2El> = Vector::new(F32MipoType {});

    gen_const_vector!(ConstVector6F2; F2El; V0, V1, V2, V3, V4, V5);
    pub type F64MipoType = ConstVector6F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}>;
    pub const F64_MIPO: Vector<F64MipoType, F2El> = Vector::new(F64MipoType {});

    gen_const_vector!(ConstVector7F2; F2El; V0, V1, V2, V3, V4, V5, V6);
    pub type F128MipoType = ConstVector7F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F128_MIPO: Vector<F128MipoType, F2El> = Vector::new(F128MipoType {});

    gen_const_vector!(ConstVector8F2; F2El; V0, V1, V2, V3, V4, V5, V6, V7);
    pub type F256MipoType = ConstVector8F2<{F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F256_MIPO: Vector<F256MipoType, F2El> = Vector::new(F256MipoType {});

    gen_const_vector!(ConstVector9F2; F2El; V0, V1, V2, V3, V4, V5, V6, V7, V8);
    pub type F512MipoType = ConstVector9F2<{F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F512_MIPO: Vector<F512MipoType, F2El> = Vector::new(F512MipoType {});

    gen_const_vector!(ConstVector10F2; F2El; V0, V1, V2, V3, V4, V5, V6, V7, V8, V9);
    pub type F1024MipoType = ConstVector10F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F1024_MIPO: Vector<F1024MipoType, F2El> = Vector::new(F1024MipoType {});

    gen_const_vector!(ConstVector12F2; F2El; V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11);
    pub type F4096MipoType = ConstVector12F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F4096_MIPO: Vector<F4096MipoType, F2El> = Vector::new(F4096MipoType {});

    gen_const_vector!(ConstVector16F2; F2El; V0, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15);
    pub type F65536MipoType = ConstVector16F2<{F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(-1)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}, {F2El::project(0)}>;
    pub const F65536_MIPO: Vector<F65536MipoType, F2El> = Vector::new(F65536MipoType {});
}

pub type F2Type = fq_small::F2Type;
pub const F2: F2Type = fq_small::F2;

pub type F4Type = fq_small::F4Type;
pub const F4: F4Type = fq_small::F4;

pub type F8Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F8MipoType, VectorArray<ZnEl<2>, 3>>;
pub const F8: F8Type = F8Type::new(F2, internal::F8_MIPO, "α");

pub type F16Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F16MipoType, VectorArray<ZnEl<2>, 4>>;
pub const F16: F16Type = F16Type::new(F2, internal::F16_MIPO, "α");

pub type F32Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F32MipoType, VectorArray<ZnEl<2>, 5>>;
pub const F32: F32Type = F32Type::new(F2, internal::F32_MIPO, "α");

pub type F64Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F64MipoType, VectorArray<ZnEl<2>, 6>>;
pub const F64: F64Type = F64Type::new(F2, internal::F64_MIPO, "α");

pub type F128Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F128MipoType, VectorArray<ZnEl<2>, 7>>;
pub const F128: F128Type = F128Type::new(F2, internal::F128_MIPO, "α");

pub type F256Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F256MipoType, VectorArray<ZnEl<2>, 8>>;
pub const F256: F256Type = F256Type::new(F2, internal::F256_MIPO, "α");

pub type F512Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F512MipoType, VectorArray<ZnEl<2>, 9>>;
pub const F512: F512Type = F512Type::new(F2, internal::F512_MIPO, "α");

pub type F1024Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F1024MipoType, VectorArray<ZnEl<2>, 10>>;
pub const F1024: F1024Type = F1024Type::new(F2, internal::F1024_MIPO, "α");

pub type F4096Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F4096MipoType, VectorArray<ZnEl<2>, 12>>;
pub const F4096: F4096Type = F4096Type::new(F2, internal::F4096_MIPO, "α");

pub type F65536Type = FiniteExtensionImpl<StaticRing::<ZnEl<2>>, internal::F65536MipoType, VectorArray<ZnEl<2>, 16>>;
pub const F65536: F65536Type = F65536Type::new(F2, internal::F65536_MIPO, "α");
