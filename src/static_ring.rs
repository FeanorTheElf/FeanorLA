
#[macro_export]
macro_rules! static_ring_ref {
    ($name:ident, $static_ring:ident: $static_ring_type:ty) => {
        pub struct $name;

        impl RingDecorator for $name {

            type DecoratedRing = $static_ring_type;

            fn decorated_ring(&self) -> &Self::DecoratedRing {
                &$static_ring
            }
        }
    };
}