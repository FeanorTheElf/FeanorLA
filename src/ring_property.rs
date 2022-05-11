use std::ops::BitAnd;
use std::sync::atomic::{AtomicU8, Ordering};
use std::default::Default;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingPropValue {
    True, False, Unknown
}

impl RingPropValue {

    pub fn can_use(&self) -> bool {
        *self == RingPropValue::True
    }
}

impl BitAnd for RingPropValue {

    type Output = RingPropValue;

    fn bitand(self, rhs: RingPropValue) -> RingPropValue {
        match (self, rhs) {
            (RingPropValue::False, _) => RingPropValue::False,
            (_, RingPropValue::False) => RingPropValue::False,
            (RingPropValue::Unknown, _) => RingPropValue::Unknown,
            (_, RingPropValue::Unknown) => RingPropValue::Unknown,
            (RingPropValue::True, RingPropValue::True) => RingPropValue::True
        }
    }
}

impl BitAnd<bool> for RingPropValue {

    type Output = RingPropValue;

    fn bitand(self, rhs: bool) -> RingPropValue {
        if rhs {
            self & RingPropValue::True
        } else {
            self & RingPropValue::False
        }
    }
}

impl std::fmt::Display for RingPropValue {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RingPropValue::True => write!(f, "true"),
            RingPropValue::False => write!(f, "false"),
            RingPropValue::Unknown => write!(f, "???")
        }
    }
}

#[derive(Debug)]
pub struct RingPropValueCache {
    content: AtomicU8
}

impl RingPropValueCache {

    pub const fn new() -> Self {
        RingPropValueCache {
            content: AtomicU8::new(0)
        }
    }

    pub fn is_computed(&self) -> bool {
        self.content.load(Ordering::Relaxed) != 0
    }

    pub fn get(&self) -> RingPropValue {
        match self.content.load(Ordering::Relaxed) {
            0 => panic!("Value not computed!"),
            1 => RingPropValue::False,
            2 => RingPropValue::True,
            3 => RingPropValue::Unknown,
            _ => unreachable!()
        }
    }

    pub fn set(&self, value: RingPropValue) {
        let new_value = match value {
            RingPropValue::False => 1,
            RingPropValue::True => 2,
            RingPropValue::Unknown => 3
        };
        let old_value = self.content.swap(new_value, Ordering::Relaxed);
        assert!(old_value == 0 || old_value == new_value);
    }
}

impl Clone for RingPropValueCache {

    fn clone(&self) -> Self {
        RingPropValueCache {
            content: AtomicU8::new(self.content.load(Ordering::Relaxed))
        }
    }
}

impl Default for RingPropValueCache {

    fn default() -> Self {
        Self::new()
    }
}