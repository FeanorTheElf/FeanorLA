#![allow(unused)]
use super::ring::*;
use super::primitive::*;
use super::bigint::*;
use super::wrapper::*;
use super::embedding::*;
use super::la::mat::*;
use super::algebra::poly::*;
use super::algebra::rat::*;
use super::algebra::eea::*;
use super::algebra::ring_ext::*;
use super::algebra::fractions::*;
use super::algebra::primes::*;
use super::algebra::fq::*;
use super::algebra::elliptic_curve::*;
use super::combinatorics::iters::*;

type QType = WrappingRing<FieldOfFractions<BigIntRing>>;
type ZType = WrappingRing<BigIntRing>;
type ZEl = <ZType as Ring>::El;
type FType = <zn_small::ZnEl<7> as RingEl>::RingType;
