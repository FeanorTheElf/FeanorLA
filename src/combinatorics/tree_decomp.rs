use super::graph::*;
use super::flow::*;

use std::collections::HashSet;

///
/// Tries to calculate a tree decomposition of the graph of size at most 4k + 3. 
/// If no result is returned, then the graph does not allow a tree decomposition of
/// size at most k.
/// 
pub fn calc_tree_decomp<G>(graph: &G, k: usize)
    where G: Graph
{

}