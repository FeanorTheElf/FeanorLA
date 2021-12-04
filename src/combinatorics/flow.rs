use super::graph::*;

use std::cmp::min;
use std::collections::{BinaryHeap, HashMap, HashSet};

type Node<'a, G> = <G as LifetimedGraph<'a>>::Node;

struct SizedElement<T>(u32, T);

impl<T> PartialEq for SizedElement<T> {
    fn eq(&self, rhs: &Self) -> bool {
        self.0 == rhs.0
    }
}

impl<T> Eq for SizedElement<T> {}

impl<T> PartialOrd for SizedElement<T> {
    fn partial_cmp(&self, rhs: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(rhs))
    }
}

impl<T> Ord for SizedElement<T> {
    fn cmp(&self, rhs: &Self) -> std::cmp::Ordering {
        self.0.cmp(&rhs.0)
    }
}

///
/// Finds an augmenting path for the given flow network or the nodes
/// that are reachable from source.
/// 
/// The iterator over the augmenting path does so backwards.  
/// 
fn find_augmenting_path<'a, G>(
    graph: &'a G, 
    s: &HashSet<Node<'a, G>>, 
    t: &HashSet<Node<'a, G>>, 
    capacity: &HashMap<Node<'a, G>, u32>
) -> Result<(u32, impl 'a + Iterator<Item = Node<'a, G>>), impl 'a + Iterator<Item = Node<'a, G>>> 
    where G: Graph, 
        <G as LifetimedGraph<'a>>::Node: Clone + std::hash::Hash + Eq
{
    let mut open: BinaryHeap<SizedElement<(Node<'a, G>, Option<Node<'a, G>>)>> = BinaryHeap::new();
    let mut closed: HashMap<Node<'a, G>, (u32, Option<Node<'a, G>>)> = HashMap::new();
    open.extend(s.iter()
        .map(|n| SizedElement(*capacity.get(n).unwrap(), (n.clone(), None)))
        .filter(|x| x.0 > 0)
    );

    while let Some(SizedElement(flow, (node, pred))) = open.pop() {
        if t.contains(&node) {
            let mut current = (Some(node), pred);
            return Ok((flow, std::iter::from_fn(move || {
                let mut tmp = current.1.as_ref().map(|x| closed.remove(x).unwrap()).and_then(|x| x.1);
                std::mem::swap(&mut current.0, &mut current.1);
                std::mem::swap(&mut current.1, &mut tmp);
                tmp
            })));
        }
        for n in graph.outgoing_edges(node.clone()) {
            let free_capacity = *capacity.get(&n).unwrap();
            let path_capacity = min(free_capacity, flow);

            if path_capacity > 0 {
                let succ_entry = (n, Some(node.clone()));
                // reopen the node if a better path was found
                if let Some((last_flow, _)) = closed.get(&succ_entry.0) {
                    if *last_flow < path_capacity {
                        closed.remove(&succ_entry.0);
                        open.push(SizedElement(path_capacity, succ_entry));
                    }
                } else {
                    open.push(SizedElement(path_capacity, succ_entry));
                }
            }
        }
        closed.insert(node, (flow, pred));
    }
    return Err(closed.into_iter().map(|x| x.0));
}

///
/// Calculates the maximum-flow resp min s-t-cut problem in the given
/// (directed, positive-integer-edge-weighted) graph.
/// 
/// !THIS IS CURRENTLY WRONG!
/// 
/// Returns the free edge capacities after having reached a maximum flow
/// and a compatible min s-t-cut
/// 
pub fn ford_fulkerson<'a, G>(
    graph: &'a G, 
    s: &HashSet<Node<'a, G>>, 
    t: &HashSet<Node<'a, G>>, 
    mut capacity: HashMap<Node<'a, G>, u32>
) -> (impl Iterator<Item = Node<'a, G>>, HashMap<Node<'a, G>, u32>)
    where G: Graph,
        <G as LifetimedGraph<'a>>::Node: Clone + std::hash::Hash + Eq
{
    assert!(s.len() > 0);
    assert!(t.len() > 0);
    loop {
        match find_augmenting_path(graph, s, t, &capacity) {
            Ok((flow, path)) => {
                for node in path {
                    let capacity_entry = capacity.get_mut(&node).unwrap();
                    assert!(*capacity_entry >= flow);
                    *capacity_entry -= flow;
                }
            },
            Err(reachable_from_source) => {
                return (reachable_from_source, capacity);
            }
        }
    }
}

#[test]
fn test_ford_fulkerson() {
    // Consider the graph
    // 4 --- 2 --- 4 --- 6
    //  \         /
    //   --- 1 --
    let mut graph = MemoryGraph::new();
    let s = graph.add_node(("s", 4));
    let a = graph.add_node(("a", 2));
    let b = graph.add_node(("b", 4));
    let c = graph.add_node(("c", 1));
    let t = graph.add_node(("t", 6));
    graph.add_edge(s, a);
    graph.add_edge(a, b);
    graph.add_edge(b, t);
    graph.add_edge(s, c);
    graph.add_edge(c, b);
    let mut capacity = HashMap::new();
    for n in &[s, t, a, b, c] {
        capacity.insert(*n, graph.get(*n).1);
    }

    let (cut, max_flow) = ford_fulkerson(
        &graph, 
        &[s].iter().copied().collect(), 
        &[t].iter().copied().collect(), 
        capacity
    );
    assert_eq!(vec![s], cut.collect::<Vec<_>>());
    assert_eq!(1, *max_flow.get(&s).unwrap());
    assert_eq!(0, *max_flow.get(&a).unwrap());
    assert_eq!(0, *max_flow.get(&c).unwrap());
    assert_eq!(1, *max_flow.get(&b).unwrap());
    assert_eq!(3, *max_flow.get(&t).unwrap());
}
