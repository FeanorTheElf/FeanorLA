use std::collections::HashMap;

pub trait LifetimedGraph<'a> {
    type Node: 'a;
    type IncomingEdgeIter: Iterator<Item = Self::Node>;
    type OutgoingEdgeIter: Iterator<Item = Self::Node>;
    type NodeIter: Iterator<Item = Self::Node>;

    fn outgoing_edges(&'a self, node: Self::Node) -> Self::IncomingEdgeIter;
    fn incoming_edges(&'a self, node: Self::Node) -> Self::OutgoingEdgeIter;
    fn nodes(&'a self) -> Self::NodeIter;
}

pub trait Graph: for<'a> LifetimedGraph<'a> {}

impl<T> Graph for T
    where T: for<'a> LifetimedGraph<'a>
{}

pub struct MemoryGraph<T> {
    nodes: HashMap<usize, (T, Vec<usize>, Vec<usize>)>,
    counter: usize
}

impl<T> MemoryGraph<T> {

    pub fn new() -> MemoryGraph<T> {
        MemoryGraph {
            nodes: HashMap::new(),
            counter: 0
        }
    }

    pub fn get(&self, node: usize) -> &T {
        &self.nodes.get(&node).unwrap().0
    }

    pub fn get_mut(&mut self, node: usize) -> &mut T {
        &mut self.nodes.get_mut(&node).unwrap().0
    }

    pub fn add_node(&mut self, content: T) -> usize {
        let result = self.counter;
        self.counter += 1;
        self.nodes.insert(result, (content, Vec::new(), Vec::new()));
        return result;
    }

    pub fn delete_node(&mut self, node: usize) -> Option<T> {
        if let Some((res, _, incoming)) = self.nodes.remove(&node) {
            for u in &incoming {
                let edges = &mut self.nodes.get_mut(&u).unwrap().2;
                edges.remove(edges.iter().enumerate().find(|(_i, x)| **x == node).unwrap().0);
            }
            return Some(res);
        } else{
            return None;
        }
    }

    pub fn add_edge(&mut self, from: usize, to: usize) {
        self.remove_edge(from, to);
        self.nodes.get_mut(&from).unwrap().1.push(to);
    }

    pub fn remove_edge(&mut self, from: usize, to: usize) {
        let outgoing = &mut self.nodes.get_mut(&from).unwrap().1;
        if let Some((index, _)) = outgoing.iter().enumerate().find(|(_i, x)| **x == to) {
            outgoing.remove(index);
            
            let incoming = &mut self.nodes.get_mut(&to).unwrap().2;
            incoming.remove(incoming.iter().enumerate().find(|(_i, x)| **x == from).unwrap().0);
        }
    }
}

impl<'a, T> LifetimedGraph<'a> for MemoryGraph<T> {
    type Node = usize;
    type NodeIter = std::ops::Range<usize>;
    type OutgoingEdgeIter = std::iter::Copied<std::slice::Iter<'a, usize>>;
    type IncomingEdgeIter = std::iter::Copied<std::slice::Iter<'a, usize>>;

    fn outgoing_edges(&'a self, node: Self::Node) -> Self::OutgoingEdgeIter {
        self.nodes.get(&node).unwrap().1.iter().copied()
    }

    fn incoming_edges(&'a self, node: Self::Node) -> Self::IncomingEdgeIter {
        self.nodes.get(&node).unwrap().1.iter().copied()
    }

    fn nodes(&'a self) -> Self::NodeIter {
        0..self.nodes.len()
    }
}