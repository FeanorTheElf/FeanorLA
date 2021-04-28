use std::collections::HashMap;
use std::collections::hash_map;
use std::marker::PhantomData;
use std::hash::Hash;
use std::iter::Peekable;

#[derive(Debug)]
pub enum SymbolMask<S> {
    Wildcard, Symbol(S)
}

#[derive(Debug)]
pub struct PrefixDatabaseNode<S> {
    is_word: bool,
    children: HashMap<S, Box<PrefixDatabaseNode<S>>>
}

pub struct PrefixDatabaseIter<'a, S, F, T> {
    mask: &'a [SymbolMask<S>],
    offset: usize,
    word: Vec<S>,
    stack: Vec<&'a PrefixDatabaseNode<S>>,
    iters: Vec<Peekable<hash_map::Iter<'a, S, Box<PrefixDatabaseNode<S>>>>>,
    converter: F,
    result: PhantomData<T>
}

impl<'a, S, F, T> PrefixDatabaseIter<'a, S, F, T> 
    where F: FnMut(&[S]) -> T, S: Eq + Hash + Clone + std::fmt::Debug
{
    fn should_continue(&mut self) -> bool {
        match &self.mask[self.offset] {
            SymbolMask::Symbol(_) => false,
            SymbolMask::Wildcard => self.iters.last_mut().unwrap().peek().is_some()
        }
    }

    fn do_continue(&mut self) {
        let (symbol, child) = self.iters.last_mut().unwrap().next().unwrap();
        self.word.push(symbol.clone());
        self.stack.push(child);
        self.offset += 1;
    }

    fn step_up(&mut self) -> Result<(), ()> {
        if self.offset == 0 {
            return Err(());
        }
        if let Some(SymbolMask::Wildcard) = self.mask.get(self.offset) {
            self.iters.pop();
        }
        self.offset -= 1;
        self.word.pop();
        self.stack.pop();
        return Ok(());
    }

    fn step(&mut self) -> Result<(), ()> {
        match self.mask.get(self.offset) {
            Some(SymbolMask::Symbol(s)) => {
                if let Some(child) = self.stack.last().unwrap().children.get(s) {
                    self.offset += 1;
                    self.word.push(s.clone());
                    self.stack.push(child);
                    return Ok(());
                }
            },
            Some(SymbolMask::Wildcard) => {
                let children_iter = self.stack.last().unwrap().children.iter().peekable();
                self.iters.push(children_iter);
                if self.iters.last_mut().unwrap().peek().is_some() {
                    self.do_continue();
                    return Ok(());
                }
            },
            None => {}
        };
        self.step_up()?;
        while !self.should_continue() {
            self.step_up()?;
        }
        self.do_continue();
        return Ok(());
    }
}

impl<'a, S, F, T> Iterator for PrefixDatabaseIter<'a, S, F, T> 
    where F: FnMut(&[S]) -> T, S: Eq + Hash + Clone + std::fmt::Debug
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.step().ok()?;
        while !self.stack.last().unwrap().is_word {
            self.step().ok()?;
        }
        return Some((self.converter)(&self.word[..]));
    }
}

#[derive(Debug)]
pub struct PrefixDatabase<S> {
    root: PrefixDatabaseNode<S>
}

impl<S> PrefixDatabase<S> 
    where S: Hash + Eq + Clone
{

    pub fn new() -> Self {
        PrefixDatabase {
            root: PrefixDatabaseNode {
                children: HashMap::new(),
                is_word: false
            }
        }
    }

    pub fn insert(&mut self, word: &[S]) {
        let mut current = &mut self.root;
        for s in word.iter() {
            current = current.children.entry(s.clone()).or_insert(Box::new(PrefixDatabaseNode {
                children: HashMap::new(),
                is_word: false
            }));
        }
        current.is_word = true;
    }

    ///
    /// Returns an iterator over all stored words that match the beginning of mask.
    /// 
    /// # Details
    /// 
    /// The iterator will contain an n-character word w, if and only if mask has at least
    /// n characters and the first n characters of mask either match exactly the characters
    /// of word, or are a wildcard. The iteration order is undefined.
    /// 
    pub fn query<'a, F, T>(&'a self, mask: &'a [SymbolMask<S>], converter: F) -> PrefixDatabaseIter<'a, S, F, T>
        where F: FnMut(&[S]) -> T
    {
        PrefixDatabaseIter {
            converter: converter,
            iters: Vec::new(),
            mask: mask,
            offset: 0,
            result: PhantomData,
            stack: vec![ &self.root ],
            word: Vec::new()
        }
    }
}

#[cfg(test)]
fn create_mask(s: &str) -> Vec<SymbolMask<char>> {
    s.chars().map(|c| if c == '_' { 
        SymbolMask::Wildcard 
    } else { 
        SymbolMask::Symbol(c) 
    }).collect()
}

#[cfg(test)]
use std::collections::HashSet;

#[test]
fn test_prefix_database() {
    let mut db: PrefixDatabase<char> = PrefixDatabase::new();
    db.insert(&"live".chars().collect::<Vec<_>>());
    db.insert(&"life".chars().collect::<Vec<_>>());
    db.insert(&"light".chars().collect::<Vec<_>>());
    db.insert(&"apple".chars().collect::<Vec<_>>());
    dbg!(&db);

    let converter = |x: &[char]| x.iter().collect::<String>();

    assert_eq!(vec![
        "live".to_owned(),
        "life".to_owned(),
        "light".to_owned()
    ].into_iter().collect::<HashSet<_>>(), db.query(&create_mask("li__t"), &converter).collect::<HashSet<String>>());

    assert_eq!(vec![
        "live".to_owned(),
        "life".to_owned(),
        "light".to_owned()
    ].into_iter().collect::<HashSet<_>>(), db.query(&create_mask("li__t..."), &converter).collect::<HashSet<String>>());

    assert_eq!(vec![
        "live".to_owned(),
        "life".to_owned()
    ].into_iter().collect::<HashSet<_>>(), db.query(&create_mask("li_e"), &converter).collect::<HashSet<String>>());

    assert_eq!(vec![].into_iter().collect::<HashSet<_>>(), db.query(&create_mask("a"), &converter).collect::<HashSet<String>>());
    assert_eq!(vec![
        "apple".to_owned()
    ].into_iter().collect::<HashSet<_>>(), db.query(&create_mask("app_e..."), &converter).collect::<HashSet<String>>());
}