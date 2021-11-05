use std::collections::HashMap;
use std::collections::hash_map;
use std::marker::PhantomData;
use std::hash::Hash;
use std::iter::Peekable;
use std::cmp::max;
use std::fs::File;
use std::io::prelude::*;

fn read_words_file() -> Result<String, std::io::Error> {
    let mut file = File::open("src/combinatorics/words.txt")?;
    let mut content = "".to_owned();
    file.read_to_string(&mut content)?;
    return Ok(content);
}

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
    where F: FnMut(&[S]) -> T, S: Eq + Hash + Clone
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
    where F: FnMut(&[S]) -> T, S: Eq + Hash + Clone
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
    where S: Eq + Hash + Clone
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

fn evaluate<S>(text: &[SymbolMask<S>], words: &PrefixDatabase<S>) -> u64 
    where S: Eq + Hash + Clone + std::fmt::Debug
{
    let mut partial_sols = Vec::new();
    partial_sols.resize(text.len() + 1, 0u64);
    for i in 0..text.len() {
        words.query(&text[i..], |word| {
            partial_sols[i + word.len()] = max(
                partial_sols[i + word.len()],
                partial_sols[i] + word.len() as u64
            );
            return ();
        }).for_each(|_| {});
        partial_sols[i + 1] = max(
            partial_sols[i + 1],
            partial_sols[i]
        );
    }
    return partial_sols[text.len()];
}

fn get_word_db() -> PrefixDatabase<char> {
    let mut result = PrefixDatabase::new();
    let words_file = read_words_file().unwrap();
    for w in words_file.split("\r\n") {
        result.insert(&w.chars().map(|c| c.to_uppercase().next().unwrap()).collect::<Vec<_>>()[..]);
    }
    return result;
}

fn decrypt(ciphertext: &str) -> String {
    let english_words = get_word_db();

    fn letter_to_number(c: char) -> u8 {
        ((c.to_uppercase().next().unwrap() as u32) - ('A' as u32)) as u8
    }

    fn number_to_letter(c: u8) -> char {
        char::from_u32(c as u32 + 'A' as u32).unwrap()
    }

    fn apply<F: FnMut(char) -> char>(string: &str, f: F) -> Vec<SymbolMask<char>> {
        string.chars().map(f).map(SymbolMask::Symbol).collect()
    }

    let mut best_solution = (0, Vec::new());
    for shift in 0..26 {
        let shifted = apply(ciphertext, |c| number_to_letter((letter_to_number(c) + shift) % 26));
        let evaluation = evaluate(&shifted[..], &english_words);
        if evaluation > best_solution.0 {
            best_solution = (evaluation, shifted);
        }
    }
    return best_solution.1.into_iter().map(|m| match m { SymbolMask::Symbol(c) => c, SymbolMask::Wildcard => '_' }).collect();
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

#[test]
fn test_evaluate() {
    let mut db: PrefixDatabase<char> = PrefixDatabase::new();
    db.insert(&"live".chars().collect::<Vec<_>>());
    db.insert(&"life".chars().collect::<Vec<_>>());
    db.insert(&"light".chars().collect::<Vec<_>>());
    db.insert(&"apple".chars().collect::<Vec<_>>());

    let mask = create_mask("____v__l__e__");
    assert_eq!(8, evaluate(&mask, &db));
}
