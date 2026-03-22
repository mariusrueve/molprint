use crate::mol::atom::Element;

/// Primitive atom predicates for SMARTS matching.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomPrimitive {
    /// Match any atom: `*`.
    Any,
    /// Match by atomic number: `#n`.
    AtomicNum(u8),
    /// Match a specific element as aliphatic (uppercase symbol, e.g. `C`, `N`).
    ElementAliphatic(Element),
    /// Match a specific element as aromatic (lowercase symbol, e.g. `c`, `n`).
    ElementAromatic(Element),
    /// Match any aromatic atom: `a`.
    Aromatic,
    /// Match any aliphatic atom: `A`.
    Aliphatic,
    /// Match total H count: `H` (= 1) or `Hn`.
    HCount(u8),
    /// Match formal charge: `+n` or `-n`.
    Charge(i8),
    /// Atom is in any ring: `R`.  `Some(n)` = in ring of size n.
    Ring(Option<usize>),
    /// Atom is not in any ring: `R0`.
    NotInRing,
    /// Heavy-atom degree: `Dn`.
    Degree(u8),
}

/// Atom expression tree supporting logical composition.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomExpr {
    /// A single primitive predicate.
    Primitive(AtomPrimitive),
    /// Logical NOT.
    Not(Box<AtomExpr>),
    /// Logical AND (both explicit `&` and implicit adjacency inside `[…]`).
    And(Box<AtomExpr>, Box<AtomExpr>),
    /// Logical OR (`,`).
    Or(Box<AtomExpr>, Box<AtomExpr>),
}

/// SMARTS bond predicate.
#[derive(Debug, Clone, PartialEq)]
pub enum SmartsBond {
    /// `-` single bond.
    Single,
    /// `=` double bond.
    Double,
    /// `#` triple bond.
    Triple,
    /// `:` aromatic bond.
    Aromatic,
    /// `~` any bond.
    Any,
    /// No explicit bond token (matches single or aromatic, context-dependent).
    Unspecified,
}

/// A single atom slot in a SMARTS pattern.
#[derive(Debug, Clone)]
pub struct SmartsAtom {
    /// Predicate that a molecule atom must satisfy to match this slot.
    pub expr: AtomExpr,
}

/// A compiled SMARTS pattern.
#[derive(Debug, Clone, Default)]
pub struct SmartsPattern {
    /// Atoms in the order they appear in the SMARTS string.
    pub atoms: Vec<SmartsAtom>,
    /// Edges: `(from_idx, to_idx, bond_predicate)`.
    pub edges: Vec<(usize, usize, SmartsBond)>,
}

impl SmartsPattern {
    /// Create an empty pattern.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of atom slots.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Neighbours of pattern atom `idx`: iterator of `(neighbour_idx, bond_predicate)`.
    pub fn neighbors(&self, idx: usize) -> impl Iterator<Item = (usize, &SmartsBond)> {
        self.edges.iter().filter_map(move |(a, b, bond)| {
            if *a == idx {
                Some((*b, bond))
            } else if *b == idx {
                Some((*a, bond))
            } else {
                None
            }
        })
    }
}
