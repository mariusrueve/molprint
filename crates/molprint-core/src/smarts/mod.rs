//! SMARTS substructure query language — minimal v1.
//!
//! # Supported primitives
//!
//! **Atom** (inside `[…]`):
//! - `#n`  — atomic number
//! - `C`, `N`, `O`, … — aliphatic element (organic subset)
//! - `c`, `n`, `o`, … — aromatic element
//! - `a` / `A`         — any aromatic / any aliphatic
//! - `Hn`              — H count (`H` = 1, `H2` = 2, …)
//! - `+n` / `-n`       — formal charge
//! - `R`               — in any ring; `R0` = not in ring; `Rn` = ring of size n
//! - `Dn`              — heavy-atom degree
//! - `*`               — any atom
//! - `!`, `&`, `;`, `,` — logical NOT, AND (high/low precedence), OR
//!
//! **Bond** (between atoms):
//! - `-` `=` `#` `:` — single, double, triple, aromatic
//! - `~`             — any bond
//! - (none)          — matches single or aromatic
//!
//! **Not supported in v1**: recursive SMARTS `$()`, stereo `@`/`@@`, component grouping.
//!
//! # Example
//!
//! ```
//! use molprint_core::smarts;
//! use molprint_core::smiles::parse_smiles;
//!
//! let pattern = smarts::compile("[#6]~[#7]").unwrap();
//! let mol = parse_smiles("c1ccncc1").unwrap(); // pyridine
//! assert!(smarts::has_match(&mol, &pattern));
//! ```

mod ast;
mod lexer;
mod matcher;

pub use ast::{AtomExpr, AtomPrimitive, SmartsAtom, SmartsBond, SmartsPattern};
pub use matcher::{has_match, subgraph_match};

use thiserror::Error;

/// Errors that can occur while compiling a SMARTS string.
#[derive(Error, Debug, PartialEq)]
pub enum SmartsError {
    #[error("unexpected character '{0}' at position {1}")]
    UnexpectedChar(char, usize),

    #[error("unterminated bracket '[' at position {0}")]
    UnterminatedBracket(usize),

    #[error("unmatched ring closure {0}")]
    UnmatchedRing(u8),

    #[error("unmatched parenthesis")]
    UnmatchedParen,

    #[error("empty SMARTS string")]
    Empty,

    #[error("invalid atomic number: '{0}'")]
    InvalidAtomicNum(String),
}

/// Compile a SMARTS string into a matchable [`SmartsPattern`].
///
/// # Errors
///
/// Returns [`SmartsError`] if the SMARTS string is syntactically invalid.
pub fn compile(smarts: &str) -> Result<SmartsPattern, SmartsError> {
    lexer::parse_smarts(smarts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::graph::MolGraph;
    use crate::smiles::parse_smiles;

    fn mol(smi: &str) -> MolGraph {
        parse_smiles(smi).unwrap()
    }

    // ── compile / syntax ──────────────────────────────────────────────────

    #[test]
    fn test_compile_empty() {
        assert!(matches!(compile(""), Err(SmartsError::Empty)));
    }

    #[test]
    fn test_compile_wildcard() {
        let p = compile("*").unwrap();
        assert_eq!(p.num_atoms(), 1);
    }

    #[test]
    fn test_compile_atomic_num() {
        let p = compile("[#6]").unwrap();
        assert_eq!(p.num_atoms(), 1);
    }

    #[test]
    fn test_compile_chain() {
        let p = compile("[#6]~[#7]").unwrap();
        assert_eq!(p.num_atoms(), 2);
        assert_eq!(p.edges.len(), 1);
    }

    #[test]
    fn test_compile_ring() {
        let p = compile("c1ccccc1").unwrap();
        assert_eq!(p.num_atoms(), 6);
        assert_eq!(p.edges.len(), 6); // 5 chain + 1 ring-closure
    }

    #[test]
    fn test_compile_branch() {
        let p = compile("C(C)C").unwrap();
        assert_eq!(p.num_atoms(), 3);
        assert_eq!(p.edges.len(), 2);
    }

    #[test]
    fn test_compile_not_in_bracket() {
        let p = compile("[!#6]").unwrap();
        assert_eq!(p.num_atoms(), 1);
    }

    #[test]
    fn test_compile_implicit_and() {
        // [NH] = nitrogen AND H-count=1
        let p = compile("[NH]").unwrap();
        assert_eq!(p.num_atoms(), 1);
        assert!(matches!(p.atoms[0].expr, AtomExpr::And(_, _)));
    }

    #[test]
    fn test_compile_or() {
        let p = compile("[N,O]").unwrap();
        assert_eq!(p.num_atoms(), 1);
        assert!(matches!(p.atoms[0].expr, AtomExpr::Or(_, _)));
    }

    #[test]
    fn test_unmatched_ring_error() {
        assert!(matches!(
            compile("C1CC"),
            Err(SmartsError::UnmatchedRing(_))
        ));
    }

    #[test]
    fn test_unmatched_paren_error() {
        assert!(matches!(compile("C(C"), Err(SmartsError::UnmatchedParen)));
    }

    // ── has_match ─────────────────────────────────────────────────────────

    #[test]
    fn test_match_any_atom() {
        let p = compile("*").unwrap();
        assert!(has_match(&mol("C"), &p));
        assert!(has_match(&mol("O"), &p));
    }

    #[test]
    fn test_match_atomic_number() {
        let p = compile("[#6]").unwrap();
        assert!(has_match(&mol("C"), &p));
        assert!(!has_match(&mol("O"), &p));
    }

    #[test]
    fn test_match_aromatic_carbon() {
        let benzene = mol("c1ccccc1");
        let ethane = mol("CC");
        assert!(has_match(&benzene, &compile("c").unwrap()));
        assert!(!has_match(&ethane, &compile("c").unwrap()));
        assert!(has_match(&ethane, &compile("C").unwrap()));
    }

    #[test]
    fn test_match_any_aromatic() {
        let p = compile("[a]").unwrap();
        assert!(has_match(&mol("c1ccccc1"), &p));
        assert!(!has_match(&mol("CCCC"), &p));
    }

    #[test]
    fn test_match_ring_membership() {
        let ring_p = compile("[R]").unwrap();
        let no_ring_p = compile("[R0]").unwrap();
        assert!(has_match(&mol("C1CC1"), &ring_p));
        assert!(!has_match(&mol("CCC"), &ring_p));
        assert!(has_match(&mol("CCC"), &no_ring_p));
        assert!(!has_match(&mol("C1CCC1"), &no_ring_p));
    }

    #[test]
    fn test_match_charge() {
        let pos_p = compile("[+]").unwrap();
        assert!(has_match(&mol("[NH4+]"), &pos_p));
        assert!(!has_match(&mol("N"), &pos_p));
    }

    #[test]
    fn test_match_h_count() {
        let p = compile("[H2]").unwrap();
        // Ethanol C1: has 3H; C2: has 2H — should match
        assert!(has_match(&mol("CCO"), &p));
    }

    #[test]
    fn test_match_or_expression() {
        let p = compile("[N,O]").unwrap();
        assert!(has_match(&mol("CCO"), &p));
        assert!(has_match(&mol("CCN"), &p));
        assert!(!has_match(&mol("CCC"), &p));
    }

    #[test]
    fn test_match_not_expression() {
        let p = compile("[!#6]").unwrap();
        // ethanol has O, not just C
        assert!(has_match(&mol("CCO"), &p));
        // ethane has only C
        assert!(!has_match(&mol("CC"), &p));
    }

    #[test]
    fn test_match_bond_types() {
        let double = compile("C=O").unwrap();
        assert!(has_match(&mol("CC(=O)O"), &double));
        assert!(!has_match(&mol("CCO"), &double));

        let triple = compile("C#N").unwrap();
        assert!(has_match(&mol("CC#N"), &triple));
        assert!(!has_match(&mol("CCN"), &triple));
    }

    #[test]
    fn test_match_any_bond() {
        let p = compile("[#6]~[#7]").unwrap();
        assert!(has_match(&mol("c1ccncc1"), &p)); // pyridine: aromatic C-N
        assert!(has_match(&mol("CCN"), &p)); // single C-N
    }

    #[test]
    fn test_match_benzene_ring_pattern() {
        let p = compile("c1ccccc1").unwrap();
        assert!(has_match(&mol("c1ccccc1"), &p));
        assert!(!has_match(&mol("CCCCCC"), &p));
    }

    #[test]
    fn test_match_substructure_in_larger_mol() {
        // carbonyl in aspirin
        let p = compile("C=O").unwrap();
        assert!(has_match(&mol("CC(=O)Oc1ccccc1C(=O)O"), &p));
    }

    #[test]
    fn test_match_degree() {
        // Nitrogen with 3 heavy-atom neighbours
        let p = compile("[ND3]").unwrap();
        assert!(has_match(&mol("CN(C)C"), &p));
        assert!(!has_match(&mol("CN"), &p));
    }

    // ── subgraph_match ────────────────────────────────────────────────────

    #[test]
    fn test_subgraph_match_count() {
        let p = compile("[#6]").unwrap();
        let matches = subgraph_match(&mol("CCO"), &p);
        // 2 carbons in ethanol
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_subgraph_match_empty_pattern() {
        let p = compile_empty();
        let matches = subgraph_match(&mol("CCO"), &p);
        assert_eq!(matches.len(), 1);
        assert!(matches[0].is_empty());
    }

    fn compile_empty() -> SmartsPattern {
        SmartsPattern::new()
    }
}
