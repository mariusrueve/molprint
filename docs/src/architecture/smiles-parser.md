# SMILES Parser

The SMILES parser converts a SMILES string into a `MolGraph`. It is split into a lexer (tokenizer) and a parser (graph builder), following the standard compiler front-end pattern.

## Lexer (`smiles/lexer.rs`)

The lexer reads the SMILES string character by character and produces a stream of `Token` values.

```rust
pub enum Token {
    Atom { element: Element, aromatic: bool },
    BracketAtom { isotope, element, aromatic, h_count, charge, map_num, chirality },
    Bond(BondType),
    BondStereo(char),     // '/' or '\'
    RingClosure(u8),      // ring closure digit or %nn
    OpenBranch,           // '('
    CloseBranch,          // ')'
    Dot,                  // '.'
}
```

Key tokenization rules:

- **Organic subset atoms** (`B`, `C`, `N`, `O`, `P`, `S`, `F`, `Cl`, `Br`, `I`) and their lowercase aromatic counterparts (`c`, `n`, `o`, `s`, `p`) are emitted as `Token::Atom`.
- **Bracket atoms** (`[NH3]`, `[2H]`, `[Fe+2]`, `[13C@@H]`) are parsed as `Token::BracketAtom` with all fields decoded inline. Two-character element symbols (`Cl`, `Br`, `Se`, etc.) require look-ahead.
- **Ring closures** are `%NN` for two-digit ring numbers (e.g., `%10`) or single digits (0–9).
- **Stereo tokens** (`/`, `\`) are consumed but stored only as `BondStereo`; stereo is not used by fingerprinting algorithms.

## Parser (`smiles/parser.rs`)

The parser takes the token stream and builds a `MolGraph` incrementally.

State maintained during parsing:

```rust
let mut current: Option<NodeIndex> = None;
let mut branch_stack: Vec<NodeIndex> = Vec::new();
let mut pending_bond: Option<BondType> = None;
let mut ring_map: HashMap<u8, (NodeIndex, Option<BondType>)> = HashMap::new();
```

The token dispatch:

| Token | Action |
|---|---|
| `Atom` / `BracketAtom` | Add node; connect to `current` with `pending_bond` or default bond; set `current` |
| `Bond` | Store in `pending_bond` |
| `BondStereo` | Ignored (no pending_bond change) |
| `OpenBranch` | Push `current` onto `branch_stack` |
| `CloseBranch` | Pop `branch_stack`; restore `current`; clear `pending_bond` |
| `RingClosure` | If ring_id seen before: add closure edge; else store in ring_map |
| `Dot` | Set `current = None` (start new fragment) |

### Default bond inference

When no explicit bond token precedes an atom, the parser infers the bond type:

- If both the previous atom and the new atom are aromatic → `BondType::Aromatic`
- Otherwise → `BondType::Single`

### Ring closure bond type

When a ring closure digit appears twice, the bond type is resolved in this priority order:

1. The `pending_bond` at the closing digit (if set explicitly)
2. The bond type stored when the ring was opened (if set explicitly)
3. Inferred from the aromaticity of the two ring-closure atoms

### Error handling

The parser returns `Result<MolGraph, ParseError>` where `ParseError` covers:
- `Lexer(LexerError)` — invalid characters in the SMILES string
- `UnmatchedRing(u8)` — a ring closure digit that was never closed
- `UnmatchedBranch` — mismatched `(` or `)`
- `Empty` — empty input string
- `DanglingBond` — a bond token with no subsequent atom

### Post-processing

After the token loop, the parser runs two passes:

1. `graph.assign_implicit_hydrogens()` — computes `h_count` for every atom
2. `assign_ring_info(&mut graph)` — runs SSSR and sets `in_ring` and `ring_sizes` on atoms

These are called unconditionally, so every `MolGraph` produced by `parse_smiles` has fully-populated ring and hydrogen information.

## Public API

```rust
// smiles/mod.rs
pub use parser::parse as parse_smiles;
```

Usage:

```rust
use molprint_core::smiles::parse_smiles;

let mol = parse_smiles("CC(=O)O").unwrap(); // acetic acid
```

## Fuzz testing

The SMILES lexer and parser are covered by a libFuzzer fuzz target at `fuzz/fuzz_targets/fuzz_smiles.rs`. The target feeds arbitrary byte strings into `parse_smiles` and checks that it never panics — only returns `Ok` or `Err`. This uncovered several edge cases during development, particularly around multi-digit ring closures and unusual bracket atom encodings.
