# SMARTS Engine

The SMARTS engine supports substructure queries against `MolGraph`. It is used internally by the MACCS-166 fingerprint for structural key evaluation and is also available as a public API.

## Supported features

**Atom primitives** (inside `[…]`):

| Primitive | Meaning |
|---|---|
| `#n` | Atomic number (e.g., `[#6]` = carbon) |
| `C`, `N`, `O`, … | Aliphatic element (organic subset) |
| `c`, `n`, `o`, … | Aromatic element |
| `a` / `A` | Any aromatic / any aliphatic atom |
| `Hn` | Hydrogen count (`H` = 1, `H2` = 2) |
| `+n` / `-n` | Formal charge |
| `R` | In any ring; `R0` = not in ring |
| `Dn` | Heavy-atom degree |
| `*` | Any atom |

**Logical operators**: `!` (NOT), `&` (high-priority AND), `;` (low-priority AND), `,` (OR).

**Bond primitives**: `-` `=` `#` `:` `~` (any), and unspecified (matches single or aromatic).

**Not supported in v1**: recursive SMARTS `$()`, stereo `@`/`@@`, component-level grouping.

## Architecture

The engine has three layers:

1. **AST** (`smarts/ast.rs`) — data structures for atom and bond expressions
2. **Lexer/Parser** (`smarts/lexer.rs`) — parses a SMARTS string into a `SmartsPattern`
3. **Matcher** (`smarts/matcher.rs`) — VF2-style subgraph isomorphism

### AST

```rust
pub enum AtomPrimitive {
    AtomicNum(u8),
    Aromatic(bool),
    AnyAromatic,
    AnyAliphatic,
    HCount(u8),
    Charge(i8),
    InRing,
    NotInRing,
    RingSize(u8),
    Degree(u8),
    Any,
}

pub enum AtomExpr {
    Primitive(AtomPrimitive),
    Not(Box<AtomExpr>),
    And(Box<AtomExpr>, Box<AtomExpr>),
    Or(Box<AtomExpr>, Box<AtomExpr>),
}
```

`SmartsBond` stores the bond constraint (specific type, any, or default). `SmartsPattern` is a graph of `SmartsAtom` nodes and `SmartsBond` edges.

### Matching

The matcher uses a recursive VF2-style algorithm. It maintains a partial mapping from pattern atoms to molecule atoms and extends it one atom at a time.

```rust
pub fn has_match(mol: &MolGraph, pattern: &SmartsPattern) -> bool;
pub fn subgraph_match(mol: &MolGraph, pattern: &SmartsPattern) -> Vec<Vec<(usize, usize)>>;
```

`has_match` returns early on the first valid mapping. `subgraph_match` collects all mappings.

The atom compatibility check evaluates the `AtomExpr` tree recursively:

```rust
fn atom_matches(mol_atom: &Atom, expr: &AtomExpr) -> bool {
    match expr {
        AtomExpr::Primitive(p) => primitive_matches(mol_atom, p),
        AtomExpr::Not(inner) => !atom_matches(mol_atom, inner),
        AtomExpr::And(a, b) => atom_matches(mol_atom, a) && atom_matches(mol_atom, b),
        AtomExpr::Or(a, b) => atom_matches(mol_atom, a) || atom_matches(mol_atom, b),
    }
}
```

## Usage

```rust
use molprint_core::smarts;
use molprint_core::smiles::parse_smiles;

// Carbon bonded to nitrogen (any bond type)
let pattern = smarts::compile("[#6]~[#7]").unwrap();
let pyridine = parse_smiles("c1ccncc1").unwrap();
assert!(smarts::has_match(&pyridine, &pattern));

// Hydroxyl group
let oh = smarts::compile("[OH]").unwrap();
let ethanol = parse_smiles("CCO").unwrap();
assert!(smarts::has_match(&ethanol, &oh));

// Count amine nitrogens
let amine = smarts::compile("[NH2]").unwrap();
let diamine = parse_smiles("NCCN").unwrap();
let matches = smarts::subgraph_match(&diamine, &amine);
assert_eq!(matches.len(), 2);
```

## Compile errors

`smarts::compile` returns `Result<SmartsPattern, SmartsError>`:

```rust
pub enum SmartsError {
    UnexpectedChar(char, usize),
    UnterminatedBracket(usize),
    UnmatchedRing(u8),
    UnmatchedParen,
    Empty,
    InvalidAtomicNum(String),
}
```
