//! SMARTS parser: converts a SMARTS string into a [`SmartsPattern`].
//!
//! Top-level syntax mirrors SMILES (atoms, bonds, branches, ring closures).
//! Atom expressions inside `[…]` support primitives, `!`, `&`, `;`, `,`.

use super::ast::{AtomExpr, AtomPrimitive, SmartsAtom, SmartsBond, SmartsPattern};
use super::SmartsError;
use crate::mol::atom::Element;
use std::collections::HashMap;

/// Parse a SMARTS string into a [`SmartsPattern`].
pub(super) fn parse_smarts(smarts: &str) -> Result<SmartsPattern, SmartsError> {
    if smarts.is_empty() {
        return Err(SmartsError::Empty);
    }

    let chars: Vec<char> = smarts.chars().collect();
    let mut pos = 0usize;
    let mut pattern = SmartsPattern::new();
    let mut current: Option<usize> = None;
    let mut branch_stack: Vec<usize> = Vec::new();
    let mut pending_bond: Option<SmartsBond> = None;
    let mut ring_map: HashMap<u8, (usize, Option<SmartsBond>)> = HashMap::new();

    while pos < chars.len() {
        match chars[pos] {
            // ── Bracket atom ───────────────────────────────────────────
            '[' => {
                let bracket_start = pos;
                pos += 1; // consume '['
                let (expr, new_pos) = parse_atom_expr(&chars, pos)?;
                if new_pos >= chars.len() || chars[new_pos] != ']' {
                    return Err(SmartsError::UnterminatedBracket(bracket_start));
                }
                pos = new_pos + 1; // consume ']'
                let idx = push_atom(
                    &mut pattern,
                    &mut current,
                    &mut pending_bond,
                    SmartsAtom { expr },
                );
                current = Some(idx);
            }

            // ── Branches ───────────────────────────────────────────────
            '(' => {
                if let Some(cur) = current {
                    branch_stack.push(cur);
                }
                pos += 1;
            }
            ')' => {
                match branch_stack.pop() {
                    Some(prev) => current = Some(prev),
                    None => return Err(SmartsError::UnmatchedParen),
                }
                pending_bond = None;
                pos += 1;
            }

            // ── Disconnected component ─────────────────────────────────
            '.' => {
                current = None;
                pending_bond = None;
                pos += 1;
            }

            // ── Bond tokens ────────────────────────────────────────────
            '-' => {
                pending_bond = Some(SmartsBond::Single);
                pos += 1;
            }
            '=' => {
                pending_bond = Some(SmartsBond::Double);
                pos += 1;
            }
            ':' => {
                pending_bond = Some(SmartsBond::Aromatic);
                pos += 1;
            }
            '~' => {
                pending_bond = Some(SmartsBond::Any);
                pos += 1;
            }

            // ── Wildcard atom ──────────────────────────────────────────
            '*' => {
                let atom = SmartsAtom {
                    expr: AtomExpr::Primitive(AtomPrimitive::Any),
                };
                let idx = push_atom(&mut pattern, &mut current, &mut pending_bond, atom);
                current = Some(idx);
                pos += 1;
            }

            // ── Aliphatic organic-subset atoms ─────────────────────────
            'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' => {
                let (elem, new_pos) = parse_aliphatic_atom(&chars, pos);
                let atom = SmartsAtom {
                    expr: AtomExpr::Primitive(AtomPrimitive::ElementAliphatic(elem)),
                };
                let idx = push_atom(&mut pattern, &mut current, &mut pending_bond, atom);
                current = Some(idx);
                pos = new_pos;
            }

            // ── Aromatic organic-subset atoms ──────────────────────────
            'b' | 'c' | 'n' | 'o' | 'p' | 's' => {
                let elem = match chars[pos] {
                    'b' => Element::B,
                    'c' => Element::C,
                    'n' => Element::N,
                    'o' => Element::O,
                    'p' => Element::P,
                    's' => Element::S,
                    _ => unreachable!(),
                };
                let atom = SmartsAtom {
                    expr: AtomExpr::Primitive(AtomPrimitive::ElementAromatic(elem)),
                };
                let idx = push_atom(&mut pattern, &mut current, &mut pending_bond, atom);
                current = Some(idx);
                pos += 1;
            }

            // ── Ring-closure digits ────────────────────────────────────
            '0'..='9' => {
                let ring_id = (chars[pos] as u8) - b'0';
                handle_ring(
                    &mut pattern,
                    &mut ring_map,
                    &mut current,
                    &mut pending_bond,
                    ring_id,
                    pos,
                )?;
                pos += 1;
            }
            '%' => {
                pos += 1;
                if pos + 1 >= chars.len()
                    || !chars[pos].is_ascii_digit()
                    || !chars[pos + 1].is_ascii_digit()
                {
                    return Err(SmartsError::UnexpectedChar('%', pos - 1));
                }
                let ring_id = (chars[pos] as u8 - b'0') * 10 + (chars[pos + 1] as u8 - b'0');
                pos += 2;
                handle_ring(
                    &mut pattern,
                    &mut ring_map,
                    &mut current,
                    &mut pending_bond,
                    ring_id,
                    pos,
                )?;
            }

            // ── Triple bond (top-level '#') ─────────────────────────────
            '#' => {
                pending_bond = Some(SmartsBond::Triple);
                pos += 1;
            }

            ch => return Err(SmartsError::UnexpectedChar(ch, pos)),
        }
    }

    if let Some((&ring_id, _)) = ring_map.iter().next() {
        return Err(SmartsError::UnmatchedRing(ring_id));
    }
    if !branch_stack.is_empty() {
        return Err(SmartsError::UnmatchedParen);
    }

    Ok(pattern)
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Add a new atom to the pattern, wiring it to `current` with `pending_bond`.
/// Returns the new atom's index.
fn push_atom(
    pattern: &mut SmartsPattern,
    current: &mut Option<usize>,
    pending_bond: &mut Option<SmartsBond>,
    atom: SmartsAtom,
) -> usize {
    let idx = pattern.atoms.len();
    pattern.atoms.push(atom);
    if let Some(prev) = *current {
        let bond = pending_bond.take().unwrap_or(SmartsBond::Unspecified);
        pattern.edges.push((prev, idx, bond));
    }
    idx
}

fn handle_ring(
    pattern: &mut SmartsPattern,
    ring_map: &mut HashMap<u8, (usize, Option<SmartsBond>)>,
    current: &mut Option<usize>,
    pending_bond: &mut Option<SmartsBond>,
    ring_id: u8,
    pos: usize,
) -> Result<(), SmartsError> {
    let cur = current.ok_or(SmartsError::UnexpectedChar('0', pos))?;
    if let Some((other, ring_bond)) = ring_map.remove(&ring_id) {
        let bond = pending_bond
            .take()
            .or(ring_bond)
            .unwrap_or(SmartsBond::Unspecified);
        pattern.edges.push((other, cur, bond));
    } else {
        ring_map.insert(ring_id, (cur, pending_bond.take()));
    }
    Ok(())
}

/// Parse an aliphatic element from the organic subset, handling `Br` and `Cl`.
fn parse_aliphatic_atom(chars: &[char], pos: usize) -> (Element, usize) {
    match chars[pos] {
        'B' if pos + 1 < chars.len() && chars[pos + 1] == 'r' => (Element::Br, pos + 2),
        'B' => (Element::B, pos + 1),
        'C' if pos + 1 < chars.len() && chars[pos + 1] == 'l' => (Element::Cl, pos + 2),
        'C' => (Element::C, pos + 1),
        'N' => (Element::N, pos + 1),
        'O' => (Element::O, pos + 1),
        'P' => (Element::P, pos + 1),
        'S' => (Element::S, pos + 1),
        'F' => (Element::F, pos + 1),
        'I' => (Element::I, pos + 1),
        _ => unreachable!(),
    }
}

// ── Bracket atom expression parser ───────────────────────────────────────────

/// Parse the expression inside `[…]`.  Returns `(expr, pos_of_closing_bracket)`.
fn parse_atom_expr(chars: &[char], start: usize) -> Result<(AtomExpr, usize), SmartsError> {
    parse_or_expr(chars, start)
}

/// `or_expr := and_low_expr (',' and_low_expr)*`
fn parse_or_expr(chars: &[char], pos: usize) -> Result<(AtomExpr, usize), SmartsError> {
    let (mut lhs, mut pos) = parse_and_low_expr(chars, pos)?;
    while pos < chars.len() && chars[pos] == ',' {
        pos += 1;
        let (rhs, new_pos) = parse_and_low_expr(chars, pos)?;
        lhs = AtomExpr::Or(Box::new(lhs), Box::new(rhs));
        pos = new_pos;
    }
    Ok((lhs, pos))
}

/// `and_low_expr := and_expr (';' and_expr)*`
fn parse_and_low_expr(chars: &[char], pos: usize) -> Result<(AtomExpr, usize), SmartsError> {
    let (mut lhs, mut pos) = parse_and_expr(chars, pos)?;
    while pos < chars.len() && chars[pos] == ';' {
        pos += 1;
        let (rhs, new_pos) = parse_and_expr(chars, pos)?;
        lhs = AtomExpr::And(Box::new(lhs), Box::new(rhs));
        pos = new_pos;
    }
    Ok((lhs, pos))
}

/// `and_expr := not_expr ('&' not_expr | <implicit-and> not_expr)*`
///
/// Implicit AND: adjacent primitives (not followed by `]`, `,`, `;`).
fn parse_and_expr(chars: &[char], pos: usize) -> Result<(AtomExpr, usize), SmartsError> {
    let (mut lhs, mut pos) = parse_not_expr(chars, pos)?;
    loop {
        if pos >= chars.len() || matches!(chars[pos], ']' | ',' | ';') {
            break;
        }
        if chars[pos] == '&' {
            pos += 1; // consume explicit '&'
        }
        // If the next character cannot start a primitive, stop.
        if pos >= chars.len() || matches!(chars[pos], ']' | ',' | ';' | '&') {
            break;
        }
        let (rhs, new_pos) = parse_not_expr(chars, pos)?;
        lhs = AtomExpr::And(Box::new(lhs), Box::new(rhs));
        pos = new_pos;
    }
    Ok((lhs, pos))
}

/// `not_expr := '!' not_expr | primitive`
fn parse_not_expr(chars: &[char], pos: usize) -> Result<(AtomExpr, usize), SmartsError> {
    if pos < chars.len() && chars[pos] == '!' {
        let (inner, new_pos) = parse_not_expr(chars, pos + 1)?;
        Ok((AtomExpr::Not(Box::new(inner)), new_pos))
    } else {
        parse_primitive(chars, pos)
    }
}

fn parse_primitive(chars: &[char], pos: usize) -> Result<(AtomExpr, usize), SmartsError> {
    if pos >= chars.len() || chars[pos] == ']' {
        return Err(SmartsError::UnexpectedChar(
            chars.get(pos).copied().unwrap_or(']'),
            pos,
        ));
    }

    match chars[pos] {
        // Atomic number: #n
        '#' => {
            let (num_str, new_pos) = collect_digits(chars, pos + 1);
            let atomic_num = num_str
                .parse::<u8>()
                .map_err(|_| SmartsError::InvalidAtomicNum(num_str))?;
            Ok((
                AtomExpr::Primitive(AtomPrimitive::AtomicNum(atomic_num)),
                new_pos,
            ))
        }

        // Aromatic/aliphatic class
        'a' => Ok((AtomExpr::Primitive(AtomPrimitive::Aromatic), pos + 1)),
        'A' => Ok((AtomExpr::Primitive(AtomPrimitive::Aliphatic), pos + 1)),

        // Wildcard
        '*' => Ok((AtomExpr::Primitive(AtomPrimitive::Any), pos + 1)),

        // H count: H (= 1) or Hn
        'H' => {
            let new_pos = pos + 1;
            if new_pos < chars.len() && chars[new_pos].is_ascii_digit() {
                let count = (chars[new_pos] as u8) - b'0';
                Ok((
                    AtomExpr::Primitive(AtomPrimitive::HCount(count)),
                    new_pos + 1,
                ))
            } else {
                Ok((AtomExpr::Primitive(AtomPrimitive::HCount(1)), new_pos))
            }
        }

        // Positive charge: + (= +1) or +n
        '+' => {
            let new_pos = pos + 1;
            if new_pos < chars.len() && chars[new_pos].is_ascii_digit() {
                let charge = (chars[new_pos] as i8) - b'0' as i8;
                Ok((
                    AtomExpr::Primitive(AtomPrimitive::Charge(charge)),
                    new_pos + 1,
                ))
            } else {
                Ok((AtomExpr::Primitive(AtomPrimitive::Charge(1)), new_pos))
            }
        }

        // Negative charge: - (= -1) or -n
        '-' => {
            let new_pos = pos + 1;
            if new_pos < chars.len() && chars[new_pos].is_ascii_digit() {
                let charge = -((chars[new_pos] as i8) - b'0' as i8);
                Ok((
                    AtomExpr::Primitive(AtomPrimitive::Charge(charge)),
                    new_pos + 1,
                ))
            } else {
                Ok((AtomExpr::Primitive(AtomPrimitive::Charge(-1)), new_pos))
            }
        }

        // Ring membership: R (any ring), R0 (not in ring), Rn (ring of size n)
        'R' => {
            let new_pos = pos + 1;
            if new_pos < chars.len() && chars[new_pos].is_ascii_digit() {
                let n = (chars[new_pos] as usize) - b'0' as usize;
                if n == 0 {
                    Ok((AtomExpr::Primitive(AtomPrimitive::NotInRing), new_pos + 1))
                } else {
                    Ok((
                        AtomExpr::Primitive(AtomPrimitive::Ring(Some(n))),
                        new_pos + 1,
                    ))
                }
            } else {
                Ok((AtomExpr::Primitive(AtomPrimitive::Ring(None)), new_pos))
            }
        }

        // Heavy-atom degree: Dn
        'D' => {
            let new_pos = pos + 1;
            if new_pos < chars.len() && chars[new_pos].is_ascii_digit() {
                let degree = (chars[new_pos] as u8) - b'0';
                Ok((
                    AtomExpr::Primitive(AtomPrimitive::Degree(degree)),
                    new_pos + 1,
                ))
            } else {
                // D without digit = any degree (treat as D >= 1, represented as Any)
                Ok((AtomExpr::Primitive(AtomPrimitive::Any), new_pos))
            }
        }

        // Element symbol (aliphatic if uppercase, aromatic if lowercase)
        c if c.is_ascii_uppercase() => parse_element_prim(chars, pos, false),
        c if c.is_ascii_lowercase() => parse_element_prim(chars, pos, true),

        c => Err(SmartsError::UnexpectedChar(c, pos)),
    }
}

/// Collect consecutive ASCII digits starting at `pos`.
fn collect_digits(chars: &[char], pos: usize) -> (String, usize) {
    let mut s = String::new();
    let mut p = pos;
    while p < chars.len() && chars[p].is_ascii_digit() {
        s.push(chars[p]);
        p += 1;
    }
    (s, p)
}

/// Parse an element symbol inside brackets (handles two-character symbols).
fn parse_element_prim(
    chars: &[char],
    pos: usize,
    aromatic: bool,
) -> Result<(AtomExpr, usize), SmartsError> {
    let upper_first = chars[pos].to_ascii_uppercase();

    // Try two-character symbol first (e.g. Cl, Br, Si, …)
    if pos + 1 < chars.len() {
        let second = chars[pos + 1];
        if second.is_ascii_lowercase() {
            let mut sym = String::new();
            sym.push(upper_first);
            sym.push(second);
            if let Some(elem) = Element::from_symbol(&sym) {
                let prim = if aromatic {
                    AtomPrimitive::ElementAromatic(elem)
                } else {
                    AtomPrimitive::ElementAliphatic(elem)
                };
                return Ok((AtomExpr::Primitive(prim), pos + 2));
            }
        }
    }

    // Single-character symbol
    let sym = upper_first.to_string();
    if let Some(elem) = Element::from_symbol(&sym) {
        let prim = if aromatic {
            AtomPrimitive::ElementAromatic(elem)
        } else {
            AtomPrimitive::ElementAliphatic(elem)
        };
        Ok((AtomExpr::Primitive(prim), pos + 1))
    } else {
        Err(SmartsError::UnexpectedChar(chars[pos], pos))
    }
}
