//! Subgraph isomorphism (VF2-style recursive backtracking) for SMARTS matching.

use super::ast::{AtomExpr, AtomPrimitive, SmartsBond, SmartsPattern};
use crate::mol::atom::Atom;
use crate::mol::bond::BondType;
use crate::mol::graph::MolGraph;
use petgraph::graph::NodeIndex;

/// Return `true` if the molecule contains at least one match for `pattern`.
pub fn has_match(mol: &MolGraph, pattern: &SmartsPattern) -> bool {
    if pattern.num_atoms() == 0 {
        return true; // empty pattern matches everything
    }
    let mol_atoms: Vec<NodeIndex> = mol.node_indices().collect();
    let pat_adj = build_pat_adj(pattern);
    let mut mapping = vec![None::<NodeIndex>; pattern.num_atoms()];
    let mut used = vec![false; mol_atoms.len()];
    search_any(
        mol,
        pattern,
        &mol_atoms,
        &pat_adj,
        0,
        &mut mapping,
        &mut used,
    )
}

/// Return all subgraph matches of `pattern` in `mol`.
///
/// Each match is a `Vec<NodeIndex>` of length `pattern.num_atoms()`,
/// where `result[i]` is the molecule atom that maps to pattern atom `i`.
pub fn subgraph_match(mol: &MolGraph, pattern: &SmartsPattern) -> Vec<Vec<NodeIndex>> {
    if pattern.num_atoms() == 0 {
        return vec![vec![]];
    }
    let mol_atoms: Vec<NodeIndex> = mol.node_indices().collect();
    if pattern.num_atoms() > mol_atoms.len() {
        return vec![];
    }
    let pat_adj = build_pat_adj(pattern);
    let mut results = Vec::new();
    let mut mapping = vec![None::<NodeIndex>; pattern.num_atoms()];
    let mut used = vec![false; mol_atoms.len()];
    search_all(
        mol,
        pattern,
        &mol_atoms,
        &pat_adj,
        0,
        &mut mapping,
        &mut used,
        &mut results,
    );
    results
}

// ── Adjacency list for pattern ────────────────────────────────────────────────

/// `pat_adj[i]` = list of `(neighbour_pat_idx, edge_idx_in_pattern.edges)`.
type PatAdj = Vec<Vec<(usize, usize)>>;

fn build_pat_adj(pattern: &SmartsPattern) -> PatAdj {
    let n = pattern.num_atoms();
    let mut adj: PatAdj = vec![vec![]; n];
    for (eidx, (a, b, _)) in pattern.edges.iter().enumerate() {
        adj[*a].push((*b, eidx));
        adj[*b].push((*a, eidx));
    }
    adj
}

// ── Recursive search (early-exit variant) ────────────────────────────────────

fn search_any(
    mol: &MolGraph,
    pattern: &SmartsPattern,
    mol_atoms: &[NodeIndex],
    pat_adj: &PatAdj,
    depth: usize,
    mapping: &mut Vec<Option<NodeIndex>>,
    used: &mut Vec<bool>,
) -> bool {
    if depth == pattern.num_atoms() {
        return true;
    }
    for (mol_pos, &mol_node) in mol_atoms.iter().enumerate() {
        if used[mol_pos] {
            continue;
        }
        if feasible(mol, pattern, pat_adj, mapping, depth, mol_node) {
            mapping[depth] = Some(mol_node);
            used[mol_pos] = true;
            if search_any(mol, pattern, mol_atoms, pat_adj, depth + 1, mapping, used) {
                return true;
            }
            mapping[depth] = None;
            used[mol_pos] = false;
        }
    }
    false
}

// ── Recursive search (collect-all variant) ───────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn search_all(
    mol: &MolGraph,
    pattern: &SmartsPattern,
    mol_atoms: &[NodeIndex],
    pat_adj: &PatAdj,
    depth: usize,
    mapping: &mut Vec<Option<NodeIndex>>,
    used: &mut Vec<bool>,
    results: &mut Vec<Vec<NodeIndex>>,
) {
    if depth == pattern.num_atoms() {
        results.push(mapping.iter().map(|m| m.unwrap()).collect());
        return;
    }
    for (mol_pos, &mol_node) in mol_atoms.iter().enumerate() {
        if used[mol_pos] {
            continue;
        }
        if feasible(mol, pattern, pat_adj, mapping, depth, mol_node) {
            mapping[depth] = Some(mol_node);
            used[mol_pos] = true;
            search_all(
                mol,
                pattern,
                mol_atoms,
                pat_adj,
                depth + 1,
                mapping,
                used,
                results,
            );
            mapping[depth] = None;
            used[mol_pos] = false;
        }
    }
}

// ── Feasibility check ────────────────────────────────────────────────────────

fn feasible(
    mol: &MolGraph,
    pattern: &SmartsPattern,
    pat_adj: &PatAdj,
    mapping: &[Option<NodeIndex>],
    pat_idx: usize,
    mol_node: NodeIndex,
) -> bool {
    let atom = &mol[mol_node];

    // 1. Atom predicate
    if !matches_atom_expr(&pattern.atoms[pat_idx].expr, atom, mol, mol_node) {
        return false;
    }

    // 2. Bond consistency with already-mapped neighbours
    for &(pat_nbr, edge_idx) in &pat_adj[pat_idx] {
        if pat_nbr >= mapping.len() {
            continue;
        }
        if let Some(mol_nbr) = mapping[pat_nbr] {
            match mol.find_edge(mol_node, mol_nbr) {
                None => return false,
                Some(eid) => {
                    if !matches_bond(&pattern.edges[edge_idx].2, mol[eid]) {
                        return false;
                    }
                }
            }
        }
    }

    true
}

// ── Atom-expression matching ──────────────────────────────────────────────────

fn matches_atom_expr(expr: &AtomExpr, atom: &Atom, mol: &MolGraph, node: NodeIndex) -> bool {
    match expr {
        AtomExpr::Primitive(p) => matches_primitive(p, atom, mol, node),
        AtomExpr::Not(inner) => !matches_atom_expr(inner, atom, mol, node),
        AtomExpr::And(a, b) => {
            matches_atom_expr(a, atom, mol, node) && matches_atom_expr(b, atom, mol, node)
        }
        AtomExpr::Or(a, b) => {
            matches_atom_expr(a, atom, mol, node) || matches_atom_expr(b, atom, mol, node)
        }
    }
}

fn matches_primitive(prim: &AtomPrimitive, atom: &Atom, mol: &MolGraph, node: NodeIndex) -> bool {
    match prim {
        AtomPrimitive::Any => true,
        AtomPrimitive::AtomicNum(n) => atom.element.atomic_number() == *n,
        AtomPrimitive::ElementAliphatic(elem) => atom.element == *elem && !atom.aromatic,
        AtomPrimitive::ElementAromatic(elem) => atom.element == *elem && atom.aromatic,
        AtomPrimitive::Aromatic => atom.aromatic,
        AtomPrimitive::Aliphatic => !atom.aromatic,
        AtomPrimitive::HCount(n) => atom.h_count == *n,
        AtomPrimitive::Charge(c) => atom.charge == *c,
        AtomPrimitive::Ring(None) => atom.in_ring,
        AtomPrimitive::Ring(Some(size)) => atom.ring_sizes.contains(&(*size as u8)),
        AtomPrimitive::NotInRing => !atom.in_ring,
        AtomPrimitive::Degree(n) => mol.neighbors(node).count() == *n as usize,
    }
}

// ── Bond matching ─────────────────────────────────────────────────────────────

fn matches_bond(pat: &SmartsBond, mol_bond: BondType) -> bool {
    match pat {
        SmartsBond::Any => true,
        SmartsBond::Unspecified => matches!(mol_bond, BondType::Single | BondType::Aromatic),
        SmartsBond::Single => mol_bond == BondType::Single,
        SmartsBond::Double => mol_bond == BondType::Double,
        SmartsBond::Triple => mol_bond == BondType::Triple,
        SmartsBond::Aromatic => mol_bond == BondType::Aromatic,
    }
}
