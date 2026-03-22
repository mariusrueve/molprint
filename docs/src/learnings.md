# Learnings

These are honest technical insights from building molprint — things that took longer than expected, bugs that required careful thought to fix, and places where the problem was harder than it looked.

## RDKit MACCS semantics: uniquify=True and min_count

The MACCS keys look straightforward on paper — 166 SMARTS patterns, one bit per pattern. The bit is set if the pattern has at least one match in the molecule.

But "at least one match" hides complexity. RDKit's `GetMACCSKeysFingerprint` uses `uniquify=True` in its SMARTS matching, which means each atom can only be used once per match. For some keys, the bit is set only if there are at least N matches (some keys test for "more than 1 ring" or "more than 2 atoms of type X").

The practical consequence: if you implement MACCS as "compile pattern, call `has_match`, set bit", you'll get the wrong answer for about 10–15 keys on any reasonably complex molecule. The correct implementation requires knowing the intended count semantics for each key, which comes from reading the RDKit source code for `MACCSkeys.py` directly.

Example: MACCS key 166 (`?`) is always 0 in RDKit's implementation. Key 0 is also always 0 (it's a placeholder). Several keys in the 140–150 range test for the presence of multiple ring systems or specific ring counts, not just "does this pattern appear at all".

Running the validation suite (`validate_against_rdkit.rs`) against a ChEMBL subset was the only reliable way to catch all of these edge cases.

## Aromatic implicit hydrogen calculation

Computing implicit H counts for aromatic atoms is the most subtle part of the SMILES parser.

Consider pyridine (`c1ccncc1`): the nitrogen is aromatic. It has two aromatic bonds (bond order sum = 2 from the sigma bonds). Its normal valence is 3. So naively: `3 - 2 = 1 H`. But pyridine nitrogen has **zero** hydrogens.

The issue: aromatic bonds each contribute 1 to the sigma bond sum, but the pi system takes up the remaining valence capacity. The standard formula adds +1 to the effective bond sum to account for the pi electron: `effective_sum = sigma_sum + 1 = 3`, then `h = 3 - 3 = 0`. Correct.

Now consider pyrrole nitrogen (`[nH]`): it contributes 2 pi electrons to the aromatic system (it's the nitrogen with the lone pair). Its sigma bonds are 2 (two ring C-N bonds). Effective sum = 3, valence = 3, h = 0. But we wrote `[nH]` explicitly — the bracket forces H=1.

The problem emerges for organic-subset aromatic nitrogen without explicit bracket notation. In practice, most aromatic nitrogens in the organic subset that need an H are written as `[nH]` in valid SMILES. But the implicit H formula still needs to be correct for cases like `n` in a ring where no H is present.

The final rule: apply the +1 pi adjustment only when the atom has available valence remaining after its sigma bonds. If the sigma bonds already consume all valence capacity (bond_sum == some_valence), don't add +1 — the atom is already fully bonded and gets 0 implicit H. This handles the trimethylamine-oxide case and similar edge cases.

## Spiro atom ring-bond exclusion bugs

The first SSSR implementation had a subtle bug with spiro compounds. A spiro compound has two rings sharing exactly one atom (the spiro center). The SSSR should report two rings.

An earlier version of the candidate selection excluded edges whose endpoints were already both "in a ring". The reasoning was: if both atoms are already in a ring, the edge is redundant. This is wrong for spiro compounds — the spiro center is in one ring, but the edges connecting it to the second ring are not redundant.

The fix was to generate candidates based purely on BFS reachability (for each edge, can you get from one endpoint to the other without using that edge?) and rely entirely on the GF(2) independence test to filter out duplicates. This is slower in theory (more candidates) but correct for all cases including spiro, fused, and bridged ring systems.

The specific test that caught this: `C1CCC2(CC1)CCCC2` (spiro[5.5]undecane) should have 2 rings. The broken version reported 1.

## SSSR vs. DFS cycle enumeration

An early temptation was to use DFS to find all simple cycles and then deduplicate. This is conceptually simple but produces too many rings for fused systems.

For naphthalene (`c1ccc2ccccc2c1`), DFS enumeration finds 3 simple cycles: the two 6-membered rings and the 10-membered outer ring. The SSSR finds 2 (just the two 6-membered rings). RDKit uses SSSR, so molprint must match.

The GF(2) basis selection elegantly handles this: the 10-membered outer ring of naphthalene is linearly dependent on the two 6-membered rings in the edge-vector basis, so it is automatically excluded.

## Hash function must be deterministic

During early development, the Morgan fingerprint implementation used `HashMap` for collecting environment hashes. `HashMap` in Rust uses `AHasher` by default since Rust 1.36, and `AHasher::default()` uses a process-level random seed for hash-flooding resistance.

This made fingerprints non-deterministic across processes. A database built in one run couldn't be searched in another because the bit positions of all hash-folded environments were different each time.

The fix was to replace all uses of `HashSet<u32>` (for deduplication) with `BTreeSet` or, better, to use a custom hash function with a fixed seed. The hash functions in `morgan.rs` now use FNV-1a constants and the golden-ratio-based combine from boost::hash_combine — both fully deterministic.

The rule: never use `HashMap`/`HashSet`/`AHasher` in any code path that contributes to fingerprint output. Only use them for intermediate data structures where the output ordering doesn't matter (and then only where you've verified it doesn't affect the final bit vector).

## FPS byte ordering: least-significant byte first

The FPS format stores fingerprints as hexadecimal with bytes in little-endian order (least-significant byte first). This is chemfp's convention, and getting it wrong means fingerprints written by molprint can't be read by chemfp, and vice versa.

For a 64-bit fingerprint with only bit 0 set:
- The internal `u64` value is `0x0000000000000001`
- The little-endian byte representation is `01 00 00 00 00 00 00 00`
- The hex encoding is `0100000000000000`

A naive implementation might write the hex of the `u64` as a big-endian integer: `0000000000000001`. That's wrong for FPS interoperability.

The `to_hex` and `from_hex` methods in `FingerprintBits` iterate over bytes in little-endian order explicitly to get this right.

## The `pending_bond = None` on CloseBranch

The SMILES parser has a subtlety: when a `CloseBranch` token is encountered, the `pending_bond` must be cleared. Consider `CC(=O)C`:

- `C` — atom 0
- `C` — atom 1, bonded to atom 0
- `(` — push atom 1 onto branch stack
- `=` — set pending_bond = Double
- `O` — atom 2, bonded to atom 1 with Double bond
- `)` — pop branch stack, restore current = atom 1, **clear pending_bond**
- `C` — atom 3, bonded to atom 1 with Single bond (not Double!)

Without clearing `pending_bond` on `)`, the final carbon would incorrectly receive a double bond. This was caught by the acetic acid test (`CC(=O)O`).
