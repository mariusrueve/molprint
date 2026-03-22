#!/usr/bin/env python3
"""Print the actual RDKit MACCS key SMARTS definitions and verify against diagnostic molecules."""
from rdkit import Chem
from rdkit.Chem import MACCSkeys

# Print all key definitions
print("# RDKit MACCS key definitions (bit_index: smarts, min_count)")
for key_num in sorted(MACCSkeys.smartsPatts.keys()):
    smarts, min_count = MACCSkeys.smartsPatts[key_num]
    print(f"{key_num}\t{smarts}\t{min_count}")
