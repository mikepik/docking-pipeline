#!/usr/bin/env python3
"""
convert_to_pdbqt.py
-------------------
Split a multi-ligand SDF into individual ligands and convert each to PDBQT.

Usage:
    python3 convert_to_pdbqt.py combined.sdf output_folder/

Requires:
    Open Babel (obabel command in PATH)

Example:
    python3 convert_to_pdbqt.py combined_20251009_142021.sdf pdbqt_out/
"""

import os
import sys
import subprocess
from rdkit import Chem

def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_to_pdbqt.py input.sdf output_folder/")
        sys.exit(1)

    input_sdf = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(input_sdf):
        print(f"❌ Error: file not found: {input_sdf}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Count and split ligands using RDKit
    suppl = Chem.SDMolSupplier(input_sdf)
    mols = [m for m in suppl if m is not None]

    print(f"✅ Found {len(mols)} molecules in {input_sdf}")
    sdf_list = []

    for i, mol in enumerate(mols, start=1):
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{i}"
        safe_name = "".join(c if c.isalnum() or c in ("_", "-") else "_" for c in name)
        sdf_file = os.path.join(output_dir, f"{safe_name}.sdf")
        pdbqt_file = os.path.join(output_dir, f"{safe_name}.pdbqt")
        sdf_list.append((sdf_file, pdbqt_file))

        # Write each ligand to its own SDF
        writer = Chem.SDWriter(sdf_file)
        writer.write(mol)
        writer.close()

    print(f"🧩 Split into {len(sdf_list)} individual SDF files.")

    # Convert each SDF → PDBQT using Open Babel
    converted = 0
    for sdf_file, pdbqt_file in sdf_list:
        cmd = ["obabel", sdf_file, "-O", pdbqt_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            converted += 1
        else:
            print(f"⚠️  Conversion failed for {sdf_file}:\n{result.stderr}")

    print(f"✨ Successfully converted {converted} ligands to PDBQT in '{output_dir}'")

if __name__ == "__main__":
    main()

