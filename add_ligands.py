#!/usr/bin/env python3
"""
add_ligands.py
---------------
Combine multiple ligand files (SDF/MOL/MOL2/PDB or PDB IDs)
into a single multi-ligand SDF file with hydrogens and 3D coordinates.

Usage:
    python3 add_ligands.py output.sdf ligand1.sdf ligand2.sdf ...
"""

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
from io import StringIO

def fetch_from_pdb(ligand_id):
    """Fetch a ligand from RCSB by its three-letter PDB ID."""
    url = f"https://files.rcsb.org/ligands/view/{ligand_id.upper()}_model.sdf"
    print(f"📥 Fetching ligand {ligand_id.upper()} from RCSB...")
    r = requests.get(url)
    if r.status_code != 200:
        print(f"⚠️  Could not fetch {ligand_id.upper()} from RCSB (HTTP {r.status_code}).")
        return None
    suppl = Chem.SDMolSupplier(StringIO(r.text))
    mols = [m for m in suppl if m is not None]
    return mols

def read_ligand_file(path):
    """Read a ligand file in SDF/MOL/MOL2/PDB format."""
    ext = os.path.splitext(path)[1].lower()
    mols = []
    try:
        if ext == ".sdf":
            suppl = Chem.SDMolSupplier(path)
            mols = [m for m in suppl if m is not None]
        elif ext in (".mol", ".mol2", ".pdb"):
            m = Chem.MolFromMolFile(path, removeHs=False)
            if m:
                mols = [m]
        else:
            print(f"⚠️  Unsupported format for {path}. Skipping.")
    except Exception as e:
        print(f"⚠️  Error reading {path}: {e}")
    return mols

def prepare_mol(mol):
    """Add hydrogens and generate 3D coordinates."""
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        print(f"⚠️  3D generation failed: {e}")
    return mol

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 add_ligands.py output.sdf ligand1.sdf ligand2.sdf ...")
        sys.exit(1)

    out_sdf = sys.argv[1]
    ligand_sources = sys.argv[2:]

    # Create the output file if it doesn't exist
    if not os.path.exists(out_sdf):
        print(f"🪄 Creating new output file: {out_sdf}")
        open(out_sdf, "w").close()

    writer = Chem.SDWriter(out_sdf)
    total = 0

    for src in ligand_sources:
        if len(src) == 3 and src.isalpha():
            mols = fetch_from_pdb(src)
        else:
            mols = read_ligand_file(src)

        if not mols:
            print(f"⚠️  No valid molecules found in {src}")
            continue

        for mol in mols:
            mol = prepare_mol(mol)
            name = os.path.splitext(os.path.basename(src))[0]
            mol.SetProp("_Name", name)
            writer.write(mol)
            total += 1

    writer.close()
    print(f"\n✅ Combined {total} ligands into {out_sdf}")

if __name__ == "__main__":
    main()

