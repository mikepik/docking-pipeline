#!/usr/bin/env python3
import os
import sys
import subprocess
from time import time

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 run_vina_batch.py receptor.pdbqt ligands_dir results_dir config.txt")
        sys.exit(1)

    receptor = sys.argv[1]
    ligands_dir = sys.argv[2]
    results_dir = sys.argv[3]
    config_file = sys.argv[4]

    os.makedirs(results_dir, exist_ok=True)
    ligands = [os.path.join(ligands_dir, f) for f in os.listdir(ligands_dir) if f.endswith(".pdbqt")]

    total = len(ligands)
    print(f"\n🧩 Docking {total} ligands using AutoDock Vina...\n")

    start_time = time()

    for i, ligand_file in enumerate(ligands):
        ligand_name = os.path.basename(ligand_file)
        out_file = os.path.join(results_dir, ligand_name.replace(".pdbqt", "_out.pdbqt"))

        cmd = [
            "vina",
            "--receptor", receptor,
            "--ligand", ligand_file,
            "--config", config_file,
            "--out", out_file
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"⚠️  Vina failed for {ligand_name}: {e}")

        elapsed = time() - start_time
        avg_time = elapsed / (i + 1)
        est_remaining = avg_time * (total - (i + 1))
        print(f"⏱️  Completed {i + 1}/{total} ligands "
              f"(avg {avg_time:.1f}s/ligand, est {est_remaining/60:.1f} min remaining)\n")

    total_time = (time() - start_time) / 60
    print(f"✅ Docking completed for {total} ligands in {total_time:.1f} minutes.")
    print(f"✅ Results stored in: {results_dir}")

if __name__ == "__main__":
    main()

