#!/usr/bin/env python3
import os, sys, subprocess, argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Batch runner for AutoDock Vina over a directory of ligand PDBQT files."
    )
    p.add_argument("receptor", help="Path to receptor .pdbqt")
    p.add_argument("ligands_dir", help="Directory containing ligand .pdbqt files")
    p.add_argument("results_dir", help="Output directory for docked poses")
    p.add_argument("config", help="Vina config file (e.g., vina_config.txt)")
    p.add_argument("--cpu", type=int, default=8, help="Number of CPU cores for Vina (default: 8)")
    return p.parse_args()

def fail(msg, code=2):
    print(f"❌ {msg}")
    sys.exit(code)

def main():
    args = parse_args()

    # Validate inputs
    if not os.path.isfile(args.receptor):
        fail(f"Receptor not found: {args.receptor}")
    if not os.path.isfile(args.config):
        fail(f"Vina config not found: {args.config}")
    if not os.path.isdir(args.ligands_dir):
        fail(f"Ligands directory not found: {args.ligands_dir}")

    ligands = [f for f in os.listdir(args.ligands_dir) if f.endswith(".pdbqt")]
    ligands.sort()
    if not ligands:
        fail(f"No .pdbqt ligands found in {args.ligands_dir}")

    os.makedirs(args.results_dir, exist_ok=True)

    total = len(ligands)
    completed = 0
    for i, lig in enumerate(ligands, 1):
        lig_path = os.path.join(args.ligands_dir, lig)
        out_name = os.path.splitext(lig)[0] + "_out.pdbqt"
        out_path = os.path.join(args.results_dir, out_name)

        cmd = [
            "vina",
            "--receptor", args.receptor,
            "--ligand", lig_path,
            "--config", args.config,
            "--out", out_path,
            "--cpu", str(args.cpu),
        ]

        try:
            subprocess.run(cmd, check=True)
            completed += 1
        except subprocess.CalledProcessError as e:
            print(f"⚠️  Vina failed for {lig}: {e}")

        print(f"⏱️  Completed {i}/{total} ligands")

    print(f"\n✅ Docking completed for {completed}/{total} ligands.")
    print(f"✅ Results stored in: {args.results_dir}")

if __name__ == "__main__":
    main()

