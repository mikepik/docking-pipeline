#!/usr/bin/env python3
import os
import csv
import re
import subprocess
import sys

def parse_vina_output(pdbqt_path):
    """Extract the best binding energy (kcal/mol) from a Vina .pdbqt file."""
    try:
        with open(pdbqt_path, "r") as f:
            for line in f:
                if "REMARK VINA RESULT" in line:
                    parts = line.split()
                    # Example: ['REMARK', 'VINA', 'RESULT:', '-5.830', '0.000', '0.000']
                    if len(parts) >= 4:
                        return float(parts[3])
    except Exception:
        return None
    return None

def main(results_dir):
    scores = []

    for file in os.listdir(results_dir):
        if file.endswith("_out.pdbqt"):
            ligand_name = file.replace("_out.pdbqt", "")
            pdbqt_path = os.path.join(results_dir, file)
            score = parse_vina_output(pdbqt_path)
            if score is not None:
                scores.append((ligand_name, score))

    if not scores:
        print("⚠️  No valid docking result files or scores found.")
        sys.exit(1)

    # Sort best (lowest energy) → worst (highest)
    scores.sort(key=lambda x: x[1])

    csv_path = os.path.join(results_dir, "binding_scores.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Ligand", "Best ΔG (kcal/mol)"])
        writer.writerows(scores)

    print(f"✅ CSV written to: {csv_path}")
    print("📊 Top 5 ligands:")
    for lig, sc in scores[:5]:
        print(f"  {lig:30s} {sc:8.2f}")

    # Try to auto-open the CSV file
    try:
        if sys.platform.startswith("linux"):
            subprocess.run(["xdg-open", csv_path])
        elif sys.platform == "darwin":
            subprocess.run(["open", csv_path])
        elif sys.platform.startswith("win"):
            os.startfile(csv_path)
    except Exception as e:
        print(f"⚠️  Could not auto-open CSV: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 summarize_vina_scores.py <results_dir>")
        sys.exit(1)
    results_dir = sys.argv[1]
    print("🚀 Summarizing docking results...")
    main(results_dir)

