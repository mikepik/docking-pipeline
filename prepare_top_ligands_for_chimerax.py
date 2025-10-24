#!/usr/bin/env python3
import csv, os, sys, shutil

def main(receptor, results_dir="results/", top_n=5):
    results_dir = os.path.abspath(results_dir)
    csv_path = os.path.join(results_dir, "binding_scores.csv")
    top_dir = os.path.join(results_dir, "top_ligands")
    os.makedirs(top_dir, exist_ok=True)

    receptor = os.path.abspath(receptor)

    if not os.path.exists(receptor):
        print(f"❌ Receptor file not found: {receptor}")
        sys.exit(1)
    if not os.path.exists(csv_path):
        print(f"❌ No binding_scores.csv found in {results_dir}")
        sys.exit(1)

    # Read and sort binding scores
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        rows = sorted(reader, key=lambda r: float(r["Best ΔG (kcal/mol)"]))

    top_rows = rows[:top_n]

    print(f"📦 Preparing top {top_n} ligands for ChimeraX display...")
    cxc_path = os.path.join(results_dir, "view_top_ligands.cxc")
    with open(cxc_path, "w") as cxc:
        cxc.write(f"open {receptor}\n")
        for row in top_rows:
            ligand_name = row["Ligand"]
            pdbqt_file = os.path.join(results_dir, f"{ligand_name}_out.pdbqt")
            if os.path.exists(pdbqt_file):
                dest = os.path.join(top_dir, f"{ligand_name}_out.pdbqt")
                shutil.copy(pdbqt_file, dest)
                abs_dest = os.path.abspath(dest)
                cxc.write(f"open {abs_dest}\n")
                print(f"  ✅ {ligand_name}")
            else:
                print(f"  ⚠️ Missing docked file for {ligand_name}")
        cxc.write("preset apply pub 3\n")
        cxc.write("view\n")

    abs_cxc = os.path.abspath(cxc_path)
    print(f"\n✅ Created ChimeraX command file: {abs_cxc}")
    print(f"💡 To open in ChimeraX, copy and paste this command:")
    print(f"\n   chimerax @{abs_cxc}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 prepare_top_ligands_for_chimerax.py <receptor.pdbqt> [results_dir] [top_n]")
        sys.exit(1)

    receptor = sys.argv[1]
    # detect whether the 2nd argument is a directory name or a number
    results_dir = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].isdigit() else "results/"
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else (
        int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 5
    )

    main(receptor, results_dir, top_n)




