#!/usr/bin/env python3
"""
run_pipeline.py
---------------
Fully automated AutoDock Vina docking workflow.

Usage:
    python3 run_pipeline.py docking_config.yml
"""

import os
import sys
import subprocess
import yaml
import requests
import shutil
import datetime

print(f"🧬 Using Python from: {sys.executable}")

# ---------------- Utility ---------------- #

def run_cmd(command, desc):
    print(f"\n🚀 {desc}...")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"✅ Done: {desc}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during: {desc}\n{e}")
        sys.exit(1)

# ---------------- Step 0: Receptor prep ---------------- #

def prepare_receptor(pdb_id, chains=None, keep_cofactors=False):
    prepped_file = f"{pdb_id}_prepped.pdbqt"

    # ✅ Skip download if receptor already exists
    if os.path.exists(prepped_file):
        print(f"📦 Found existing {prepped_file}, skipping download and prep.")
        return prepped_file

    pdb_file = f"{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    print(f"📥 Downloading {pdb_id} from RCSB...")
    r = requests.get(url)
    if r.status_code != 200:
        sys.exit(f"❌ Failed to download {pdb_id} (HTTP {r.status_code})")

    with open(pdb_file, "w") as f:
        f.write(r.text)
    print(f"✅ Downloaded {pdb_file}")

    # Clean the file
    clean_file = f"{pdb_id}_clean.pdb"
    print(f"🧹 Cleaning {pdb_file} → {clean_file}")

    with open(pdb_file) as f:
        lines = f.readlines()

    keep_chains = set(chains.split(",")) if chains else None
    cleaned = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        resname = line[17:20].strip()
        chain = line[21].strip()

        # Skip solvent
        if resname in {"HOH", "WAT", "SOL"}:
            continue
        # Skip cofactors unless requested
        if not keep_cofactors and line.startswith("HETATM"):
            continue
        # Keep selected chains only
        if keep_chains and chain not in keep_chains:
            continue
        # Keep only main conformation
        alt = line[16].strip()
        if alt not in ("", "A"):
            continue
        # Normalize MSE -> MET
        if resname == "MSE":
            line = line[:17] + "MET" + line[20:]
        cleaned.append(line)

    with open(clean_file, "w") as f:
        f.writelines(cleaned)
        f.write("END\n")
    print(f"✅ Cleaned receptor saved: {clean_file}")

    # Add hydrogens and charges
    prepped_file = f"{pdb_id}_prepped.pdbqt"
    print(f"🔄 Adding hydrogens and charges via Open Babel...")
    cmd = f"obabel {clean_file} -O {prepped_file} -xh -p 7.4 --partialcharge gasteiger"
    subprocess.run(cmd, shell=True, check=True)
    print(f"✅ Receptor prepared: {prepped_file}")

    # --- Cleanup: remove ligand-style torsion tags ---
    tmp = prepped_file + ".tmp"
    with open(prepped_file) as fin, open(tmp, "w") as fout:
        for ln in fin:
            if not ln.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                fout.write(ln)
    os.replace(tmp, prepped_file)
    print("🧹 Removed ligand-style torsion tags from receptor file.")

    return prepped_file

# ---------------- Main pipeline ---------------- #

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 run_pipeline.py docking_config.yml")
        sys.exit(1)

    config_file = sys.argv[1]
    if not os.path.exists(config_file):
        sys.exit(f"❌ Config file not found: {config_file}")

    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    # ✅ Make key directories absolute before changing into the run directory
    cfg["ligand_dir"] = os.path.abspath(cfg["ligand_dir"])
    cfg["pdbqt_dir"] = os.path.abspath(cfg["pdbqt_dir"])
    cfg["results_dir"] = os.path.abspath(cfg["results_dir"])

    # ✅ Create a timestamped subdirectory for this run
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = os.path.join(os.getcwd(), f"docking_run_{timestamp}")
    src_config = os.path.abspath(config_file)
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    print(f"📁 Working directory: {run_dir}")
    
        # ✅ Redirect all pipeline outputs into this run's folders
    pdbqt_dir = os.path.join(run_dir, "pdbqt_files")
    results_dir = os.path.join(run_dir, "results")
    os.makedirs(pdbqt_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # Keep ligand_dir absolute (ligands are shared and not rewritten)
    ligand_dir = os.path.abspath(cfg["ligand_dir"])

    # ✅ Copy the YAML config into the run directory for reproducibility
    shutil.copy(src_config, os.path.join(run_dir, "config_used.yml"))
    print(f"🧾 Copied configuration to: {run_dir}/config_used.yml")


    # ✅ Copy the YAML config into the run directory for reproducibility
    shutil.copy(src_config, os.path.join(run_dir, "config_used.yml"))
    print(f"🧾 Copied configuration to: {run_dir}/config_used.yml")

    # ✅ Define helper script paths (absolute)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    add_ligands_script = os.path.join(base_dir, "add_ligands.py")
    convert_script = os.path.join(base_dir, "convert_to_pdbqt.py")
    box_script = os.path.join(base_dir, "make_vina_box_from_residues.py")
    vina_script = os.path.join(base_dir, "run_vina_batch.py")
    summary_script = os.path.join(base_dir, "summarize_vina_scores.py")

    # Load config values
    pdb_id = cfg.get("pdb_id", "")
    chains = cfg.get("chains", "")
    keep_cofactors = cfg.get("keep_cofactors", False)
    ligand_dir = cfg["ligand_dir"]
    combined_sdf = cfg["combined_sdf"]
    residues = cfg["residues"]
    padding = cfg["padding"]
    config_out = cfg["config_out"]
    cpu = cfg.get("cpu", 8)
    top_ligands = cfg.get("top_ligands", 5)

    # Ensure directories exist
    for d in [ligand_dir, pdbqt_dir, results_dir]:
        os.makedirs(d, exist_ok=True)

    # Step 0: Receptor prep
    receptor = cfg.get("receptor")
    if receptor and os.path.exists(receptor):
        print(f"🧬 Using pre-prepared receptor: {receptor}")
    else:
        receptor = prepare_receptor(pdb_id, chains, keep_cofactors)

    # Step 1: Combine ligands
    ligand_files_list = [
        os.path.join(ligand_dir, f)
        for f in os.listdir(ligand_dir)
        if f.endswith((".sdf", ".mol", ".mol2"))
    ]
    ligand_files = " ".join(ligand_files_list)
    cmd1 = f"python3 {add_ligands_script} {combined_sdf} {ligand_files}"
    run_cmd(cmd1, "Combining ligand files")

    # Step 2: Convert ligands to PDBQT
    cmd2 = f"python3 {convert_script} {combined_sdf} {pdbqt_dir}"
    run_cmd(cmd2, "Converting ligands to PDBQT")

    # Step 3: Create docking box
    cmd3 = f"python3 {box_script} {receptor} \"{residues}\" --padding {padding} --out {config_out}"
    run_cmd(cmd3, "Generating vina_config.txt")

    # Step 4: Recenter ligands
    print("\n🧭 Recentering ligand coordinates to docking box center...")

    center_x = center_y = center_z = None
    if os.path.exists(cfg["config_out"]):
        with open(cfg["config_out"]) as f:
            for line in f:
                if "center_x" in line:
                    center_x = float(line.split("=")[1].strip())
                elif "center_y" in line:
                    center_y = float(line.split("=")[1].strip())
                elif "center_z" in line:
                    center_z = float(line.split("=")[1].strip())

    if None in (center_x, center_y, center_z):
        print("⚠️  Could not read box center; skipping ligand recentering.")
    else:
        for lig in os.listdir(pdbqt_dir):
            if not lig.endswith(".pdbqt"):
                continue
            lig_path = os.path.join(pdbqt_dir, lig)
            with open(lig_path) as f:
                lines = f.readlines()

            xs, ys, zs = [], [], []
            for ln in lines:
                if ln.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(ln[30:38]))
                        ys.append(float(ln[38:46]))
                        zs.append(float(ln[46:54]))
                    except ValueError:
                        pass

            if not xs:
                continue

            cx0 = sum(xs) / len(xs)
            cy0 = sum(ys) / len(ys)
            cz0 = sum(zs) / len(zs)
            dx, dy, dz = center_x - cx0, center_y - cy0, center_z - cz0

            out = []
            for ln in lines:
                if ln.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(ln[30:38]) + dx
                        y = float(ln[38:46]) + dy
                        z = float(ln[46:54]) + dz
                        ln = f"{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}"
                    except ValueError:
                        pass
                out.append(ln)

            with open(lig_path, "w") as f:
                f.writelines(out)

        print("✅ Ligands recentered to box center.")

    # Step 5: Run docking
    cmd4 = f"python3 run_vina_batch.py {receptor} {pdbqt_dir} {results_dir} {config_out} --cpu {cpu}"
    run_cmd(cmd4, f"Running AutoDock Vina batch docking with {cpu} CPUs")

    # Step 6: Summarize results
    cmd5 = f"python3 {summary_script} {results_dir}"
    run_cmd(cmd5, "Summarizing docking results")

    print("\n🎉 Pipeline completed successfully!")
    print(f"Results stored in: {results_dir}\n")
    
    # Step 7: Prepare top ligands for ChimeraX visualization
    prepare_top_script = os.path.join(base_dir, "prepare_top_ligands_for_chimerax.py")
    cmd6 = f"python3 {prepare_top_script} {receptor} {results_dir} {top_ligands}"
    run_cmd(cmd6, f"Preparing top {top_ligands} ligands for ChimeraX")


if __name__ == "__main__":
    main()

