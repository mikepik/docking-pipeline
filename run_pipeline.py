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
import shutil
import datetime
import subprocess
import yaml
import requests

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


def resolve_path(base_dir, path_value):
    """Resolve a possibly-relative path against base_dir."""
    if not path_value:
        return path_value
    if os.path.isabs(path_value):
        return path_value
    return os.path.abspath(os.path.join(base_dir, path_value))


def write_vina_config(config_path, receptor, center_x, center_y, center_z, size_x, size_y, size_z):
    """Write vina_config.txt directly from box coordinates."""
    with open(config_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"center_x = {center_x}\n")
        f.write(f"center_y = {center_y}\n")
        f.write(f"center_z = {center_z}\n")
        f.write(f"size_x = {size_x}\n")
        f.write(f"size_y = {size_y}\n")
        f.write(f"size_z = {size_z}\n")
        f.write("exhaustiveness = 8\n")
        f.write("num_modes = 9\n")
        f.write("energy_range = 3\n")


# ---------------- Step 0: Receptor prep ---------------- #

def prepare_receptor(pdb_id, chains=None, keep_cofactors=False):
    """
    Download receptor from RCSB, clean it, and convert to PDBQT.
    Works inside the current run directory.
    """
    prepped_file = f"{pdb_id}_prepped.pdbqt"

    # Skip if already present in the current working directory
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

    # Convert to PDBQT
    print("🔄 Adding hydrogens and charges via Open Babel...")
    cmd = f"obabel {clean_file} -O {prepped_file} -xh -p 7.4 --partialcharge gasteiger"
    subprocess.run(cmd, shell=True, check=True)
    print(f"✅ Receptor prepared: {prepped_file}")

    # Remove ligand-style torsion tags from receptor PDBQT
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

    config_file = os.path.abspath(sys.argv[1])
    if not os.path.exists(config_file):
        sys.exit(f"❌ Config file not found: {config_file}")

    config_dir = os.path.dirname(config_file)
    base_dir = os.path.dirname(os.path.abspath(__file__))

    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    # ---------------- Config values ---------------- #

    pdb_id = cfg.get("pdb_id", "")
    receptor_cfg = cfg.get("receptor", "")
    chains = cfg.get("chains", "")
    keep_cofactors = cfg.get("keep_cofactors", False)

    ligand_dir = resolve_path(config_dir, cfg["ligand_dir"])
    combined_sdf = cfg["combined_sdf"]

    # These are subdirectory/file names used INSIDE the run folder
    pdbqt_subdir = cfg.get("pdbqt_dir", "pdbqt_files/")
    results_subdir = cfg.get("results_dir", "results/")
    config_out = cfg.get("config_out", "vina_config.txt")

    center_x = float(cfg["center_x"])
    center_y = float(cfg["center_y"])
    center_z = float(cfg["center_z"])

    size_x = float(cfg["size_x"])
    size_y = float(cfg["size_y"])
    size_z = float(cfg["size_z"])

    cpu = int(cfg.get("cpu", 8))
    top_ligands = int(cfg.get("top_ligands", 5))

    # ---------------- Create run directory ---------------- #

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    runs_dir = os.path.join(base_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    run_dir = os.path.join(runs_dir, f"docking_run_{timestamp}")
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    print(f"📁 Working directory: {run_dir}")

    # Run-local outputs
    pdbqt_dir = os.path.join(run_dir, os.path.basename(os.path.normpath(pdbqt_subdir)))
    results_dir = os.path.join(run_dir, os.path.basename(os.path.normpath(results_subdir)))
    os.makedirs(pdbqt_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # Copy YAML into run directory
    shutil.copy(config_file, os.path.join(run_dir, "config_used.yml"))
    print(f"🧾 Copied configuration to: {run_dir}/config_used.yml")

    # ---------------- Helper script paths ---------------- #

    add_ligands_script = os.path.join(base_dir, "add_ligands.py")
    convert_script = os.path.join(base_dir, "convert_to_pdbqt.py")
    vina_script = os.path.join(base_dir, "run_vina_batch.py")
    summary_script = os.path.join(base_dir, "summarize_vina_scores.py")
    prepare_top_script = os.path.join(base_dir, "prepare_top_ligands_for_chimerax.py")

    # ---------------- Sanity checks ---------------- #

    if not os.path.isdir(ligand_dir):
        sys.exit(f"❌ Ligand directory not found: {ligand_dir}")

    ligand_files_list = [
        os.path.join(ligand_dir, f)
        for f in os.listdir(ligand_dir)
        if f.lower().endswith((".sdf", ".mol", ".mol2"))
    ]
    ligand_files_list.sort()

    if not ligand_files_list:
        sys.exit(f"❌ No ligand files found in: {ligand_dir}")

    # ---------------- Step 0: receptor ---------------- #

    receptor_cfg = resolve_path(config_dir, receptor_cfg) if receptor_cfg else ""

    if receptor_cfg and os.path.exists(receptor_cfg):
        print(f"🧬 Using pre-prepared receptor: {receptor_cfg}")

        # Copy into run dir so the run stays self-contained
        receptor_name = os.path.basename(receptor_cfg)
        receptor = os.path.join(run_dir, receptor_name)
        if not os.path.exists(receptor):
            shutil.copy(receptor_cfg, receptor)
    else:
        if not pdb_id:
            sys.exit("❌ No receptor provided: set either 'pdb_id' or 'receptor' in docking_config.yml")
        receptor = prepare_receptor(pdb_id, chains, keep_cofactors)

    # ---------------- Step 1: combine ligands ---------------- #

    ligand_files = " ".join(ligand_files_list)
    cmd1 = f"python3 {add_ligands_script} {combined_sdf} {ligand_files}"
    run_cmd(cmd1, "Combining ligand files")

    # ---------------- Step 2: convert ligands to PDBQT ---------------- #

    cmd2 = f"python3 {convert_script} {combined_sdf} {pdbqt_dir}"
    run_cmd(cmd2, "Converting ligands to PDBQT")

    # ---------------- Step 3: write vina_config.txt ---------------- #

    print(f"\n🚀 Generating {config_out}...")
    write_vina_config(
        config_out,
        receptor,
        center_x, center_y, center_z,
        size_x, size_y, size_z
    )
    print(f"✅ Done: Generating {config_out}")

    # ---------------- Step 4: recenter ligands ---------------- #

    print("\n🧭 Recentering ligand coordinates to docking box center...")

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

        dx = center_x - cx0
        dy = center_y - cy0
        dz = center_z - cz0

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

    # ---------------- Step 5: run docking ---------------- #

    cmd4 = f"python3 {vina_script} {receptor} {pdbqt_dir} {results_dir} {config_out} --cpu {cpu}"
    run_cmd(cmd4, f"Running AutoDock Vina batch docking with {cpu} CPUs")

    # ---------------- Step 6: summarize ---------------- #

    cmd5 = f"python3 {summary_script} {results_dir}"
    run_cmd(cmd5, "Summarizing docking results")

    # ---------------- Step 7: prepare top ligands ---------------- #

    cmd6 = f"python3 {prepare_top_script} {receptor} {results_dir} {top_ligands}"
    run_cmd(cmd6, f"Preparing top {top_ligands} ligands for ChimeraX")

    print("\n🎉 Pipeline completed successfully!")
    print(f"Results stored in: {results_dir}\n")


if __name__ == "__main__":
    main()
