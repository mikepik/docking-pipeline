#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import re

def convert_sdf_to_pdbqt(sdf_path: Path, out_dir: Path, obabel_bin: str = "obabel"):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        obabel_bin,
        "-isdf", str(sdf_path),
        "-O", str(out_dir / "ligand_.pdbqt"),  # base name; "-m" splits automatically
        "-m"  # split multi-ligand SDF into individual files
    ]

    print("🔄 Running Open Babel:")
    print("   " + " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Count the number of .pdbqt files created
    pdbqt_files = list(out_dir.glob("*.pdbqt"))
    num_files = len(pdbqt_files)

    # Optionally check that at least one file exists
    if num_files == 0:
        print(f"⚠️  No PDBQT files found in {out_dir}. Check that your SDF file contained valid molecules.")
    else:
        # Try to get the number of ligands in the SDF file (for reference)
        try:
            sdf_text = sdf_path.read_text(errors="ignore")
            num_in_sdf = len(re.findall(r"\$\$\$\$", sdf_text))
            print(f"✅ Converted {num_files} ligands (of {num_in_sdf} entries in SDF) successfully.")
        except Exception:
            print(f"✅ Converted {num_files} ligands successfully.")

    print(f"📂 Output directory: {out_dir.resolve()}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert multi-ligand SDF into individual PDBQT files using Open Babel.")
    p.add_argument("--sdf", required=True, help="Input multi-ligand .sdf file")
    p.add_argument("--out_dir", required=True, help="Output directory for .pdbqt ligands")
    p.add_argument("--obabel_bin", default="obabel", help="Path to obabel executable (default: obabel in PATH)")
    args = p.parse_args()

    convert_sdf_to_pdbqt(Path(args.sdf).resolve(), Path(args.out_dir).resolve(), args.obabel_bin)
