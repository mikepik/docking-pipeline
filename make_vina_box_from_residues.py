#!/usr/bin/env python3
import os, sys, argparse
from statistics import mean

def parse_residue_list(reslist):
    residues = []
    for item in reslist.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            chain, resid = item.split(":")
            residues.append((chain.strip(), int(resid)))
        except ValueError:
            sys.exit(f"❌ Invalid residue specifier: {item}. Use format Chain:Resid")
    return residues

def get_residue_coords(pdbfile, residues):
    xs, ys, zs = [], [], []
    targets = set(residues)
    with open(pdbfile, "r") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip()
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            if (chain, resid) in targets:
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    pass
    if not xs:
        sys.exit("❌ No matching residues found in receptor file.")
    return xs, ys, zs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor", help="Receptor PDBQT file")
    ap.add_argument("residues", help="Comma-separated residues (Chain:Resid)")
    ap.add_argument("--padding", type=float, default=5.0, help="Padding Å")
    ap.add_argument("--out", default="vina_config.txt", help="Output config name")
    args = ap.parse_args()

    if not os.path.exists(args.receptor):
        sys.exit(f"❌ Receptor file not found: {args.receptor}")

    residues = parse_residue_list(args.residues)
    xs, ys, zs = get_residue_coords(args.receptor, residues)

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)

    cx, cy, cz = mean([xmin, xmax]), mean([ymin, ymax]), mean([zmin, zmax])
    sx, sy, sz = (xmax - xmin) + 2*args.padding, (ymax - ymin) + 2*args.padding, (zmax - zmin) + 2*args.padding

    with open(args.out, "w") as f:
        f.write(f"receptor = {args.receptor}\n")
        f.write(f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n")
        f.write(f"size_x = {sx:.3f}\nsize_y = {sy:.3f}\nsize_z = {sz:.3f}\n")
        f.write("exhaustiveness = 8\nnum_modes = 9\nenergy_range = 3\n")

    print("📦 Vina box parameters:")
    print(f"  center_x = {cx:.3f}")
    print(f"  center_y = {cy:.3f}")
    print(f"  center_z = {cz:.3f}")
    print(f"  size_x   = {sx:.3f}")
    print(f"  size_y   = {sy:.3f}")
    print(f"  size_z   = {sz:.3f}")
    print(f"✅ Config written to: {args.out}")

if __name__ == "__main__":
    main()
