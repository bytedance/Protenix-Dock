#!/usr/bin/env python
"""Fetch PDB structures, extract protein + ligand, compute pocket box.

Usage:
    python scripts/fetch_examples.py
"""
import gzip
import io
import os
import shutil
import subprocess
import tempfile
import urllib.request

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.join(SCRIPT_DIR, "..", "examples")

SYSTEMS = [
    ("8AY3", "OE3"),
    ("7BMI", "U4B"),
    ("5SB2", "1K2"),
    ("7T0D", "FPP"),
]

BUFFER = 4.0  # Angstrom padding around ligand for box


def fetch_pdb(pdb_id):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz"
    print(f"  Downloading {url}")
    req = urllib.request.urlopen(url)
    return gzip.decompress(req.read()).decode("utf-8")


def extract_protein_pdb(pdb_text, pdb_id):
    """Extract ATOM records, add hydrogens with reduce, write clean PDB."""
    lines = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            lines.append(line)
    lines.append("END")

    # Write raw ATOM-only PDB to temp file
    tmp_dir = tempfile.mkdtemp()
    raw_pdb = os.path.join(tmp_dir, "raw.pdb")
    with open(raw_pdb, "w") as f:
        f.write("\n".join(lines) + "\n")

    # Add hydrogens with reduce (needed for HIS protonation detection)
    out_path = os.path.join(EXAMPLES_DIR, f"{pdb_id.lower()}_protein_prepared.pdb")
    try:
        result = subprocess.run(
            ["reduce", "-BUILD", "-Quiet", raw_pdb],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0 and not result.stdout:
            raise RuntimeError(f"reduce failed: {result.stderr[:200]}")
        # Filter to only ATOM + END lines (reduce adds USER MOD headers)
        out_lines = []
        for line in result.stdout.splitlines():
            if line.startswith("ATOM") or line.startswith("END"):
                out_lines.append(line)
        with open(out_path, "w") as f:
            f.write("\n".join(out_lines) + "\n")
        n_heavy = sum(1 for l in lines if l.startswith("ATOM"))
        n_total = sum(1 for l in out_lines if l.startswith("ATOM"))
        print(f"  Wrote protein: {out_path} ({n_heavy} heavy + {n_total - n_heavy} H atoms)")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    return out_path


def extract_ligand_sdf(pdb_text, pdb_id, lig_id):
    """Extract ligand HETATM records, convert to SDF via RDKit."""
    # Build a minimal PDB with just the ligand
    lig_lines = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname == lig_id:
                lig_lines.append(line)
    if not lig_lines:
        raise ValueError(f"Ligand {lig_id} not found in {pdb_id}")

    lig_lines.append("END")
    lig_pdb = "\n".join(lig_lines) + "\n"

    mol = Chem.MolFromPDBBlock(lig_pdb, removeHs=True, sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit failed to parse ligand {lig_id} from {pdb_id}")

    # Keep only first conformer / first occurrence
    out_path = os.path.join(EXAMPLES_DIR, f"{pdb_id.lower()}_ligand_prepared.sdf")
    writer = Chem.SDWriter(out_path)
    writer.write(mol)
    writer.close()
    print(f"  Wrote ligand: {out_path} ({mol.GetNumAtoms()} heavy atoms)")
    return out_path, mol


def compute_pocket(mol, buffer=BUFFER):
    """Compute pocket center and box size from ligand conformer."""
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0) + 2 * buffer
    return center, span


def write_pocket_config(pdb_id, center, size):
    out_path = os.path.join(EXAMPLES_DIR, f"{pdb_id.lower()}_pocket.config")
    with open(out_path, "w") as f:
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n")
    print(f"  Wrote pocket: {out_path}")
    return out_path


def main():
    os.makedirs(EXAMPLES_DIR, exist_ok=True)
    for pdb_id, lig_id in SYSTEMS:
        print(f"\n=== {pdb_id} / {lig_id} ===")
        pdb_text = fetch_pdb(pdb_id)
        extract_protein_pdb(pdb_text, pdb_id)
        _, mol = extract_ligand_sdf(pdb_text, pdb_id, lig_id)
        center, size = compute_pocket(mol)
        write_pocket_config(pdb_id, center, size)
    print("\nDone!")


if __name__ == "__main__":
    main()
