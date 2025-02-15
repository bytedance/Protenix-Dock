# Copyright (C) 2025 ByteDance and/or its affiliates

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import re

import numpy as np
from rdkit import Chem


def read_coords_from_sdf_file(sdf_path):
    file = open(sdf_path, mode="r")
    content = file.read()
    file.close()
    match = re.search(r"^ {2,}-?\d.*(?:\r?\n {3,}-?\d.*)*", content, re.M)
    datasplit = []
    if match:
        for line in match.group().splitlines():
            datasplit.append([float(part) for part in line.split()[:3]])

    return np.array(datasplit)


def compute_pocket_box(file_name, buffer=5):
    # Find the pocket according to the ref ligand
    if file_name.endswith(".pdb") or file_name.endswith(".pdbqt"):
        with open(file_name, "r") as f:
            lines = [
                l
                for l in f.readlines()
                if l.startswith("ATOM") or l.startswith("HETATM")
            ]
            xs = [float(l[30:38]) for l in lines]
            ys = [float(l[38:46]) for l in lines]
            zs = [float(l[46:54]) for l in lines]
            pocket_center = [
                round((max(xs) + min(xs)) / 2, 3),
                round((max(ys) + min(ys)) / 2, 3),
                round((max(zs) + min(zs)) / 2, 3),
            ]
            box_size = [
                round((max(xs) - min(xs)) + buffer, 3),
                round((max(ys) - min(ys)) + buffer, 3),
                round((max(zs) - min(zs)) + buffer, 3),
            ]
    elif (
        file_name.endswith(".sdf")
        or file_name.endswith(".mol")
        or file_name.endswith(".mol2")
    ):
        if file_name.endswith(".mol2"):
            mol = Chem.MolFromMol2File(file_name, removeHs=False)
            coords = mol.GetConformers()[0].GetPositions()
            # TODO: if failed, read from mol2 file
        else:
            try:
                mol = Chem.MolFromMolFile(file_name, removeHs=False)
                coords = mol.GetConformers()[0].GetPositions()
            except:
                coords = read_coords_from_sdf_file(file_name)
        pocket_center = list(
            np.round((np.max(coords, axis=0) + np.min(coords, axis=0)) / 2, 3)
        )
        pocket_center = [item.item() for item in pocket_center]

        box_size = list(
            np.round((np.max(coords, axis=0) - np.min(coords, axis=0)) + buffer, 3)
        )
        box_size = [item.item() for item in box_size]

    return pocket_center, box_size


def list_pdbbind_dir(path):
    file_lists = os.listdir(path)
    file_lists = [x for x in file_lists if len(x) == 4 and x[0].isdigit()]
    return file_lists


def write_pocket_conf(file, center, pocket):
    with open(file, "w") as f:
        f.writelines("center:" + ",".join(center) + "\n")
        f.writelines("box:" + ",".join(pocket) + "\n")


def read_pocket_conf(file):
    with open(file, "r") as f:
        lines = f.readlines()
        center = list(map(float, lines[0].strip("center:").split(",")))
        box = list(map(float, lines[1].strip("box:").split(",")))
    return center, box


if __name__ == "__main__":
    targets = list_pdbbind_dir("./")
    print(targets)
    for target in targets:
        ref_pdb_path = f"./{target}/{target}_ligand_prepared.pdb"
        pocket_conf = f"./{target}/{target}_pocket.txt"
        if os.path.exists(ref_pdb_path):
            pocket_center, box_size = compute_pocket_box(ref_pdb_path)
            print(pocket_center, box_size)
            pocket_center, box_size = compute_pocket_box(
                ref_pdb_path.replace(".pdb", ".sdf")
            )
            # pocket_center = [str(round(x, 3)) for x in pocket_center]
            # box_size = [str(round(x, 3)) for x in box_size]
            print(pocket_center, box_size)
            # write_pocket_conf(pocket_conf, pocket_center, box_size)
