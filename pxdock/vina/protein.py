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
import subprocess

import numpy as np
import rdkit.Chem as Chem

try:
    import AutoDockTools as ADT

    ADTPath = ADT.__path__[0]
except ImportError:
    raise f"import AutoDockTools error, please install `AutoDockTools`"


class PrepProt(object):
    def __init__(self, pdb_file: str):
        self.prot = pdb_file
        self.prot_pqr = None

    def del_water(self, output: str):  # optional
        if not os.path.exists(output):
            with open(self.prot) as f:
                lines = [
                    l
                    for l in f.readlines()
                    if l.startswith("ATOM") or l.startswith("HETATM")
                ]
                dry_lines = [l for l in lines if not "HOH" in l]
            with open(output, "w") as f:
                f.write("".join(dry_lines))
        self.prot = output

    def del_metal(self, output: str):
        metals = set(
            "ZN",
            "CA",
            "HG",
            "MG",
            "AU",
            "MN",
            "NA",
            "CU",
            "CO",
            "NI",
            "FE",
            "K",
            "GA",
            "SR",
            "BA",
            "CS",
            "LI",
            "RB",
            "IN",
            "SE",
            "CD",
        )

        if not os.path.exists(output):
            with open(self.prot) as f:
                lines = [
                    l
                    for l in f.readlines()
                    if l.startswith("ATOM") or l.startswith("HETATM")
                ]

            non_metal_lines = [l for l in lines if l[17:20].strip() not in metals]
            with open(output, "w") as f:
                f.writelines(non_metal_lines)

        self.prot = output

    def addH(self, prot_pqr: str, prot_pdb: str):  # call pdb2pqr
        self.prot_pqr = prot_pqr
        self.prot_pdb = prot_pdb
        if not os.path.exists(prot_pqr):
            subprocess.Popen(
                [
                    "pdb2pqr30",
                    "--ff=AMBER",
                    self.prot,
                    prot_pqr,
                    "--pdb-output",
                    prot_pdb,
                ]
            ).communicate()

    def get_pdbqt(self, prot_pdbqt: str):
        if self.prot_pqr is None:
            self.prot_pqr = self.prot
        if not os.path.exists(prot_pdbqt):
            prepare_receptor = f"{ADTPath}/Utilities24/prepare_receptor4.py"
            subprocess.Popen(
                ["python3", prepare_receptor, "-r", self.prot_pdb, "-o", prot_pdbqt]
            ).communicate()

    # def get_box(self, ref_path, buf=0):
    #     self.ref_path = ref_path
    #     with open(ref_path, 'r') as f:
    #         lines = [l for l in f.readlines() if l.startswith('ATOM') or l.startswith('HETATM')]

    #         xs = [float(l[30:38]) for l in lines]
    #         ys = [float(l[38:46]) for l in lines]
    #         zs = [float(l[46:54]) for l in lines]
    #         self.pocket_center = [round((max(xs) + min(xs)) / 2, 3), round((max(ys) + min(ys)) / 2, 3),
    #                               round((max(zs) + min(zs)) / 2, 3)]
    #         self.box_size = [round((max(xs) - min(xs)) + buf, 3), round((max(ys) - min(ys)) + buf, 3),
    #                          round((max(zs) - min(zs)) + buf, 3)]

    def get_box(self, ref_path: str, buf=0):
        # Find the pocket according to the ref ligand
        if ref_path.endswith(".pdb") or ref_path.endswith(".pdbqt"):
            with open(ref_path, "r") as f:
                lines = [
                    l
                    for l in f.readlines()
                    if l.startswith("ATOM") or l.startswith("HETATM")
                ]

            xs = [float(l[30:38]) for l in lines]
            ys = [float(l[38:46]) for l in lines]
            zs = [float(l[46:54]) for l in lines]

            coords = np.stack((xs, ys, zs), axis=-1)

        elif (
            ref_path.endswith(".sdf")
            or ref_path.endswith(".mol")
            or ref_path.endswith(".mol2")
        ):
            if ref_path.endswith(".mol2"):
                mol = Chem.MolFromMol2File(ref_path, removeHs=False)
            else:
                mol = Chem.MolFromMolFile(ref_path, removeHs=False)
            coords = mol.GetConformers()[0].GetPositions()

        self.pocket_center = tuple(
            np.round((np.max(coords, axis=0) + np.min(coords, axis=0)) / 2, 3).tolist()
        )
        self.box_size = tuple(
            np.round((np.max(coords, axis=0) - np.min(coords, axis=0)) + buf, 3).tolist()
        )
