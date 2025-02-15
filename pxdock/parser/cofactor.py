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

import MDAnalysis as mda
from rdkit import Chem

from pxdock.common import get_logger, kWorkDir, my_random_string
from pxdock.common.amber_forcefield import ION_RESNAME_TO_ATOMNAME
from pxdock.common.constants import ATOMIC_NUMBERS_TO_ELEMENTS
from pxdock.common.utilities import convert_value_to_list
from pxdock.parser.ffdata import parm_to_ffdata, pdb_to_parm
from pxdock.parser.ligand import LigandParser

logger = get_logger(__name__)


class CofactorParser(LigandParser):
    def __init__(self, *args, ionforcefield: str = "FF14SB", **kwargs):
        self.ionforcefield = ionforcefield
        _tmp_working_dir = os.path.join(kWorkDir, my_random_string(10))
        self._working_dir = _tmp_working_dir
        os.makedirs(_tmp_working_dir, exist_ok=True)
        super().__init__(*args, **kwargs)

    def pdb_to_ffdata(self, ligand):
        working_dir = os.path.abspath(self._working_dir)
        pdbfile = os.path.abspath(ligand)
        # step 0: rename metal ions
        u = mda.Universe(pdbfile)
        if u.residues[0].resname in ION_RESNAME_TO_ATOMNAME:
            assert len(u.residues) == 1
            assert len(u.atoms) == 1
            u.atoms[0].name = ION_RESNAME_TO_ATOMNAME[u.residues[0].resname]
            step0_pdb = os.path.join(working_dir, "step0_cleaned.pdb")
            with mda.Writer(step0_pdb) as w:
                w.write(u)
        elif u.residues[0].resname == "HEM":
            step0_pdb = pdbfile
        else:
            raise ValueError(".pdb format only supports ions and hem now.")

        # step 1, 2 convert to parmed object
        parm = pdb_to_parm(working_dir, step0_pdb, self.ionforcefield)

        # step 3: convert to mol2 format
        mol2 = os.path.join(working_dir, "step3.mol2")
        for atom in parm.atoms:
            atom.atom_type.atomic_number = atom.atomic_number
            atom.type = ATOMIC_NUMBERS_TO_ELEMENTS[atom.atomic_number]
        parm.save(mol2, overwrite=True)
        # step 4: parse ffdata from parmed object
        logger.info("parsing ffdata and xyz from parmed object")
        ffdata, xyz = parm_to_ffdata(parm)
        # transform to list
        ffdata = convert_value_to_list(ffdata)
        logger.info("parse ffdata from parmed object finished")
        rkmol = Chem.MolFromMol2File(
            mol2, removeHs=False, sanitize=False, cleanupSubstructures=False
        )
        molecule = self._from_iterable_rkmol([rkmol], one_conformer_each_rkmol=True)
        return molecule, ffdata, xyz

    @staticmethod
    def rdmol_to_gro(rkmol, working_dir):
        """
        Convert an RDKit molecule object to a GROMACS GRO file format string and save the GRO file.

        Args:
            rkmol (Chem.rdchem.Mol): rdkit mol object
            working_dir (String): path to save the gro file

        Returns:
            string: gro string of rkmol
        """
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        cofacator_gro_path = os.path.join(working_dir, "cofactor.gro")
        u = mda.Universe(rkmol)
        with mda.Writer(cofacator_gro_path, reindex=False, format="GRO") as gro_write:
            gro_write.write(u.atoms)
        with open(cofacator_gro_path, "r", encoding="utf-8") as gro_file:
            gro_string = gro_file.read()
        return gro_string

    def get_data(self):
        cofactor_data = {
            "gmx_gro_string": self.rdmol_to_gro(self.rkmol, self._working_dir),
            "ffdata": self.ffdata,
            "xyz": self.xyz,
            "bond_index": self.bond_index,
            "is_rotatable": self.is_rotatable,
        }
        return convert_value_to_list(cofactor_data)


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fpath", type=str, required=True)
    parser.add_argument("--output_fpath", type=str, required=True)
    parser.add_argument("--parm_save_path", type=str, default=None)
    args = parser.parse_args()
    cofactor = CofactorParser(
        ligand=args.input_fpath,
        parm_save_path=args.parm_save_path,
    )
    json.dump(
        cofactor.get_data(),
        open(args.output_fpath, "w", encoding="utf-8"),
        indent=4,
    )
