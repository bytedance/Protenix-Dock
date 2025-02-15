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

import copy
import time
import traceback
from collections.abc import Iterable
from typing import List, Union

import numpy as np
from byteff.mol import Molecule, rkutil
from byteff.mol.conformer import Conformer, prefix_prop_key
from func_timeout import exceptions as func_timeout_exceptions
from func_timeout import func_set_timeout
from rdkit import Chem

from pxdock.common import get_logger
from pxdock.common.pdbqt_utils import MyPDBQTMolecule
from pxdock.common.utilities import (
    convert_value_to_list,
    run_with_timeout,
    save_sdf_supplier,
)
from pxdock.dock.utils import kScoreConfigFile
from pxdock.parser.base import MolBase, TorsionStrainParams
from pxdock.parser.ligand_ffdata import calculate_molecular_weight, get_ffdata
from pxdock.parser.ligand_rbonds import get_rdkit_rotatable_bonds

logger = get_logger(__name__)


@func_set_timeout(5)
def _export_pdbqt_to_rkmol(pmol):
    return pmol.export_rdkit_mol(exclude_hydrogens=False, return_index_map=False)


def export_pdbqt_to_rkmol_with_timeout(pdbqt_mol):
    try:
        rkmol = _export_pdbqt_to_rkmol(pdbqt_mol)
    except func_timeout_exceptions.FunctionTimedOut as exc:
        tb = traceback.format_exc()
        logger.warning(f"[ERROR export_pdbqt_to_rkmol]: {str(exc)} \n {tb}")
        return None
    return rkmol


def read_sdf_with_timeout(sdf_file, timeout_duration=5):

    def _read_sdf(x):
        mol_list = []
        for mol in save_sdf_supplier(x, skip_failed_mol=False):
            mol_list.append(mol)
        return mol_list

    return run_with_timeout(
        _read_sdf, args=(sdf_file,), timeout_duration=timeout_duration
    )


def to_iterable_rkmol(ligand: str):
    if ligand.endswith(".sdf"):
        return read_sdf_with_timeout(ligand)
    elif ligand.endswith(".mol2"):
        return [
            Chem.MolFromMol2File(
                ligand, removeHs=False, sanitize=False, cleanupSubstructures=False
            )
        ]
    else:
        if ligand.endswith(".pdbqt"):
            with open(ligand, "r", encoding="utf-8") as input_file:
                pdbqt_string = input_file.read()
            ligand = pdbqt_string
        # get rkmol iterator
        pdbqt_mol = MyPDBQTMolecule(ligand)
        rkmol_list = []
        for i, t_mol in enumerate(pdbqt_mol):
            rkmol = export_pdbqt_to_rkmol_with_timeout(t_mol)
            rkmol_list.append(rkmol)
        logger.info(f"Got a iteratble rkmol list with {len(rkmol_list)} rkmol.")
        return rkmol_list


def generate_candidate_confs(mapped_smiles, num_candidate_confs, ffopt=True):
    logger.info(f"generate {num_candidate_confs} candidate conformers")
    mol = Molecule.from_mapped_smiles(
        mapped_smiles, nconfs=num_candidate_confs, ffopt=ffopt
    )
    return [conf.coords.tolist() for conf in mol.conformers]


class LigandParser:

    def __init__(
        self,
        ligand: Union[
            str, Chem.Mol, Molecule
        ],  # pdbqt_string/sdf_filename/rdkit mol/bytemol
        forcefield: str = "gaff2",
        charge_model: str = "am1bcc",
        charge_conformer: str = "rdkit",
        parm_save_path: str = None,
        **kwargs,
    ):
        """
        Args:
            ligand:
                - str:
                    - sdf filename:
                        - The file name should end with ".sdf".
                        - All molecules in the SDF file will be read. It is required that all
                            molecules in the SDF file are different poses of **the same molecule**.
                    - pdbqt filename:
                        - The file name should end with ".pdbqt".
                        - The PDBQT string in the file will be read.
                    - pdbqt_string:
                        - A PDBQT string containing multiple poses of **the same molecule**.
                    - pdb filename:
                        - PDB format files are no longer supported as ligand inputs.
                        - PDB files are only used as inputs for cofactors such as HEME and IONs.
                - Chem.Mol:
                    - A Chem.Mol object that can contain multiple poses. [It must already contain the coordinates of hydrogen atoms.]
                - byteff mol Molecule:
                    - A Molecule object that can contain multiple poses.
            forcefield: Force field parameters.
                - Currently available: "gaff2".
            parm_save_path:
                - if None, the ligand parm file will be removed.
                - Otherwise, it will be saved to the specified path.
            charge_model: am1bcc or abcg2 or path to bcc json file
            charge_conformer: "rdkit" or "vina", the conformer used for computing partial charge
        """

        ### initilize molecule ###
        self._rkmol = None
        self._mapped_smiles = None
        self.ffdata = None
        xyz = None

        self._molecule, self.ffdata = self.process_input(ligand)
        if self.molecule.nconfs < 1:
            raise ValueError("No conformer is read.")
        else:
            logger.info(f"{self.molecule.nconfs} conformers are read in LigandParser.")

        ### ffdata ###
        if self.ffdata is not None:
            pass
        elif charge_conformer == "rdkit":
            molecule = copy.deepcopy(self.molecule)
            # remove confs
            molecule._conformers = []
            # generate a rdkit conformer
            rkmol = molecule.rkmol
            rkmol.RemoveAllConformers()
            try:
                rkmol, _, _ = rkutil.generate_confs(rkmol, nconfs=1)
            except Exception as e:
                Chem.SetUseLegacyStereoPerception(True)
                rkmol, _, _ = rkutil.generate_confs(rkmol, nconfs=1)
            finally:
                Chem.SetUseLegacyStereoPerception(False)
            # append to molecule
            conformers = self.get_conformers_from_rkmol(rkmol)
            molecule.append_conformers(conformers)
            logger.info(
                f"{molecule.nconfs} rdkit conformer is generated for computing partial charge."
            )
            self.ffdata = get_ffdata(
                molecule,
                forcefield=forcefield,
                charge_model=charge_model,
                parm_save_path=parm_save_path,
                conf_id=0,
            )
        else:
            self.ffdata = get_ffdata(
                self._molecule,
                forcefield=forcefield,
                charge_model=charge_model,
                parm_save_path=parm_save_path,
                conf_id=0,
            )

        ### xyz ###
        self.xyz = self.get_coords_from_molecule(self.molecule)
        if xyz is not None:
            assert np.all(np.isclose(self.xyz, np.array([xyz]), atol=1e-2))
        self.num_poses = self.molecule.nconfs
        assert self.num_poses == len(self.xyz)

        ### bond_index ###
        self.bonds = self.rkmol.GetBonds()
        self.num_bonds = len(self.bonds)
        self.bond_index = [
            sorted([int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())])
            for bond in self.bonds
        ]

        ### is_rotatable ###
        rotatable_bonds = get_rdkit_rotatable_bonds(self.rkmol)
        self.rotatable_bonds = [sorted(list(bond)) for bond in rotatable_bonds]
        self.is_rotatable = [
            1 if set(bond) in rotatable_bonds else 0 for bond in self.bond_index
        ]

        self.score_function_json = kwargs.pop("score_function_json", None)
        if self.score_function_json is None:
            logger.info("score_function_json not provided, using bytedock default")
            self.score_function_json = kScoreConfigFile
        logger.info(
            f"user should ensure consistent config in parser and score function, score_function_json: {self.score_function_json}"
        )
        # Find specific atom types or chemical environments
        self.specific_atom_finder = MolBase(
            self.rkmol, score_function_json=self.score_function_json, **kwargs
        )
        self.hydroph_atoms = self.specific_atom_finder.get_hydrophobic_atoms
        # new add for hydrophic groups
        self.hbond_acc_atoms = self.specific_atom_finder.get_hba
        self.hbond_don_atoms = self.specific_atom_finder.get_hbd
        self.halogen_acc_atoms = self.specific_atom_finder.get_xba
        self.halogen_don_atoms = self.specific_atom_finder.get_xbd
        self.cation_atoms = self.specific_atom_finder.get_cationa
        self.anion_atoms = self.specific_atom_finder.get_aniona
        self.pi_ring_atoms = self.specific_atom_finder.get_pi_ring
        self.metal_donar_atoms = self.specific_atom_finder.get_metal_donar
        self.metal_acceptor_atoms = self.specific_atom_finder.get_metal_acceptor
        self.torsion_strain_bonds = self.specific_atom_finder.get_torsion_strain_bonds

        # ffdata: add additional information
        self.update_ffdata()
        # formal_charge = [atom.GetFormalCharge() for atom in self.rkmol.GetAtoms()]
        self.ffdata.update(
            {
                "atomic_numbers": [
                    atom.GetAtomicNum() for atom in self.rkmol.GetAtoms()
                ],
                "mapped_smiles": self.mapped_smiles,
                "bond_index": self.bond_index,
                "is_rotatable": self.is_rotatable,
                "rotatable_bond_index": self.rotatable_bonds,
            }
        )

        # remove keys with "_paraidx" except for FF_vdw_paraidx and tstrain_paraidx
        keys = list(self.ffdata.keys())
        for key in keys:
            if key.endswith("_paraidx") and key not in [
                "FF_vdW_paraidx",
                "tstrain_paraidx",
            ]:
                self.ffdata.pop(key)

        logger.info("Initialize Done!")

    def process_input(self, ligand):
        if isinstance(ligand, Molecule):
            return ligand, None
        elif isinstance(ligand, Chem.Mol):
            if len(ligand.GetConformers()) < 1:
                raise ValueError("The input Chem.Mol does not contain conformers.")
            return Molecule(ligand, keep_conformers=True), None
        elif isinstance(ligand, str):
            if ligand.endswith(".pdb"):
                molecule, ffdata, _ = self.pdb_to_ffdata(ligand)
                return molecule, ffdata
            st = time.time()
            rkmols = to_iterable_rkmol(ligand)
            mt = time.time()
            tt = round(mt - st, 3)
            logger.info(f"Convert ligands to rkmol in {tt} seconds.")
            molecule = self._from_iterable_rkmol(rkmols, one_conformer_each_rkmol=True)
            tt = round(time.time() - mt, 3)
            logger.info(f"Convert rkmol to byteff molecule in {tt} seconds.")
            return molecule, None

    def pdb_to_ffdata(self, ligand):
        raise ValueError(f"PDB file {ligand} is not supported.")

    @staticmethod
    def _from_iterable_rkmol(rkmol_iterator, one_conformer_each_rkmol=False):
        # if one_conformer_each_rkmol is True:
        # each rkmol in the iterator must contain exactly one conformer, otherwise an error will be raised.
        assert isinstance(rkmol_iterator, Iterable), "rkmol_iterator is not iterable."
        molecule = None
        all_conformers = []
        for index, rkmol in enumerate(rkmol_iterator):
            if rkmol is None:
                continue
            if molecule is None:
                # initialize the byteff Molecule
                molecule = Molecule(rkmol, keep_conformers=True)
                if one_conformer_each_rkmol:
                    assert molecule.nconfs == 1  # contains the pose in rkmol_wH
            else:
                # append the conformer of rkmol_wH to molecule
                conformers = LigandParser.get_conformers_from_rkmol(rkmol)
                if one_conformer_each_rkmol:
                    assert len(conformers) == 1
                all_conformers.extend(conformers)
        if molecule is None:
            raise ValueError("[LigandParser] No rkmol is read successfully.")

        # this append function will check whether the order of atomic symbols match.
        molecule.append_conformers(all_conformers)

        return molecule

    @property
    def molecule(self):
        return self._molecule

    @property
    def mapped_smiles(self):
        if self._mapped_smiles is None:
            self._mapped_smiles = self.molecule.get_mapped_smiles(isomeric=True)
        return self._mapped_smiles

    @property
    def rkmol(self):
        if self._rkmol is None:
            self._rkmol = self.molecule.to_rkmol()
        return self._rkmol

    @staticmethod
    def get_coords_from_molecule(molecule: Molecule):
        coords = []
        for conf in molecule.conformers:
            coords.append(conf.coords)
        return np.stack(coords, axis=0)

    @staticmethod
    def get_conformers_from_rkmol(rkmol: Chem.Mol):
        """
        parse conformers rkmol
        rkmol should have all atoms EXPLICIT
        """
        rkconfs = rkmol.GetConformers()
        if len(rkconfs) == 0:
            return []

        natoms = rkmol.GetNumAtoms(onlyExplicit=False)
        atomic_numbers = [at.GetAtomicNum() for at in rkmol.GetAtoms()]

        all_conformers = []

        for conf in rkconfs:
            coords = np.vstack([conf.GetAtomPosition(idx) for idx in range(natoms)])
            conformer = Conformer(coords, atomic_numbers)
            confprops = conf.GetPropsAsDict()
            for k, v in confprops.items():
                prop = prefix_prop_key(k)
                conformer.confdata[prop] = v
            all_conformers.append(conformer)
        return all_conformers

    @staticmethod
    def get_torsion_strain_ffdata(
        tstrain_bonds: List, score_function_json: str = None, task="pose_selection"
    ):
        tstrain_ffdata = dict()
        tstrain_paramslist = TorsionStrainParams(score_function_json, task=task)
        tstrain_paraidx = [x.paraidx for x in tstrain_bonds]
        tstrain_ffdata[f"tstrain_paraidx_{task}"] = tstrain_paraidx
        tstrain_ffdata[f"tstrain_params_{task}"] = [
            tstrain_paramslist.get_params_by_idx(idx) for idx in tstrain_paraidx
        ]
        if tstrain_bonds:
            # [ntorsion, nmatch, 4]
            tstrain_atomidx = [x.proper_atomidx for x in tstrain_bonds]
        else:
            tstrain_atomidx = []
        tstrain_ffdata[f"tstrain_atomidx_{task}"] = tstrain_atomidx
        return tstrain_ffdata

    def update_ffdata(self):
        self.ffdata["hydrophobic_atomidx"] = [x.idx for x in self.hydroph_atoms]
        self.ffdata["hbondacc_atomidx"] = [x.idx for x in self.hbond_acc_atoms]
        self.ffdata["hbonddon_atomidx"] = [x.idx for x in self.hbond_don_atoms]
        self.ffdata["halogenacc_atomidx"] = [x.idx for x in self.halogen_acc_atoms]
        self.ffdata["halogendon_atomidx"] = [x.idx for x in self.halogen_don_atoms]
        self.ffdata["cation_atomidx"] = [x.idx for x in self.cation_atoms]
        self.ffdata["anion_atomidx"] = [x.idx for x in self.anion_atoms]
        self.ffdata["piring5_atomidx"] = [
            x.idx for x in self.pi_ring_atoms if x.type == "five_mem_ring"
        ]
        self.ffdata["piring6_atomidx"] = [
            x.idx for x in self.pi_ring_atoms if x.type == "six_mem_ring"
        ]
        self.ffdata["metaldonar_atomidx"] = [x.idx for x in self.metal_donar_atoms]
        self.ffdata["metalaccep_atomidx"] = [x.idx for x in self.metal_acceptor_atoms]

        # H-bond根据是否带电细分
        self.ffdata["hbonddon_charged_atomidx"] = [
            x
            for x in self.ffdata["hbonddon_atomidx"]
            if x[0] in (self.ffdata["cation_atomidx"] + self.ffdata["anion_atomidx"])
        ]
        self.ffdata["hbondacc_charged_atomidx"] = [
            x
            for x in self.ffdata["hbondacc_atomidx"]
            if x in (self.ffdata["cation_atomidx"] + self.ffdata["anion_atomidx"])
        ]
        self.ffdata["hbonddon_neut_atomidx"] = [
            x
            for x in self.ffdata["hbonddon_atomidx"]
            if x[0]
            not in (self.ffdata["cation_atomidx"] + self.ffdata["anion_atomidx"])
        ]
        self.ffdata["hbondacc_neut_atomidx"] = [
            x
            for x in self.ffdata["hbondacc_atomidx"]
            if x not in (self.ffdata["cation_atomidx"] + self.ffdata["anion_atomidx"])
        ]

        for task, bonds in self.torsion_strain_bonds.items():
            tstrain_ffdata = self.get_torsion_strain_ffdata(
                bonds, score_function_json=self.score_function_json, task=task
            )
            self.ffdata.update(tstrain_ffdata)

        self.ffdata["molecular_weight"] = self.get_molecular_weight()

    def get_data(self):
        ligand_data = {
            "mapped_smiles": self.mapped_smiles,
            "ffdata": self.ffdata,
            "xyz": self.xyz,
            "bond_index": self.bond_index,
            "is_rotatable": self.is_rotatable,
        }
        return convert_value_to_list(ligand_data)

    def get_molecular_weight(self):
        """update molecular weight"""
        return calculate_molecular_weight(self.mapped_smiles)
