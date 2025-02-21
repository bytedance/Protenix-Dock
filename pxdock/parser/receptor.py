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

import json
import os
import shutil
from copy import deepcopy
from io import StringIO

import MDAnalysis as mda
import numpy as np

from pxdock.common import get_logger, kWorkDir, my_random_string
from pxdock.common.constants import ELEMENTS_TO_ATOMIC_NUMBERS
from pxdock.common.utilities import convert_value_to_list
from pxdock.geometry.geometry import ReceptorGeometry
from pxdock.parser.ffdata import parm_to_ffdata, pdb_to_parm
from pxdock.parser.molecule_toolbox import generate_sanitized_protein_pdb

logger = get_logger(__name__)

LIST_PARAMS = [
    "rotatable_hbd_bond",
    "aromatic_ring",
]
INDEX_PARAMS = [
    "vina_hydrophobic",
    "vina_donor",
    "vina_acceptor",
    "cation_atoms",
    "anion_atoms",
    "rotatable_h",
]


class ReceptorParser:
    """_summary_"""

    def __init__(
        self,
        pdb_file=None,
        forcefield="FF14SB",
        find_pocket=True,
        find_frozen=True,
        center=None,
        box=None,
        parm_output_dir=None,
        receptor_fpath=None,
    ):
        """
        Args:
            pdb_file (str): Path to the PDB file of the receptor.
            forcefield (str, optional): Force field to be used. Defaults to "FF14SB".
            find_pocket (bool, optional): Whether to consider the hydrogen atoms far away from pocket not to be rotatable. Defaults to True.
            find_frozen (bool, optional): Whether to remove intramolecular interactions of frozen atoms. Defaults to True.
            center (list, optional): Center coordinates of the binding pocket. Defaults to None.
            box (list, optional): Dimensions of the binding pocket. Defaults to None.
            parm_output_dir (str, optional): if not None, the parm files (prmtop, inpcrd, gro etc) will be saved to this directory.
            receptor_fpath (str, optional): if not None, load the receptor data from this path.
        """
        _tmp_working_dir = os.path.join(kWorkDir, my_random_string(10))
        self.working_dir = parm_output_dir or _tmp_working_dir

        if receptor_fpath:
            self.load_data(receptor_fpath)
            return

        self.receptor_path = pdb_file
        self.center = center
        self.box = box
        # generate gro file and ffdata
        gmx_gro, self.ffdata, xyz = self.pdb_to_ffdata(
            self.working_dir, pdb_file, forcefield
        )
        self.xyz = np.array(xyz)
        with open(gmx_gro, "r") as f:
            self.gmx_gro_string = f.read()
        self.bond_index = self.ffdata["FF_Bonds_atomidx"]

        # read residue parameters .json file
        residue_params_path = os.path.join(
            os.path.dirname(__file__), "../data/residue_params.json"
        )
        with open(residue_params_path, "r", encoding="utf-8") as param_file:
            self.residue_data_dict = json.load(param_file)

        # read protein .pdb file using mdanalysis
        self.universe = mda.Universe(gmx_gro)

        # get atomic_numbers, use first letter of atom name as element
        self.atomic_numbers = [
            ELEMENTS_TO_ATOMIC_NUMBERS[a.name[0]] for a in self.universe.atoms
        ]

        # get residues-atoms mapping and inverse mapping for binding-site(pocket) discovery
        self.residues_atoms_mapping = {
            i: list(self.universe.residues.indices[i])
            for i in range(len(self.universe.residues.indices))
        }
        self.atoms_residues_mapping = {
            k: i for (i, j) in self.residues_atoms_mapping.items() for k in j
        }
        # get rotatable hbd H mask list
        rotatable_hbd_bond, _ = self.find_interaction("rotatable_hbd_bond")
        self.rotatable_h, _ = self.find_interaction("rotatable_h")

        if find_pocket:
            assert self.center is not None, "please input the center: [x0,y0,z0]"
            assert self.box is not None, "please input the box: [x1,y1,z1]"
            self.pocket_index = self.find_binding_site()
            rotatable_hbd_bond = [
                bond_index
                for bond_index in rotatable_hbd_bond
                if bond_index[0] in self.pocket_index
                or bond_index[1] in self.pocket_index
            ]
            self.rotatable_h = [i for i in self.rotatable_h if i in self.pocket_index]
        else:
            self.pocket_index = []

        self.is_rotatable = [
            1 if bond_index in rotatable_hbd_bond else 0
            for bond_index in self.bond_index
        ]
        assert sum(self.is_rotatable) == len(rotatable_hbd_bond)

        # get movable atoms, remove any ff idx that is completely frozen
        # TODO flexible residues
        if find_frozen:
            self.update_ffdata_mm_idx(self.rotatable_h)

        self.hydroph_atoms, _ = self.find_interaction("vina_hydrophobic")
        self.hbond_acc_atoms, _ = self.find_interaction("vina_acceptor")
        self.hbond_don_h_atoms, self.hbond_don_h_atom_names = self.find_interaction(
            "vina_donor"
        )  # return H as the acceptor
        self.element_check(self.hbond_don_h_atom_names, ["H"])
        self.cation_atoms, _ = self.find_interaction("cation_atoms")
        self.anion_atoms, _ = self.find_interaction("anion_atoms")
        self.pi_ring_atoms, _ = self.find_interaction("aromatic_ring")
        self.update_ffdata()

        # add additional information
        self.ffdata.update(
            deepcopy(
                {
                    "bond_index": self.bond_index,
                    "is_rotatable": self.is_rotatable,
                    "pocket_atomidx": self.pocket_index,
                }
            )
        )
        self.geometry_data = self.get_geometry_data(self.ffdata)

        # remove parm files
        if parm_output_dir is None:
            shutil.rmtree(os.path.abspath(self.working_dir))

        print("Initialize Done")

    def get_attach_atom_with_h(self, h_list, bond_list):
        H_donar_pairs = []
        h_idx = 0
        bidx = 0
        while h_idx < len(h_list) and bidx < len(bond_list):
            if h_list[h_idx] not in bond_list[bidx]:
                bidx += 1
            else:
                if h_list[h_idx] == bond_list[bidx][0]:
                    H_donar_pairs.append([bond_list[bidx][1], bond_list[bidx][0]])
                else:
                    H_donar_pairs.append([bond_list[bidx][0], bond_list[bidx][1]])
                h_idx += 1
                bidx += 1
        return H_donar_pairs

    def update_ffdata(self):
        # Unify the interface with small molecules. The commented parts are reserved interfaces.
        self.ffdata["atomic_numbers"] = self.atomic_numbers
        self.ffdata["hydrophobic_atomidx"] = self.hydroph_atoms
        self.ffdata["hbondacc_atomidx"] = self.hbond_acc_atoms
        self.ffdata["hbonddon_atomidx"] = self.get_attach_atom_with_h(
            self.hbond_don_h_atoms, self.bond_index
        )
        # self.ffdata["halogenacc_atomidx"] = self.halogen_acc_atoms
        # self.ffdata["halogendon_atomidx"] = self.halogen_don_atoms
        self.ffdata["cation_atomidx"] = self.cation_atoms
        self.ffdata["anion_atomidx"] = self.anion_atoms
        self.ffdata["piring5_atomidx"] = [x for x in self.pi_ring_atoms if len(x) == 5]
        self.ffdata["piring6_atomidx"] = [x for x in self.pi_ring_atoms if len(x) == 6]
        # self.ffdata["metaldonar_atomidx"] = self.metal_donar_atoms
        # self.ffdata["metalaccep_atomidx"] = self.metal_acc_atoms
        # subdivide hydrogen bonds based on charge
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
        self.ffdata["rotatable_h_atomidx"] = self.rotatable_h

    def update_ffdata_mm_idx(self, movable):
        def update_mm_idx_by_key(idx_key, param_keys=[]):
            if len(param_keys) > 0:
                for param_key in param_keys:
                    self.ffdata[param_key] = [
                        self.ffdata[param_key][i]
                        for i in range(len(self.ffdata[idx_key]))
                        if any(j in movable for j in self.ffdata[idx_key][i])
                    ]
            self.ffdata[idx_key] = [
                idx for idx in self.ffdata[idx_key] if any(j in movable for j in idx)
            ]
            if len(param_keys) > 0:
                for param_key in param_keys:
                    assert len(self.ffdata[param_key]) == len(self.ffdata[idx_key])

        update_mm_idx_by_key("FF_NonbondedAll_atomidx")
        update_mm_idx_by_key("FF_Nonbonded14_atomidx")
        update_mm_idx_by_key("FF_Bonds_atomidx", ["FF_Bonds_k", "FF_Bonds_length"])
        update_mm_idx_by_key("FF_Angles_atomidx", ["FF_Angles_k", "FF_Angles_angle"])
        update_mm_idx_by_key(
            "FF_ProperTorsions_atomidx",
            [
                "FF_ProperTorsions_periodicity",
                "FF_ProperTorsions_phase",
                "FF_ProperTorsions_k",
            ],
        )
        update_mm_idx_by_key(
            "FF_ImproperTorsions_atomidx",
            [
                "FF_ImproperTorsions_periodicity",
                "FF_ImproperTorsions_phase",
                "FF_ImproperTorsions_k",
            ],
        )

    def pdb_to_ffdata(self, working_dir, pdbfile, forcefield):

        working_dir = os.path.abspath(working_dir)
        pdbfile = os.path.abspath(pdbfile)

        # step 0: clean protein
        step0_pdb = os.path.join(working_dir, "step0_cleaned.pdb")
        generate_sanitized_protein_pdb(pdbfile, working_dir, "step0_cleaned.pdb")

        # step 1, 2 convert to parmed object
        parm = pdb_to_parm(working_dir, step0_pdb, forcefield)

        # step 3 convert to gromacs format
        gmx_gro = os.path.join(working_dir, "step3.gro")
        parm.save(gmx_gro, overwrite=True)

        # step 4: parse ffdata from parmed object
        logger.info("parsing ffdata and xyz from parmed object")
        ffdata, xyz = parm_to_ffdata(parm)
        # transform to list
        ffdata = convert_value_to_list(ffdata)
        logger.info("parse ffdata from parmed object finished")
        return gmx_gro, ffdata, xyz

    def get_geometry_data(self, re_ffdata):
        receptor_geometry = ReceptorGeometry(
            bond_index=self.bond_index,
            is_rotatable=self.is_rotatable,
            num_atoms=len(re_ffdata["atomic_numbers"]),
        )
        num_atoms_value = (
            receptor_geometry.num_atoms
            if isinstance(receptor_geometry.num_atoms, int)
            else receptor_geometry.num_atoms.tolist()
        )
        geometry_data = {
            "num_atoms": num_atoms_value,
            # "flexible_atoms": receptor_geometry.flexible_atoms,
            "num_rotatable_bonds": receptor_geometry.num_rotatable_bonds,
            "data": receptor_geometry.data,
            # "frag_traverse_levels": receptor_geometry.frag_traverse_levels,
        }
        return geometry_data

    @staticmethod
    def element_check(atom_list, element_list):
        """_summary_

        Args:
            atom_list (_type_): _description_
            element (_type_): _description_
        """
        for atom_name in atom_list:
            # TODO: only supports the first character of the element
            atom_element = atom_name[0]
            assert (
                atom_element in element_list
            ), f"Element {atom_element} and atom name do not match"

    def load_data(self, receptor_fpath):
        """_summary_

        Args:
            receptor_fpath (_type_): _description_
        """
        with open(receptor_fpath, "r", encoding="utf-8") as f:
            receptor_data = json.load(f)

        self.gmx_gro_string = receptor_data["gmx_gro_string"]
        self.ffdata = receptor_data["ffdata"]
        self.xyz = np.array(receptor_data["xyz"])  # np.array
        self.bond_index = receptor_data["bond_index"]
        self.is_rotatable = receptor_data["is_rotatable"]
        self.pocket_index = receptor_data["pocket_atomidx"]

    def _convert_idx_aname(self, atom_name_list, local_index_list, global_mapping_dict):
        return sorted(
            [
                global_mapping_dict[atom_name_list[local_index]]
                for local_index in local_index_list
            ]
        )

    def find_interaction(self, interaction_name):
        """_summary_

        Args:
            interaction_name (_type_): _description_

        Returns:
            _type_: _description_
        """
        interaction_idx_list, interact_atom_name_list = [], []
        for residue in self.universe.residues:
            resname = residue.resname
            atom_index_mapping_dict = dict(zip(residue.atoms.names, residue.atoms.ix))
            atom_name_mapping_dict = dict(zip(residue.atoms.ix, residue.atoms.names))
            if "H1" in atom_index_mapping_dict and (resname not in ["ACE", "NME"]):
                resname = "N" + resname

            if "OXT" in atom_index_mapping_dict:
                resname = "C" + resname
            residue_data = self.residue_data_dict[resname]
            atom_name_list = residue_data["atom_names"]
            if (
                interaction_name not in residue_data
            ):  # skip if no interaction in the params file
                continue
            if interaction_name in LIST_PARAMS:
                for atom_list in residue_data[interaction_name]:
                    global_mapped_atom_list = self._convert_idx_aname(
                        atom_name_list, atom_list, atom_index_mapping_dict
                    )
                    interaction_idx_list.append(global_mapped_atom_list)
                    interact_atom_name_list.append(
                        atom_name_mapping_dict[i] for i in global_mapped_atom_list
                    )
            elif interaction_name in INDEX_PARAMS:
                data_list = residue_data[interaction_name]
                atom_list = [i for i in range(len(data_list)) if data_list[i]]
                global_mapped_atom_list = self._convert_idx_aname(
                    atom_name_list, atom_list, atom_index_mapping_dict
                )
                interaction_idx_list += global_mapped_atom_list
                interact_atom_name_list += [
                    atom_name_mapping_dict[i] for i in global_mapped_atom_list
                ]
        return interaction_idx_list, interact_atom_name_list

    def find_binding_site(self):
        atom_to_center_vector = self.xyz - np.expand_dims(np.array(self.center), axis=0)
        atom_in_box_flag = np.all(
            abs(atom_to_center_vector)
            < np.expand_dims(np.array(self.box), axis=0) / 2.0,
            axis=1,
        )
        atom_in_box_index = list(np.where(atom_in_box_flag)[0])
        residues_in_box_index = set(
            [self.atoms_residues_mapping[x] for x in atom_in_box_index]
        )
        refine_atom_in_box_index = [
            j for x in residues_in_box_index for j in self.residues_atoms_mapping[x]
        ]
        return refine_atom_in_box_index

    def add_cofactor(self, cofactor_data):
        # merge two .gro file
        cofactor_gro_string = cofactor_data["gmx_gro_string"]
        cofactor_universe = mda.Universe(StringIO(cofactor_gro_string), format="GRO")
        receptor_universe = mda.Universe(StringIO(self.gmx_gro_string), format="GRO")
        merged_universe = mda.Merge(receptor_universe.atoms, cofactor_universe.atoms)

        # update merged universe
        merged_gro_path = os.path.join(self.working_dir, "receptor_cofactor_merged.gro")
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        with mda.Writer(merged_gro_path, reindex=False, format="GRO") as gro_write:
            gro_write.write(merged_universe.atoms)
        with open(merged_gro_path, "r", encoding="utf-8") as gro_file:
            self.gmx_gro_string = gro_file.read()

        receptor_atom_num = receptor_universe.atoms.n_atoms
        receptor_atomidx = sorted(
            list(set([atom[0] for atom in self.ffdata["FF_vdW_atomidx"]]))
        )
        cofactor_atom_num = cofactor_universe.atoms.n_atoms

        old_receptor_data = deepcopy(self.get_data())
        # update ff data
        overlap_keys = set(old_receptor_data["ffdata"].keys()) & set(
            cofactor_data["ffdata"].keys()
        )
        assert (
            len(overlap_keys) == 35
        ), f"overlap keys' number error: {len(overlap_keys)} v.s. 35"

        # update FF_vdW_paraidx
        self.ffdata["FF_vdW_paraidx"].extend(cofactor_data["ffdata"]["FF_vdW_paraidx"])

        # updata atomidx
        # TODO: update rotatable_h_atomidx later
        update_keys = [key for key in overlap_keys if "atomidx" in key]

        assert (
            len(update_keys) == 18
        ), f"atomidx keys' number error : {len(update_keys)} v.s. 18"
        for update_key in update_keys:
            shift_num = receptor_atom_num
            new_data = np.arrray(cofactor_data["ffdata"][update_key]) + shift_num
            self.ffdata[update_key].extend(new_data.tolist())

        # update pair between old and new
        cofactor_atomidx = sorted(
            list(
                set(
                    [
                        atom[0] + receptor_atom_num
                        for atom in cofactor_data["ffdata"]["FF_vdW_atomidx"]
                    ]
                )
            )
        )
        new_nobond_pairs = [
            [receptor_idx, cofactor_idx]
            for receptor_idx in receptor_atomidx
            for cofactor_idx in cofactor_atomidx
        ]
        self.ffdata["FF_NonbondedAll_atomidx"].extend(new_nobond_pairs)

        # update FF parameters
        update_keys = [key for key in overlap_keys if "idx" not in key]
        assert (
            len(update_keys) == 16
        ), f"FF paramters number error : {len(update_keys)} v.s. 16"
        for update_key in update_keys:
            new_data = cofactor_data["ffdata"][update_key]
            self.ffdata[update_key].extend(new_data)

        # update xyz
        self.xyz = np.concatenate([self.xyz, cofactor_data["xyz"][0]], axis=0)

        # update bond_index
        shift_num = receptor_atom_num
        new_bond_index = np.array(cofactor_data["bond_index"]) + shift_num
        self.bond_index.extend(new_bond_index.tolist())

        # update is_rotatable
        self.is_rotatable.extend([0] * len(cofactor_data["is_rotatable"]))

        # update pocket index
        if self.pocket_index:
            new_data = [receptor_atom_num + i for i in range(1, cofactor_atom_num + 1)]
            self.pocket_index.extend(new_data)

        self.geometry_data = self.get_geometry_data(self.ffdata)

        self.ffdata.update(
            deepcopy(
                {
                    "bond_index": self.bond_index,
                    "is_rotatable": self.is_rotatable,
                    "pocket_atomidx": self.pocket_index,
                }
            )
        )

    def get_data(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        return convert_value_to_list(
            {
                "gmx_gro_string": self.gmx_gro_string,
                "ffdata": self.ffdata,
                "xyz": self.xyz,
                "bond_index": self.bond_index,
                "is_rotatable": self.is_rotatable,
                "pocket_atomidx": self.pocket_index,
                "geometry": self.geometry_data,
            }
        )
