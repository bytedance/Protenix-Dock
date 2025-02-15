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

# Some of the SMARTS patterns used in the interaction classes are derived from ProLIF
# which is licensed under Apache License 2.0.
#
# see https://github.com/chemosim-lab/ProLIF/blob/master/prolif/interactions/interactions.py

import json
import os
from collections import namedtuple

import numpy as np
from rdkit import Chem

from pxdock.dock.utils import kScoreConfigFile


def cluster_atoms(mol, atom_indices):
    # Initialize a set for visited atoms to avoid redundancy
    visited = set()
    # Dictionary to hold clusters with atom index as key and cluster number as value
    atom_cluster = {}
    # Initialize cluster count
    cluster_count = 0

    for idx in atom_indices:
        # If the atom has not been visited, start a new cluster
        if idx not in visited:
            # Increase cluster count
            cluster_count += 1
            # Depth-First Search to find all connected atoms
            stack = [idx]
            while stack:
                atom_idx = stack.pop()
                if atom_idx not in visited:
                    # Assign atom to the current cluster
                    atom_cluster[atom_idx] = cluster_count
                    # Mark this atom as visited
                    visited.add(atom_idx)
                    # Get all connected atoms
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neigh_idx = neighbor.GetIdx()
                        # If neighbor is in the initial atom indices and hasn't been visited, add to stack
                        if neigh_idx in atom_indices and neigh_idx not in visited:
                            stack.append(neigh_idx)
    # Group atom indices by cluster
    clusters = {}
    for atom_idx, cluster_id in atom_cluster.items():
        clusters.setdefault(cluster_id, []).append(atom_idx)
    return clusters


class TorsionStrainParams:
    def __init__(self, json_fpath: str = None, task="pose_selection"):
        """
        class storing and getting torsion strain parameters

        Args:
            json_fpath (str, optional): json file path. Defaults to None. If None, use default params.
        """
        assert task in ["pose_selection", "affinity_ranking"]
        json_fpath = json_fpath or kScoreConfigFile
        self.load_params_from_json(json_fpath, task=task)
        self.num_pattern = len(self.params_list)

    def load_params_from_json(self, param_fpath: str, task="pose_selection"):
        assert os.path.exists(param_fpath), f"{param_fpath} not exist"
        param_cls = namedtuple(
            "TorsionStrainPattern", ["name", "paraidx", "smarts", "params"]
        )
        with open(param_fpath, "r") as f:
            data = json.load(f)[task]["torsion_strain_penalty"]["types"]
        self.params_list = []
        for paraidx, dd in enumerate(data):
            param = param_cls(**dd)
            assert paraidx == param.paraidx
            self.params_list.append(param)

    def get_params_by_idx(self, paraidx: int):
        assert paraidx < self.num_pattern
        params = self.params_list[paraidx].params
        return params


class MolBase:
    def __init__(self, rkmol, **kwargs) -> None:
        # Find specific atom types or chemical environments
        # self.all_atoms: List[rdkit.Chem.rdchem.Atom] same order with mapped_smiles
        self.rkmol = rkmol
        self.all_atoms = list(self.rkmol.GetAtoms())
        self.hydroph_atoms = self.find_hyda(self.rkmol, self.all_atoms)
        self.hbond_acc_atoms = self.find_hba(self.rkmol, self.all_atoms)
        self.hbond_don_atoms = self.find_hbd(self.rkmol, self.all_atoms)
        self.halogen_acc_atoms = self.find_xba(self.rkmol, self.all_atoms)
        self.halogen_don_atoms = self.find_xbd(self.rkmol, self.all_atoms)
        self.cation_atoms = self.find_cationa(self.rkmol, self.all_atoms)
        self.anion_atoms = self.find_aniona(self.rkmol, self.all_atoms)
        self.pi_ring_atoms = self.find_pi_rings(self.rkmol, self.all_atoms)
        self.metal_donar_atoms = self.find_metda(self.rkmol, self.all_atoms)
        self.metal_acceptor_atoms = self.find_metaa(self.rkmol, self.all_atoms)
        score_function_json = kwargs.get("score_function_json", None)
        self.torsion_strain_bonds = {
            "pose_selection": self.find_torsion_strain_bonds(
                self.rkmol, score_function_json, "pose_selection"
            ),
            "affinity_ranking": self.find_torsion_strain_bonds(
                self.rkmol, score_function_json, "affinity_ranking"
            ),
        }

    @staticmethod
    def find_hyda(rkmol, all_atoms):
        """Select all carbon atoms which have only carbons and/or hydrogens as direct neighbors."""
        atom_set = []
        data = namedtuple("hydrophobic", "atom idx")
        hydrophobic_smarts = (
            "[c,s,Br,I,S&H0&v2," "$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6X4H0]);+0]"
        )

        patt = Chem.MolFromSmarts(hydrophobic_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(atom=all_atoms[lig_match[0]], idx=all_atoms[lig_match[0]].GetIdx())
            )
        return atom_set

    @staticmethod
    def find_hba(rkmol, all_atoms):
        """Find all possible hydrogen bond acceptors"""
        data = namedtuple("hbondacceptor", "atom idx")
        atom_set = []
        acceptor_smarts = (
            "[#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),"
            "O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0]"
        )

        patt = Chem.MolFromSmarts(acceptor_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=all_atoms[lig_match[0]],
                    idx=all_atoms[lig_match[0]].GetIdx(),
                )
            )
        return atom_set

    @staticmethod
    def find_hbd(rkmol, all_atoms):
        """Find all possible hydrogen bond donar"""
        data = namedtuple("hbonddonar", "atom idx")
        atom_set = []
        donar_smarts = "[$([O,S;+0]),$([N;v3,v4&+1]),n+0]-[H]"
        patt = Chem.MolFromSmarts(donar_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=[all_atoms[lig_match[0]], all_atoms[lig_match[1]]],
                    idx=[
                        all_atoms[lig_match[0]].GetIdx(),
                        all_atoms[lig_match[1]].GetIdx(),
                    ],
                )
            )
        return atom_set

    @staticmethod
    def find_xba(rkmol, all_atoms):
        """Find all possible Halogen bond acceptors"""
        data = namedtuple("halogenacceptor", "atom idx")
        atom_set = []
        acceptor_smarts = "[#7,#8,P,S,Se,Te,a;!+{1-}][*]"
        patt = Chem.MolFromSmarts(acceptor_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=[all_atoms[lig_match[0]], all_atoms[lig_match[1]]],
                    idx=[
                        all_atoms[lig_match[0]].GetIdx(),
                        all_atoms[lig_match[1]].GetIdx(),
                    ],
                )
            )
        return atom_set

    @staticmethod
    def find_xbd(rkmol, all_atoms):
        """Find all possible Halogen bond donar"""
        data = namedtuple("halogendonar", "atom idx")
        atom_set = []
        donar_smarts = "[#6,#7,Si,F,Cl,Br,I]-[Cl,Br,I,At]"
        patt = Chem.MolFromSmarts(donar_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=[all_atoms[lig_match[0]], all_atoms[lig_match[1]]],
                    idx=[
                        all_atoms[lig_match[0]].GetIdx(),
                        all_atoms[lig_match[1]].GetIdx(),
                    ],
                )
            )
        return atom_set

    @staticmethod
    def find_cationa(rkmol, all_atoms):
        """Find all possible cationic atom"""
        data = namedtuple("cationa", "atom idx")
        atom_set = []
        cation_smarts = "[+{1-},$([NX3&!$([NX3]-O)]-[C]=[NX3+])]"
        patt = Chem.MolFromSmarts(cation_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=all_atoms[lig_match[0]],
                    idx=all_atoms[lig_match[0]].GetIdx(),
                )
            )
        return atom_set

    @staticmethod
    def find_aniona(rkmol, all_atoms):
        """Find all possible anionaic atom"""
        data = namedtuple("aniona", "atom idx")
        atom_set = []
        anion_smarts = "[-{1-},$(O=[C,S,P]-[O-])]"
        patt = Chem.MolFromSmarts(anion_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=all_atoms[lig_match[0]],
                    idx=all_atoms[lig_match[0]].GetIdx(),
                )
            )
        return atom_set

    @staticmethod
    def find_pi_rings(rkmol, all_atoms):
        """SMARTS for aromatic rings (5 and 6 membered rings only)"""
        data = namedtuple("pi_ring", "atom idx type")
        atom_set = []
        pi_ring_smarts = (
            "[a;r6]1:[a;r6]:[a;r6]:[a;r6]:[a;r6]:[a;r6]:1",
            "[a;r5]1:[a;r5]:[a;r5]:[a;r5]:[a;r5]:1",
        )
        for idx, pi_ring_smart in enumerate(pi_ring_smarts):
            if idx == 0:
                type = "six_mem_ring"
            if idx == 1:
                type = "five_mem_ring"
            patt = Chem.MolFromSmarts(pi_ring_smart)
            hit_ats = rkmol.GetSubstructMatches(patt)
            for lig_match in hit_ats:
                atom_set.append(
                    data(
                        atom=[all_atoms[i] for i in lig_match],
                        idx=[all_atoms[i].GetIdx() for i in lig_match],
                        type=type,
                    )
                )
        return atom_set

    @staticmethod
    def find_metda(rkmol, all_atoms):
        """Find Metal donar atoms"""
        data = namedtuple("metal_donar", "atom idx")
        atom_set = []
        metal_donar_smarts = "[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]"
        patt = Chem.MolFromSmarts(metal_donar_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=all_atoms[lig_match[0]],
                    idx=all_atoms[lig_match[0]].GetIdx(),
                )
            )

        return atom_set

    @staticmethod
    def find_metaa(rkmol, all_atoms):
        """Find Metal acceptor atoms"""
        data = namedtuple("metal_acceptor", "atom idx")
        atom_set = []
        metal_acceptor_smarts = (
            "[O,#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4]),-{1-};!+{1-}]"
        )
        patt = Chem.MolFromSmarts(metal_acceptor_smarts)
        hit_ats = rkmol.GetSubstructMatches(patt)
        for lig_match in hit_ats:
            atom_set.append(
                data(
                    atom=all_atoms[lig_match[0]],
                    idx=all_atoms[lig_match[0]].GetIdx(),
                )
            )

        return atom_set

    @staticmethod
    def find_torsion_strain_bonds(
        rkmol, score_function_json=None, task="pose_selection"
    ):
        """Find all possible bonds which maybe penalized by torsion strain"""

        tstrain_ffdata = namedtuple(
            "tstrain_ffdata",
            ["pattern_name", "paraidx", "bond_tag", "proper_atomidx"],
        )
        penalty_bonds = dict()
        coords = [m.GetPositions() for m in rkmol.GetConformers()]  # [nconf, natom, 3]
        coords = np.array(coords)
        torsion_strain_patterns = TorsionStrainParams(
            score_function_json, task
        ).params_list

        for paraidx, pattern in enumerate(torsion_strain_patterns):
            patt_mol = Chem.MolFromSmarts(pattern.smarts)
            idx_map = dict()
            for atom in patt_mol.GetAtoms():
                # this checks the map y in pattern [#x:y].
                # y is zero if not defined in pattern, like [#x]
                map_index = atom.GetAtomMapNum()
                if map_index != 0:
                    idx_map[map_index - 1] = atom.GetIdx()
            map_list = [idx_map[x] for x in sorted(idx_map)]
            full_matches = rkmol.GetSubstructMatches(patt_mol, uniquify=False)
            full_matches = set(
                [tuple(match[x] for x in map_list) for match in full_matches]
            )
            # [n_match, 4]

            # Note 1: there may be multiple different proper matches correspond to same bond
            # when the bonds related with multiple propers, the penalty will be average of these proper penalties

            proper_set = set()
            bond_proper_match = dict()
            for proper_atomids in full_matches:
                bond_tag = "_".join(map(str, sorted(proper_atomids[1:3])))

                # Note 2: When uniquify in GetSubstructMatches() is set to be False,
                # sometimes identical proper will be returned in hitpropers, which depends on
                # how you write your smarts pattern. Therefore deduplication is needed here.

                proper_tag = "_".join(map(str, proper_atomids))
                if proper_tag in proper_set:
                    continue
                proper_set.add(proper_tag)
                bond_proper_match.setdefault(bond_tag, []).append(proper_atomids)

            for bond_tag, proper_atomids_list in bond_proper_match.items():
                # Note 3: Max number of different propers match to single bond is defaultly set to 4.
                # the rest padding values are zeros.

                # Note 4: Same bond will never be penaltized twise.
                # If single bonds matches multiple DIFFERENT torsion penalty patterns,
                # only the latest matched pattern (i.e., larger paraidx) will be kept.

                penalty_bonds[bond_tag] = tstrain_ffdata(
                    pattern_name=pattern.name,
                    paraidx=paraidx,
                    bond_tag=bond_tag,
                    proper_atomidx=proper_atomids_list,
                )

        return list(penalty_bonds.values())

    @property
    def get_hydrophobic_atoms(self):
        return self.hydroph_atoms

    @property
    def get_hba(self):
        return self.hbond_acc_atoms

    @property
    def get_hbd(self):
        return self.hbond_don_atoms

    @property
    def get_xba(self):
        return self.halogen_acc_atoms

    @property
    def get_xbd(self):
        return self.halogen_don_atoms

    @property
    def get_cationa(self):
        return self.cation_atoms

    @property
    def get_aniona(self):
        return self.anion_atoms

    @property
    def get_pi_ring(self):
        return self.pi_ring_atoms

    @property
    def get_metal_donar(self):
        return self.metal_donar_atoms

    @property
    def get_metal_acceptor(self):
        return self.metal_acceptor_atoms

    @property
    def get_torsion_strain_bonds(self):
        return self.torsion_strain_bonds
