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

from itertools import accumulate
from typing import List

import numpy as np

from pxdock.common import get_logger
from pxdock.geometry.rotation import (
    euler_rotation_matrix,
    rodrigues_rotation_matrix,
)
from pxdock.geometry.tree import RigidFragmentTree

logger = get_logger(__name__)


class BaseGeometry:
    def __init__(
        self,
        bond_index: List[List[int]],
        is_rotatable: List[int],
        num_atoms: int = None,
    ):
        super().__init__()
        ### build rigid fragment tree ###
        self.frag_tree = RigidFragmentTree(
            bond_index, is_rotatable, num_atoms=num_atoms
        )
        self.frag_traverse_levels = self.frag_tree.levels()
        self.num_rotatable_bonds = self.frag_tree.number_of_edges
        self.num_atoms = self.frag_tree.number_of_atoms


class ReceptorGeometry(BaseGeometry):
    def __init__(
        self,
        bond_index: List[List[int]],
        is_rotatable: List[int],
        num_atoms: int = None,
    ):
        """
        A class manipulate receptor coordinates with movements in flexible residue torsional rotations.

        Args:
        - bond_index (List[List[int]]): bond index of receptor.
        - is_rotatable (List[int]): Each element indicates whether the bond is rotatable.
        - num_atoms (int, optional): number of atoms. Defaults to None.
        """
        super().__init__(bond_index, is_rotatable, num_atoms)

        self.data = {}
        if self.num_rotatable_bonds == 0:
            self.flexible_atoms = []
            return

        assert (
            len(self.frag_traverse_levels) == 2
        ), f"The receptor fragment tree has level: {len(self.frag_traverse_levels)}"

        self.leaf_level = self.frag_traverse_levels[1]

        for frag_size in [2, 3, 4]:
            fragments = []
            rot_bond_1 = []
            rot_bond_2 = []
            rot_bond_order = []

            for i, data_i in self.leaf_level.items():
                if self.frag_tree.fragment_size(i) != frag_size:
                    continue
                fragments.append(self.frag_tree.fragment(i))
                edge = self.frag_tree.edge(data_i["parent"][0], i)
                bond_atom_1, bond_atom_2 = edge["atoms"]
                rot_bond_1.append(bond_atom_1)
                rot_bond_2.append(bond_atom_2)
                rot_bond_order.append(edge["order"])

            if fragments:
                self.data[frag_size] = {
                    "fragments": fragments,
                    "rot_bond_1": rot_bond_1,
                    "rot_bond_2": rot_bond_2,
                    "rot_bond_order": rot_bond_order,
                }
        self.flexible_atoms = sorted(
            sum([sum(data["fragments"], []) for data in self.data.values()], [])
        )

        # make sure all fragments have size in [2,3,4]
        assert len(self.leaf_level) == sum(
            [len(data["rot_bond_order"]) for data in self.data.values()]
        )

    def torsion_to_xyz(self, torsion: np.ndarray, xyz0: np.ndarray):
        """
        transform flexible receptor residue torsion coordinates to xyz coordinates.

        Args:
        - torsion (np.ndarray): torsion angles changes.
        - xyz0 (np.ndarray): original xyz coordinates of the receptor.
        """
        if self.num_rotatable_bonds == 0:
            return xyz0.copy()
        xyz_out = []

        atoms = []
        for frag_size, data in self.data.items():
            fragments = np.array(data["fragments"]).astype(int)
            rot_bond_1, rot_bond_2, rot_bond_order = (
                data["rot_bond_1"],
                data["rot_bond_2"],
                data["rot_bond_order"],
            )
            atom_index = fragments.reshape(-1).tolist()
            atoms.append(atom_index)

            xyz_frag = xyz0[..., atom_index, :]
            xyz_frag = xyz_frag.reshape(
                xyz_frag.shape[:-2] + fragments.shape + (3,)
            )  # [..., num_fragments, num_atoms_in_fragment, 3]
            rot_bond_xyz1 = xyz0[..., rot_bond_1, :]  # [..., num_fragments, 3]
            rot_bond_xyz2 = xyz0[..., rot_bond_2, :]  # [..., num_fragments, 3]
            theta = torsion[..., rot_bond_order]  # [..., num_fragments]

            u = rot_bond_xyz2 - rot_bond_xyz1
            u = u / np.clip((u**2).sum(axis=-1), 1e-5, None)**0.5

            # rotation matrix
            R = rodrigues_rotation_matrix(u, theta)  # [..., num_fragments, 3, 3]
            # torsion
            out = (
                np.einsum(
                    "...ij,...kj->...ki",
                    R,
                    (xyz_frag - rot_bond_xyz1.unsqueeze(axis=-2)),
                )
                + rot_bond_xyz1.unsqueeze(axis=-2).detach()
            )  # [..., num_fragments, num_atoms_in_fragment, 3]
            out = out.reshape(
                out.shape[:-3] + (len(atom_index), 3)
            )  # [..., num_fragments * num_atoms_in_fragment, 3]

            xyz_out.append(out)

        xyz_out = np.concatenate(xyz_out, axis=-2)
        atom_index = sum(atoms, [])
        xyz = xyz0.clone().detach()
        xyz[..., atom_index, :] = xyz_out

        return xyz


class LigandGeometry(BaseGeometry):
    def __init__(
        self,
        bond_index: List[List[int]],
        is_rotatable: List[int],
        num_atoms: int = None,
    ):
        """
        A class manipulate ligand coordinates with movements in tranlation, rotation, and torsional rotation along rotatable bonds.

        Args:
        - bond_index (List[List[int]]): bond index of ligand.
        - is_rotatable (List[int]): is_rotatable of ligand.
        - num_atoms (int, optional): number of atoms. Defaults to None.
        """
        super().__init__(bond_index, is_rotatable, num_atoms)

        self.fragments = self.frag_tree.fragments
        self.frag_split_index = list(accumulate([len(frag) for frag in self.fragments]))
        n = self.frag_split_index.pop(-1)
        assert n == self.num_atoms

        self.permute_index = sum(self.fragments, [])
        index_map = dict(zip(self.permute_index, list(range(self.num_atoms))))
        self.reverse_index = [index_map[i] for i in range(self.num_atoms)]

    def permute_xyz(self, xyz: np.ndarray, reverse: bool = False):
        """
        Permute the coordinates of the atoms according to the corresponding fragment ordering.

        Args:
        - xyz: np.ndarray of shape [..., natom, 3], coordinates of ligand.
        - reverse: bool, if True, permute the coordinates back to the original order.
        """
        n = xyz.shape[-2]
        assert n == len(self.permute_index)
        if reverse:
            index_map = dict(zip(self.permute_index, list(range(n))))
            permute_index = [index_map[i] for i in range(n)]
        else:
            permute_index = self.permute_index
        return xyz[..., permute_index, :]

    def torsion6_to_xyz(
        self,
        torsion: np.ndarray,
        translaton: np.ndarray,
        rotation: np.ndarray,
        xyz0: np.ndarray,
        rot_center: str = "first_atom",
    ):
        """
        transform ligand torsion6 coordinates to xyz coordinates.

        Args:
        - torsion: np.ndarray of shape [..., num_rbonds], torsion angles changes.
        - translaton: np.ndarray of shape [..., 3], translation vector changes.
        - rotation: np.ndarray of shape [..., 3], rotation angles changes.
        - xyz0: np.ndarray of shape [..., natom, 3], original coordinates of ligand.
        - rot_center: str, optional, the center of rotation, can be "first_atom" or "geo_center".
        """
        xyz0_permuted = self.permute_xyz(xyz0)

        ### (1) Rotation ###
        rot_matrix = euler_rotation_matrix(angles=rotation)  # [num_poses, 3, 3]
        if rot_center == "first_atom":
            rot_center = xyz0_permuted[
                ..., 0:1, :
            ]  # the first atom in the permuted xyz must belong to the root frame, and has the smallest atom index in the root frame
        elif rot_center == "geo_center":
            rot_center = np.mean(xyz0_permuted, axis=1, keepdims=True)
        else:
            rot_center = xyz0_permuted[
                ..., 0:1, :
            ]  # the first atom in the permuted xyz must belong to the root frame, and has the smallest atom index in the root frame
        xyz_centered = xyz0_permuted - rot_center
        xyz = np.einsum("...ij,...kj->...ki", rot_matrix, xyz_centered)

        ### (2) Transition ###
        if len(translaton.shape) < len(xyz.shape):
            translaton = translaton.unsqueeze(axis=-2)
        xyz = xyz + rot_center + translaton

        ### (3) Torsion ###
        xyz_fragments = np.split(xyz, self.frag_split_index, axis=-2)

        # rotation matrix
        RotationMatrix = {}
        # parent node
        parent_xyz = {}

        for _, level in enumerate(self.frag_traverse_levels):
            # TODO: Fragments at the same level can be processed in parallel
            for i, frag_data in level.items():
                children_list, parent = frag_data["children"], frag_data["parent"]
                if len(parent) > 0:
                    # rotation matrix of parent
                    parent = list(parent)[0]
                    if parent in RotationMatrix:
                        R_parent = RotationMatrix[parent]
                    else:
                        assert parent in self.frag_traverse_levels[0].keys()
                        R_parent = None

                    # rotatable bond between parent and i
                    edge_data = self.frag_tree.edge(parent, i)
                    parent_atom, atom = edge_data["atoms"]
                    atom_xyz0, parent_xyz0 = (
                        xyz[..., self.permute_index.index(atom), :],
                        xyz[..., self.permute_index.index(parent_atom), :],
                    )
                    u = atom_xyz0 - parent_xyz0
                    if R_parent is not None:
                        u = np.einsum("...ij,...j->...i", R_parent, u)
                    u = u / np.clip((u**2).sum(axis=-1), 1e-5, None)**0.5

                    # new rotation matrix
                    edge_order = edge_data["order"]
                    theta = torsion[..., edge_order]  # [...]
                    R = rodrigues_rotation_matrix(u, theta)
                    if R_parent is not None:
                        R = np.einsum("...ij,...jk", R, R_parent)
                    RotationMatrix[i] = R

                    # perform torsion
                    xyz_fragments[i] = np.einsum(
                        "...ij,...kj->...ki",
                        R,
                        xyz_fragments[i] - parent_xyz0.unsqueeze(axis=-2),
                    ) + parent_xyz[i].unsqueeze(axis=-2)

                # save location for its children
                for child in children_list:
                    # bond
                    edge_data = self.frag_tree.edge(i, child)
                    atom, _ = edge_data["atoms"]
                    parent_xyz[child] = xyz_fragments[i][
                        ..., self.fragments[i].index(atom), :
                    ]

        xyz_rotated = np.concatenate(xyz_fragments, axis=-2)
        # permute the index back
        return self.permute_xyz(xyz_rotated, reverse=True)


def vec_distance_sqr(a, b):
    """
    Calculate the square of the distance between two 3-D points.

    Args:
    - a: [3] the first point's coordinates (a list or tuple of three elements).
    - b: [3] the second point's coordinates (a list or tuple of three elements).

    Returns:
    - float: the square of the distance between the two points.
    """
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
