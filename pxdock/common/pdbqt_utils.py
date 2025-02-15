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

# This file contains portions of code derived from meeko-0.3.3,
# which is licensed under the Apache License, Version 2.0 (the "License").

# You may obtain a copy of the License at:
#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Portions of the code in this file have been modified by ByteDance on 2025-02-10.

from meeko import MoleculePreparation, PDBQTMolecule
from rdkit import Chem
from rdkit.Geometry import Point3D


class MyPDBQTMolecule(PDBQTMolecule):
    def export_rdkit_mol(self, exclude_hydrogens=False, return_index_map=False):
        """
        Part of the code below is derived from the PDBQTMolecule.export_rdkit_mol which is licensed under the Apache License, Version 2.0.

        Main changes compared to the parent class:
        1. Sort the h_parents, fix the order of h_parents. Ensure that the positions of hydrogen atoms with known coordinates are consistently overwritten each time hydrogens are added.
        2. Return pdbqt_index_to_rdkit_index, which saves the atomic index mapping.
        3. Set RDKit's removeHs to False, otherwise, the atomic order may be problematic.
        """

        smiles = self._pose_data["smiles"]  # REMARK SMILES
        if smiles is None:
            raise RuntimeError("need REMARK SMILES in PDBQT file")
        index_map = self._pose_data["smiles_index_map"]  # REMARK SMILES IDS
        h_parents = self._pose_data["smiles_h_parent"]  # REMARK H PARENT
        # sort h_parents by x[0] and then by x[1]
        h_parents = sorted(h_parents, key=lambda x: (x[0], x[1]))

        # set removeHs=False
        params = Chem.SmilesParserParams()
        params.removeHs = False
        mol = Chem.MolFromSmiles(smiles, params)
        n_atoms = mol.GetNumAtoms()
        # https://sourceforge.net/p/rdkit/mailman/message/36474923/
        conf = Chem.Conformer(n_atoms)
        n_matched_atoms = 0

        # init index map
        pdbqt_index_to_rdkit_index = {}

        for i in range(n_atoms):
            pdbqt_index = index_map[i + 1] - 1
            x, y, z = self._positions[self._current_pose][pdbqt_index, :]
            conf.SetAtomPosition(i, Point3D(x, y, z))
            pdbqt_index_to_rdkit_index[pdbqt_index] = i

        conf_id = mol.AddConformer(conf)
        if exclude_hydrogens:
            return mol
        mol = Chem.AddHs(mol, explicitOnly=False, addCoords=True)
        conf = mol.GetConformer()
        used_h = []
        for parent_rdkit_index, h_pdbqt_index in h_parents:
            h_pdbqt_index -= 1
            x, y, z = self._positions[self._current_pose][h_pdbqt_index, :]
            parent_atom = mol.GetAtomWithIdx(parent_rdkit_index - 1)
            candidate_hydrogens = [
                atom.GetIdx()
                for atom in parent_atom.GetNeighbors()
                if atom.GetAtomicNum() == 1
            ]
            if len(candidate_hydrogens) < 1:
                raise ValueError("H parent does not have H neighbours.")
            for h_rdkit_index in candidate_hydrogens:
                if h_rdkit_index not in used_h:
                    break
            used_h.append(h_rdkit_index)
            conf.SetAtomPosition(h_rdkit_index, Point3D(x, y, z))
            pdbqt_index_to_rdkit_index[h_pdbqt_index] = h_rdkit_index
        # sort by pdbqt_index
        pdbqt_index_to_rdkit_index = dict(
            sorted(pdbqt_index_to_rdkit_index.items(), key=lambda x: x[0])
        )
        if not return_index_map:
            return mol
        else:
            return mol, pdbqt_index_to_rdkit_index


class MyMoleculePreparation(MoleculePreparation):
    def write_pdbqt_string(
        self, add_index_map=None, remove_smiles=None, return_index_map=True
    ):
        if self.is_ok == False:
            raise RuntimeError(
                "Molecule not OK, refusing to write PDBQT\n\nLOG:\n%s" % self.log
            )
        if add_index_map is None:
            add_index_map = self.add_index_map
        if remove_smiles is None:
            remove_smiles = self.remove_smiles
        if self.setup is not None:
            pdbqt_string = self._writer.write_string(
                self.setup, add_index_map, remove_smiles
            )
            if return_index_map:
                return pdbqt_string, self._writer._numbering
            return pdbqt_string
        else:
            raise RuntimeError(
                "Cannot generate PDBQT file, the molecule is not prepared."
            )


def pdbqt_to_mol(
    pdbqt_file, is_string=False, exclude_hydrogens=False, return_index_map=False
):
    """
    Convert a PDBQT file or string to an RDKit molecule object and optionally add hydrogen atoms.

    Args:
        pdbqt_file (file/string):  A path to a PDBQT file or a PDBQT string.
        is_string (bool, optional): Specifies whether the input is a string. Defaults to False.
        exclude_hydrogens (bool, optional): If True, hydrogen atoms are excluded from the output molecule. Defaults to False.
        return_index_map (bool, optional): If True, returns a dictionary mapping PDBQT atom indices to RDKit molecule atom indices. Defaults to False.
    Returns:
        Union[Chem.rdchem.Mol, Tuple[Chem.rdchem.Mol, Dict[int, int]]]:
            - If return_index_map is False, returns an RDKit molecule object.
            - If return_index_map is True, returns a tuple containing the RDKit molecule object and a dictionary mapping PDBQT atom indices to RDKit molecule atom indices.
    """
    if is_string:
        pdbqt_str = pdbqt_file
    else:
        with open(pdbqt_file, "r") as f:
            pdbqt_str = f.read()

    pdbqtmol = MyPDBQTMolecule(pdbqt_str)
    rkmol = pdbqtmol.export_rdkit_mol(
        exclude_hydrogens=exclude_hydrogens, return_index_map=return_index_map
    )
    if return_index_map:
        rkmol, pdbqtIdx_to_rkmolIdx = rkmol
        return rkmol, pdbqtIdx_to_rkmolIdx
    return rkmol


def mol_to_sdf(mol, sdf_file):
    with Chem.SDWriter(sdf_file) as writer:
        writer.write(mol)


def mol_to_pdbqt(
    rkmol, pdbqt_file=None, keep_nonpolar_hydrogens=True, return_index_map=True
):
    """
    Convert an RDKit molecule object to a PDBQT string. meeko==0.3.3.

    Args:
        rkmol (Chem.rdchem.Mol): An RDKit molecule object.
        pdbqt_file (str, optional): If provided, the PDBQT content will be saved to this file. Otherwise, this function returns a string. Defaults to None.
        keep_nonpolar_hydrogens (bool, optional): Whether to keep hydrogens. If False, it will "merge hydrogens bound to carbons". Defaults to True.
        return_index_map (bool, optional): Whether to return a dictionary mapping RDKit molecule atom indices to PDBQT atom indices. Defaults to True.

    Returns:
        Union[str, Tuple[str, Dict[int, int]]]:
            - If return_index_map is False, returns a string containing the PDBQT content.
            - If return_index_map is True, returns a tuple containing the PDBQT string and a dictionary mapping RDKit molecule atom indices to PDBQT atom indices.
    """
    preparator = MyMoleculePreparation(
        keep_nonpolar_hydrogens=keep_nonpolar_hydrogens,
        add_index_map=True,
        rigid_macrocycles=True,
    )
    preparator.prepare(rkmol)
    pdbqt_str = preparator.write_pdbqt_string(return_index_map=return_index_map)
    if return_index_map:
        pdbqt_str, index_map = pdbqt_str
        rkmolIdx_to_pdbqtIdx = {k: v - 1 for k, v in index_map.items()}
    if pdbqt_file is not None:
        with open(pdbqt_file, "w") as f:
            f.write(pdbqt_str)
    if return_index_map:
        return pdbqt_str, rkmolIdx_to_pdbqtIdx
    return pdbqt_str
