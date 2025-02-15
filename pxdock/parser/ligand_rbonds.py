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

from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts


def find_pattern(mol, smarts):
    p = Chem.MolFromSmarts(smarts)
    return mol.GetSubstructMatches(p)


def get_rdkit_rotatable_bonds(rkmol: Chem.Mol):
    rot_bonds = rkmol.GetSubstructMatches(RotatableBondSmarts)
    # return sorted list of rot_bonds
    return [set(bond) for bond in rot_bonds]
