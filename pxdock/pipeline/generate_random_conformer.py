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

import argparse
import copy
import os

from byteff.mol import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def clear_all_properties(mol):
    prop_names = list(mol.GetPropNames())
    for prop_name in prop_names:
        mol.ClearProp(prop_name)
    return mol


def get_heavy_atom_indices(mol):
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]


def aligned_rmsd_of_heavy_atoms(mol1, mol2):
    mol1 = copy.deepcopy(mol1)
    mol2 = copy.deepcopy(mol2)
    heavy_atom_indices = get_heavy_atom_indices(mol1)
    # Calculate the aligned RMSD between the two molecules
    atom_map = [(idx, idx) for idx in heavy_atom_indices]
    rmsd = AllChem.AlignMol(mol1, mol2, atomMap=atom_map)
    return rmsd


def process_one_mol(
    mol, rename_ligand=None, seed=42, max_iteration_rmsd_check=3, verbose=False
):
    """
    Args:
        max_iteration_rmsd_check (int):
            Limit the maximum number of attempts to generate a conformer. Check if the aligned RMSD meets the requirements.
            If the generated conformer does not meet the requirements, regenerate it.
            If <= 1, it is equivalent to not checking and using the first generated conformer. Defaults to 3.
    """

    seed = int(seed)
    name = mol.GetProp("_Name")
    if rename_ligand:
        name = rename_ligand
    mapped_smiles = Molecule(mol).get_mapped_smiles()

    max_iteration_rmsd_check = max(1, max_iteration_rmsd_check)

    while max_iteration_rmsd_check > 0:
        # Generate a new conformer using the mapped SMILES string, instead of using the conformer in the SDF file
        molecule = Molecule.from_mapped_smiles(mapped_smiles, nconfs=1, randomSeed=seed)
        new_mol = molecule.to_rkmol()
        # Check if the aligned RMSD meets the requirements
        n_rotb = Descriptors.NumRotatableBonds(mol)
        rmsd = aligned_rmsd_of_heavy_atoms(mol, new_mol)
        if (0.3 * n_rotb <= rmsd) and rmsd <= 2.0 * (n_rotb + 1):
            break
        max_iteration_rmsd_check -= 1
        seed += 1

    if verbose:
        print(f"Aligned RMSD: {rmsd}")
        print(f"Number of rotatable bonds: {n_rotb}")

    # add properties: _Name / ori_smiles / stereo_smiles / mapped_smiles
    new_mol = clear_all_properties(new_mol)
    for prop_name in mol.GetPropNames():
        new_mol.SetProp(prop_name, mol.GetProp(prop_name))

    stereo_smiles = molecule.get_smiles()
    new_mol.SetProp("_Name", name)
    new_mol.SetProp("smiles", stereo_smiles)
    return new_mol


def process_one_sdf(
    input_filename, output_filename, rename_ligand_with_filename=False, **kwargs
):
    """
    input_filename: .sdf
    output_filename: .sdf
    rename_ligand_with_filename:
        If True, use the input file name to replace the ligand names in the SDF file.
        This is useful when some reference ligands do not have unique names. Defaults to False.
    """

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    writer = Chem.SDWriter(output_filename)

    for i, mol in enumerate(Chem.SDMolSupplier(input_filename)):
        rename_ligand = None
        if rename_ligand_with_filename:
            rename_ligand = os.path.basename(input_filename)[: -len(".sdf")]
            if i > 0:
                rename_ligand = rename_ligand + f"_{i}"
        new_mol = process_one_mol(mol, rename_ligand=rename_ligand, **kwargs)
        writer.write(new_mol)

    writer.close()


if __name__ == "__main__":
    _parser = argparse.ArgumentParser()
    _parser.add_argument("--input_fpath", type=str, required=True)
    _parser.add_argument("--output_fpath", type=str, required=True)
    _parser.add_argument("--max_iteration_rmsd_check", type=int, default=3)
    _parser.add_argument("--rename_ligand_with_filename", action="store_true")
    _parser.add_argument("--verbose", action="store_true")
    args = _parser.parse_args()
    print(args)

    process_one_sdf(
        input_filename=args.input_fpath,
        output_filename=args.output_fpath,
        rename_ligand_with_filename=args.rename_ligand_with_filename,
        max_iteration_rmsd_check=args.max_iteration_rmsd_check,
        verbose=args.verbose,
    )
