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

# This file contains code derived from and copied directly from posebusters-0.2.8,
# which is licensed under the BSD 3-Clause License.

# Original Copyright (c) 2023, Martin Buttenschoen
# All rights reserved.

# The full license text is included in the LICENSE file in this project.

import os
from copy import deepcopy
from io import StringIO
from typing import List

import MDAnalysis as mda
import numpy as np
import rdkit.Chem as Chem
import torch
from byteff.mol import Conformer, Molecule
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolAlign import CalcRMS, GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs, RemoveStereochemistry
from tqdm import tqdm

from .logger import get_logger
from .utilities import run_with_timeout, save_sdf_supplier

logger = get_logger(__name__)


def check_rmsd(
    mol_pred: Mol,
    mol_true: Mol,
    rmsd_threshold: float = 2.0,
    heavy_only: bool = True,
    choose_by: str = "rmsd",
    timeout_duration=10,
):
    """
    If the calculation takes more than {timeout_duration} seconds, then break, reduce maxMatches from 1000000 to 2000 and try again.
    """
    configs = {
        "rmsd_threshold": rmsd_threshold,
        "heavy_only": heavy_only,
        "choose_by": choose_by,
        "maxMatches": 1000000,
    }
    try:
        return run_with_timeout(
            _check_rmsd,
            args=(mol_pred, mol_true),
            kwargs=dict(configs),
            timeout_duration=timeout_duration,
        )
    except Exception as e:
        logger.warning(f"[rmsd Error] [first try] {e}")
        logger.warning(f"[rmsd retry] Will set maxMatches to 2000 and retry.")
        configs["maxMatches"] = 2000
        return run_with_timeout(
            _check_rmsd,
            args=(mol_pred, mol_true),
            kwargs=dict(configs),
            timeout_duration=timeout_duration,
        )


def _check_rmsd(
    mol_pred: Mol,
    mol_true: Mol,
    rmsd_threshold: float = 2.0,
    heavy_only: bool = True,
    choose_by: str = "rmsd",
    maxMatches: int = 1000000,
):
    """
    This function is derived from the posebusters.modules.rmsd.check_rmsd, which is licensed under the BSD 3-Clause License.

    Main changes compared to the original code:
    - pass maxMatches argument to CalcRMS and GetBestRMS.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        mol_true: Ground truth molecule (crystal ligand) with at least one conformer. If multiple conformers are
            present, the lowest RMSD will be reported.
        rmsd_threshold: Threshold in angstrom for reporting whether RMSD is within threshold. Defaults to 2.0.
        heavy_only: Whether to only consider heavy atoms for RMSD calculation. Defaults to True.
        choose_by: Metric to choose which mol_true conformation to compare to. Defaults to "rmsd".
        maxMatches: the max number of matches found in a SubstructMatch(). Defaults to 1000000.
    Returns:
        PoseBusters results dictionary.
    """

    assert isinstance(mol_true, Mol), "Ground truth molecule is missing."
    assert isinstance(mol_pred, Mol), "Predicted molecule is missing."
    num_conf = mol_true.GetNumConformers()
    assert num_conf > 0, "Crystal ligand needs at least one conformer."
    assert (
        mol_pred.GetNumConformers() == 1
    ), "Docked ligand should only have one conformer."

    rmsds = [
        robust_rmsd(
            mol_true,
            mol_pred,
            conf_id_probe=i,
            heavy_only=heavy_only,
            maxMatches=maxMatches,
        )
        for i in range(num_conf)
    ]
    kabsch_rmsds = [
        robust_rmsd(
            mol_true,
            mol_pred,
            conf_id_probe=i,
            kabsch=True,
            heavy_only=heavy_only,
            maxMatches=maxMatches,
        )
        for i in range(num_conf)
    ]
    intercentroids = [
        intercentroid(mol_true, mol_pred, conf_id_probe=i, heavy_only=heavy_only)
        for i in range(num_conf)
    ]

    if choose_by == "rmsd":
        i = np.argmin(np.nan_to_num(rmsds, nan=np.inf))
    elif choose_by == "kabsch_rmsd":
        i = np.argmin(np.nan_to_num(kabsch_rmsds, nan=np.inf))
    elif choose_by == "centroid_distance":
        i = np.argmin(np.nan_to_num(intercentroids, nan=np.inf))
    else:
        raise ValueError(
            f"Invalid value {choose_by} for choose_by. Use 'rmsd', 'kabsch_rmsd', 'centroid_distance'."
        )
    rmsd = rmsds[i]
    kabsch_rmsd = kabsch_rmsds[i]
    centroid_dist = intercentroids[i]

    rmsd_within_threshold = rmsd <= rmsd_threshold

    results = {
        "rmsd": rmsd,
        "kabsch_rmsd": kabsch_rmsd,
        "centroid_distance": centroid_dist,
        "rmsd_within_threshold": rmsd_within_threshold,
    }
    return {"results": results}


def robust_rmsd(
    mol_probe: Mol,
    mol_ref: Mol,
    conf_id_probe: int = -1,
    conf_id_ref: int = -1,
    drop_stereo: bool = False,
    heavy_only: bool = True,
    kabsch: bool = False,
    symmetrizeConjugatedTerminalGroups=True,
    **params,
) -> float:
    """
    This function is derived from the posebusters.modules.rmsd.robust_rmsd, which is licensed under the BSD 3-Clause License.

    Main changes compared to the original code:
    - do not modify molecule and do not retry when rdkit RMSD calculation fails.

    RMSD calculation for isomers.
    """
    mol_probe = deepcopy(mol_probe)
    mol_ref = deepcopy(mol_ref)  # copy mols because rdkit RMSD calculation aligns mols

    if drop_stereo:
        RemoveStereochemistry(mol_probe)
        RemoveStereochemistry(mol_ref)

    if heavy_only:
        mol_probe = RemoveHs(mol_probe)
        mol_ref = RemoveHs(mol_ref)

    # combine parameters
    params = dict(
        symmetrizeConjugatedTerminalGroups=symmetrizeConjugatedTerminalGroups,
        kabsch=kabsch,
        **params,
    )

    # calculate RMSD
    rmsd = _call_rdkit_rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    return np.nan


def _call_rdkit_rmsd(
    mol_probe: Mol, mol_ref: Mol, conf_id_probe: int, conf_id_ref: int, **params
):
    """
    Original code from the posebusters, which is licensed under the BSD 3-Clause License.
    """
    try:
        return _rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, **params)
    except RuntimeError:
        pass
    except ValueError:
        pass

    return np.nan


def _rmsd(
    mol_probe: Mol,
    mol_ref: Mol,
    conf_id_probe: int,
    conf_id_ref: int,
    kabsch: bool = False,
    **params,
):
    """
    Original code from the posebusters, which is licensed under the BSD 3-Clause License.
    """
    if kabsch is True:
        return GetBestRMS(
            prbMol=mol_probe,
            refMol=mol_ref,
            prbId=conf_id_probe,
            refId=conf_id_ref,
            **params,
        )
    return CalcRMS(
        prbMol=mol_probe,
        refMol=mol_ref,
        prbId=conf_id_probe,
        refId=conf_id_ref,
        **params,
    )


def intercentroid(
    mol_probe: Mol,
    mol_ref: Mol,
    conf_id_probe: int = -1,
    conf_id_ref: int = -1,
    heavy_only: bool = True,
) -> float:
    """
    Original code from the posebusters, which is licensed under the BSD 3-Clause License.
    Distance between centroids of two molecules.
    """
    if heavy_only:
        mol_probe = RemoveHs(mol_probe)
        mol_ref = RemoveHs(mol_ref)

    centroid_probe = mol_probe.GetConformer(conf_id_probe).GetPositions().mean(axis=0)
    centroid_ref = mol_ref.GetConformer(conf_id_ref).GetPositions().mean(axis=0)

    return float(np.linalg.norm(centroid_probe - centroid_ref))


def get_rmsd_results(mol, true_mol, timeout_duration):
    configs = {
        "rmsd_threshold": 2.0,
        "heavy_only": True,
        "timeout_duration": timeout_duration,
    }
    try:
        res = check_rmsd(mol, true_mol, **configs)
        return res["results"]
    except Exception as e:
        logger.warning(f"[rmsd Error] {e}")
        return {}


def calculate_rmsd(
    predict_sdf, true_sdf, timeout_duration_per_mol=3, progress_bar=False
):
    # read the first molecule from true_sdf
    with Chem.SDMolSupplier(true_sdf) as supplier:
        for mol in supplier:
            true_mol = mol
            break

    results = []
    supp = save_sdf_supplier(
        predict_sdf, skip_failed_mol=True, removeHs=True
    )  # save load, ignore some wrong data, removeHs since we only consider rmsd of heavy atoms
    if progress_bar:
        supp = tqdm(supp)
    for mol in supp:
        props = mol.GetPropsAsDict()
        if "ligand" not in props:
            props["ligand"] = mol.GetProp("_Name")
        props.update(
            get_rmsd_results(mol, true_mol, timeout_duration=timeout_duration_per_mol)
        )
        results.append(props)

    return results


class MyMolecule(Molecule):
    def to_sdf(
        self, filename: str, conf_id: int = None, append=False, write_all_props=True
    ) -> int:
        """
        Molecule: only write keys satisfying k.startswith('prop_') to sdf properties.
        MyMolecule: write all properties in confdata to sdf properties.
        """
        rkmol = self.to_rkmol(conf_id=conf_id)
        if conf_id is not None:
            assert -self.nconfs <= conf_id < self.nconfs
            conf_list = [self.conformers[conf_id]]
        else:
            conf_list = self.conformers
        if not filename.endswith("sdf"):
            filename = filename + ".sdf"

        mode = "a" if append else "w"
        with open(filename, mode) as fout:
            with Chem.SDWriter(fout) as w:
                for i, conf in enumerate(conf_list):
                    confdata = conf.confdata
                    # add conf prop to rkmol for writing
                    for k, v in confdata.items():
                        if k == "coords":
                            continue
                        if write_all_props or k.startswith("prop_"):
                            rkmol.SetProp(k, str(v))
                        if k == "ligand":
                            rkmol.SetProp("_Name", str(v))
                    w.write(rkmol, confId=i)
                    # remove conf prop for next
                    for k, v in confdata.items():
                        if k.startswith("prop_"):
                            rkmol.ClearProp(k)

        return len(conf_list)


class PoseWriter:
    def __init__(
        self,
        ligand_smiles_string=None,
        receptor_gro_string=None,
    ):
        if ligand_smiles_string is not None:
            # create ligand byteff mol object
            self.ligand_smiles_string = ligand_smiles_string
            self.ligand = MyMolecule.from_mapped_smiles(ligand_smiles_string)
        else:
            self.ligand = None
        if receptor_gro_string is not None:
            # create receptor mda universe object
            stringstream = StringIO(receptor_gro_string)
            self.receptor = mda.Universe(stringstream, format="GRO")
        else:
            self.receptor = None

    def clear_ligand_pose(self):
        for i in range(self.ligand.nconfs):
            self.ligand.remove_conformer(0)

    def set_ligand_pose(
        self,
        xyz: [torch.Tensor, np.ndarray, List],  # [nposes, natoms, 3]
        confdata: list = None,  # [nposes]
    ):
        self.clear_ligand_pose()
        if confdata is None:
            confdata = [None] * len(xyz)
        assert len(xyz) == len(confdata)
        if isinstance(xyz, list):
            xyz = np.array(xyz)
        elif isinstance(xyz, torch.Tensor):
            xyz = xyz.detach().numpy()
        for xyz_i, confdata_i in zip(xyz, confdata):
            conformer = Conformer(xyz_i, self.ligand.atomic_numbers, confdata_i)
            self.ligand.append_conformers(conformer)

    def set_receptor_pose(
        self,
        xyz: [torch.Tensor, np.ndarray, List],  # [nposes, natoms, 3]
        index: list = None,
    ):
        assert (
            self.receptor is not None
        ), f"receptor is None since GRO_string is not given."
        if isinstance(xyz, list):
            xyz = np.array(xyz)
        elif isinstance(xyz, torch.Tensor):
            xyz = xyz.detach().numpy()
        if index is not None:
            xyz_all = self.receptor.atoms.positions
            assert len(xyz_all) > 0
            xyz_all = xyz_all.reshape(-1, *xyz_all.shape).repeat(len(xyz), axis=0)
            xyz_all[:, index, :] = xyz
            xyz = xyz_all
        self.receptor.load_new(xyz)

    def write_ligand(self, output_dir, fname, append=False):
        if not fname.endswith(".sdf"):
            fname = fname + ".sdf"
        os.makedirs(output_dir, exist_ok=True)
        fpath = os.path.join(output_dir, fname)
        self.ligand.to_sdf(fpath, append=append)

    def write_receptor(self, output_dir, fname):
        if not fname.endswith(".pdb"):
            fname = fname + ".pdb"
        os.makedirs(output_dir, exist_ok=True)
        fpath = os.path.join(output_dir, fname)
        with mda.Writer(
            fpath, multiframe=True, bonds=None, n_atoms=self.receptor.atoms.n_atoms
        ) as PDB:
            for ts in self.receptor.trajectory:
                PDB.write(self.receptor.atoms)

    def write(self, output_dir, fname):
        self.write_ligand(output_dir, f"ligand_{fname}")
        self.write_receptor(output_dir, f"receptor_{fname}")


def write_ligand_to_sdf(
    mapped_smiles: str,
    xyz: [torch.Tensor, np.ndarray, List],
    output_fpath: str,
    confdata: List[dict] = None,
    append=False,
):
    pw = PoseWriter(
        ligand_smiles_string=mapped_smiles,
        receptor_gro_string=None,
    )
    if not isinstance(xyz, (torch.Tensor, np.ndarray)):
        xyz = np.array(xyz)
    pw.set_ligand_pose(xyz, confdata)
    output_dir = os.path.dirname(os.path.abspath(output_fpath))
    os.makedirs(output_dir, exist_ok=True)
    pw.write_ligand(output_dir, os.path.basename(output_fpath), append=append)

    return


def write_receptor_to_pdb(
    gro_string: str,
    xyz: [torch.Tensor, np.ndarray, List],
    index: [torch.Tensor, np.ndarray, List],
    output_fpath: str,
):
    pw = PoseWriter(ligand_smiles_string=None, receptor_gro_string=gro_string)
    if xyz is not None:
        pw.set_receptor_pose(xyz, index=index)
    output_dir = os.path.dirname(output_fpath)
    os.makedirs(output_dir, exist_ok=True)
    pw.write_receptor(output_dir, os.path.basename(output_fpath))

    return
