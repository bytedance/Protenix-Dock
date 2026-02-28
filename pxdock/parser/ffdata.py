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
import parmed as pmd
from parmed.amber import AmberParm

from pxdock.common import get_logger
from pxdock.common.amber_forcefield import (
    AMBERFF_ATOMTYPES,
    AMBERFF_ATOMTYPES2idx,
    FORCEFIELD,
    FORCEFIELD_WATER_ION,
)
from pxdock.common.utilities import run_cmd

logger = get_logger(__name__)
ALL_FF_KEYS = [
    "FF_Angles_atomidx",
    "FF_Bonds_atomidx",
    "FF_ProperTorsions_atomidx",
    "FF_ImproperTorsions_atomidx",
    "FF_NonbondedAll_atomidx",
    "FF_Nonbonded14_atomidx",
    "FF_Bonds_k",
    "FF_Bonds_length",
    "FF_Angles_k",
    "FF_Angles_angle",
    "FF_ProperTorsions_periodicity",
    "FF_ProperTorsions_phase",
    "FF_ProperTorsions_k",
    "FF_ImproperTorsions_periodicity",
    "FF_ImproperTorsions_phase",
    "FF_ImproperTorsions_k",
    "FF_vdW_atomidx",
    "FF_vdW_sigma",
    "FF_vdW_epsilon",
    "FF_vdW_paraidx",
    "partial_charges",
    "atomic_numbers",
]


def check_vdw_paraidx(
    vdw_sigma,
    vdw_epsilon,
    vdw_paraidx,
):
    """
    Check if the van der Waals parameter index uniquely corresponds to a set of van der Waals parameters.

    Args:
        vdw_sigma (list): List of van der Waals radii.
        vdw_epsilon (list): List of van der Waals energy parameters.
        vdw_paraidx (list): List of van der Waals parameter indices.

    Raises:
        AssertionError: If the same parameter index corresponds to different sets of van der Waals parameters.
    """
    paraidx_to_para = {}
    for sigma, epsilon, paraidx in zip(vdw_sigma, vdw_epsilon, vdw_paraidx):
        if paraidx not in paraidx_to_para:
            paraidx_to_para[paraidx] = (sigma, epsilon)
        else:
            assert (
                abs(paraidx_to_para[paraidx][0] - sigma) < 1e-3
            ), f"paraidx {paraidx} has different vdw parameters"
            assert (
                abs(paraidx_to_para[paraidx][1] - epsilon) < 1e-3
            ), f"paraidx {paraidx} has different vdw parameters"


def pdb_to_parm(working_dir, pdbfile, forcefield="FF14SB"):

    # step 1: remove Hs using pdb4amber
    amber_pdb = os.path.join(working_dir, "step1_amber.pdb")
    command = [
        "pdb4amber",
        "-i",
        pdbfile,
        "-o",
        amber_pdb,
        "--nohyd",  # remove all hydrogens
        "--dry",  # remove all water molecules
        "--no-conect",
    ]

    # step 2: generate force field in Amber format
    ret_code = subprocess.run("which pdb4amber", shell=True, capture_output=True)
    if ret_code.returncode == 0:
        run_cmd(working_dir, command, "pdb4amber", "", logger)
    else:
        import pdb4amber

        pdb4amber.run(
            arg_pdbin=pdbfile,
            arg_pdbout=amber_pdb,
            arg_nohyd=True,
            arg_dry=True,
            arg_conect=False,
        )
    prmtop = os.path.join(working_dir, "step2_protein_amber.prmtop")
    inpcrd = os.path.join(working_dir, "step2_protein_amber.inpcrd")
    tleap_in = [
        f"source {FORCEFIELD[forcefield.upper()]}",
        f"source {FORCEFIELD_WATER_ION[forcefield.upper()]}",
        f"loadamberparams frcmod.hemall",
        f"loadamberprep heme_all.in",
        f"pro = loadpdb {amber_pdb}",
        f"saveamberparm pro {prmtop} {inpcrd}",
        "quit",
    ]
    tleap_input_file = os.path.join(working_dir, "step2_tleap.in")
    with open(tleap_input_file, "w", encoding="utf-8") as fout:
        fout.write("\n".join(tleap_in))
    amberff_dir = os.path.join(os.path.dirname(__file__), "../data/amberff")
    command = [
        "tleap",
        "-f",
        tleap_input_file,
        "-I",
        amberff_dir,
    ]
    ret = run_cmd(working_dir, command, "tleap", "", logger)
    if ret != 0:
        raise RuntimeError(
            f"tleap exited with error code {ret}. "
            f"Check logs at {os.path.join(working_dir, 'tleap_stdout.log')} "
            f"and {os.path.join(working_dir, 'tleap_stderr.log')}"
        )
    if not os.path.exists(prmtop) or not os.path.exists(inpcrd):
        raise FileNotFoundError(
            f"tleap did not produce expected output files. "
            f"This usually means the PDB structure has residues or atoms "
            f"incompatible with the AMBER force field. "
            f"Check {os.path.join(working_dir, 'tleap_stdout.log')} for details."
        )
    parm = pmd.load_file(prmtop, inpcrd)

    return parm


def parm_to_ffdata(parm: AmberParm) -> tuple[dict[str, list], list[tuple[float]]]:
    """
    Parse required forcefield parameters from an AmberParm object.
    Only restricted types of forcefield terms are allowed.

    Args:
        parm (AmberParm): The input AmberParm object.

    Returns:
        dict: A dictionary containing forcefield data.
        list: The xyz coordinates of the atoms. (N, 3)
    """
    ffdata = {}
    assert not parm.adjusts
    assert not parm.impropers
    assert not parm.rb_torsions

    partial_charges = []
    atomic_numbers = []
    vdw_sigma = []
    vdw_epsilon = []
    vdw_atomidx = []
    vdw_paraidx = []
    for i, atom in enumerate(parm.atoms):
        assert not atom.children
        assert i == atom.idx
        partial_charges.append(atom.charge)
        atomic_numbers.append(atom.atomic_number)
        vdw_sigma.append(float(atom.sigma))
        vdw_epsilon.append(float(atom.epsilon))
        vdw_atomidx.append([atom.idx])
        vdw_paraidx.append(AMBERFF_ATOMTYPES2idx[atom.atom_type.name] + 1000)
    ffdata["partial_charges"] = partial_charges
    ffdata["atomic_numbers"] = atomic_numbers
    ffdata["FF_vdW_sigma"] = vdw_sigma
    ffdata["FF_vdW_epsilon"] = vdw_epsilon
    ffdata["FF_vdW_atomidx"] = vdw_atomidx
    ffdata["FF_vdW_paraidx"] = vdw_paraidx
    check_vdw_paraidx(
        ffdata["FF_vdW_sigma"], ffdata["FF_vdW_epsilon"], ffdata["FF_vdW_paraidx"]
    )
    ffdata.update(parm_to_bond_params(parm))
    ffdata.update(parm_to_angle_params(parm))
    ffdata.update(parm_to_torsion_params(parm, type="Proper"))
    ffdata.update(parm_to_torsion_params(parm, type="Improper"))
    ffdata.update(calc_nonbonded_atomidx(ffdata))
    for key in ALL_FF_KEYS:
        assert key in ffdata, key

    xyz = [(atom.xx, atom.xy, atom.xz) for atom in parm.atoms]
    return ffdata, xyz


def calc_nonbonded_atomidx(
    ffdata: dict[str, list], sort_idx: bool = True
) -> dict[str, list]:
    """
    Calculate the non-bonded interactions atom indices from the bond, angle, proper torsion atom indices.

    Args:
        ffdata (dict): A dictionary containing forcefield data.
        sort_idx (bool): Whether to sort the indices of the non-bonded interactions.

    Returns:
        dict: A dictionary containing the indices of the non-bonded interactions.
    """
    for key in [
        "FF_vdW_atomidx",
        "FF_Bonds_atomidx",
        "FF_Angles_atomidx",
        "FF_ProperTorsions_atomidx",
    ]:
        assert key in ffdata

    def _encode_pairs(pairs):
        """Encode sorted (i, j) pairs as single int64 for fast set ops."""
        arr = np.array(pairs, dtype=np.int64)
        if arr.size == 0:
            return np.empty(0, dtype=np.int64)
        a, b = arr[:, 0], arr[:, 1]
        lo = np.minimum(a, b)
        hi = np.maximum(a, b)
        return lo * n_atoms + hi

    n_atoms = len(ffdata["FF_vdW_atomidx"])

    bond_pairs = [(p[0], p[1]) for p in ffdata["FF_Bonds_atomidx"]]
    angle_pairs = [(p[0], p[2]) for p in ffdata["FF_Angles_atomidx"]]
    torsion_pairs = [(p[0], p[3]) for p in ffdata["FF_ProperTorsions_atomidx"]]

    enc12 = _encode_pairs(bond_pairs)
    enc13 = _encode_pairs(angle_pairs)
    enc_torsion = _encode_pairs(torsion_pairs)

    # 1-4 pairs = torsion pairs minus 1-2 and 1-3
    exclude_from_14 = np.union1d(enc12, enc13)
    enc14 = np.setdiff1d(np.unique(enc_torsion), exclude_from_14)

    # all pairs via triu_indices minus 1-2, 1-3, 1-4
    idx_i, idx_j = np.triu_indices(n_atoms, k=1)
    enc_all = idx_i.astype(np.int64) * n_atoms + idx_j.astype(np.int64)
    exclude_all = np.union1d(exclude_from_14, enc14)
    enc_nonbonded = np.setdiff1d(enc_all, exclude_all)

    # decode back to tuple lists
    def _decode(enc):
        if enc.size == 0:
            return []
        lo = enc // n_atoms
        hi = enc % n_atoms
        return list(zip(lo.tolist(), hi.tolist()))

    nonbondedall = _decode(enc_nonbonded)
    nonbonded14 = _decode(enc14)
    if sort_idx:
        nonbondedall.sort()
        nonbonded14.sort()

    return {
        "FF_NonbondedAll_atomidx": nonbondedall,
        "FF_Nonbonded14_atomidx": nonbonded14,
    }


def parm_to_bond_params(parm: AmberParm) -> dict[str, list]:
    """
    Parse harmonic bond forcefield parameters from an AmberParm object.

    Args:
        parm (AmberParm): The input AmberParm object.

    Returns:
        dict: A dictionary containing forcefield data.
    """
    if not parm.bonds:
        return {
            "FF_Bonds_atomidx": [],
            "FF_Bonds_k": [],
            "FF_Bonds_length": [],
        }
    ffdata = {}
    assert all([bond.funct == 1 for bond in parm.bonds])
    parm.bonds.sort(key=lambda x: x.atom2.idx)
    parm.bonds.sort(key=lambda x: x.atom1.idx)

    ffdata["FF_Bonds_atomidx"] = [
        [bond.atom1.idx, bond.atom2.idx] for bond in parm.bonds
    ]
    ffdata["FF_Bonds_k"] = [bond.type.k * 2.0 for bond in parm.bonds]
    ffdata["FF_Bonds_length"] = [bond.type.req for bond in parm.bonds]
    return ffdata


def parm_to_angle_params(parm: AmberParm) -> dict[str, list]:
    """
    Parse harmonic angle forcefield parameters from an AmberParm object.

    Args:
        parm (AmberParm): The input AmberParm object.

    Returns:
        dict: A dictionary containing forcefield data.
    """
    if not parm.angles:
        return {
            "FF_Angles_atomidx": [],
            "FF_Angles_k": [],
            "FF_Angles_angle": [],
        }
    ffdata = {}
    assert all([angle.funct == 1 for angle in parm.angles])
    parm.angles.sort(key=lambda x: x.atom3.idx)
    parm.angles.sort(key=lambda x: x.atom2.idx)
    parm.angles.sort(key=lambda x: x.atom1.idx)

    ffdata["FF_Angles_atomidx"] = [
        (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx) for angle in parm.angles
    ]
    ffdata["FF_Angles_k"] = [angle.type.k * 2.0 for angle in parm.angles]
    ffdata["FF_Angles_angle"] = [angle.type.theteq for angle in parm.angles]
    return ffdata


def parm_to_torsion_params(parm: AmberParm, type: str = "Proper") -> dict[str, list]:
    """
    Parse periodic torsion forcefield parameters from an AmberParm object.

    Args:
        parm (AmberParm): The input AmberParm object.
        type (str): The type of torsion to parse. Can be "Proper" or "Improper".

    Returns:
        dict: A dictionary containing forcefield data.
    """
    assert type in ["Proper", "Improper"]
    ffdata = {
        f"FF_{type}Torsions_atomidx": [],
        f"FF_{type}Torsions_periodicity": [],
        f"FF_{type}Torsions_phase": [],
        f"FF_{type}Torsions_k": [],
    }
    if not parm.dihedrals:
        return ffdata
    parm.dihedrals.sort(key=lambda x: x.atom4.idx)
    parm.dihedrals.sort(key=lambda x: x.atom3.idx)
    parm.dihedrals.sort(key=lambda x: x.atom2.idx)
    parm.dihedrals.sort(key=lambda x: x.atom1.idx)
    prev_idx = None
    for dihedral in parm.dihedrals:
        if type == "Proper" and dihedral.improper:
            continue
        if type == "Improper" and not dihedral.improper:
            continue
        if dihedral.improper:
            assert dihedral.funct == 4
        else:
            assert dihedral.funct == 1
            assert abs(dihedral.type.scee - 1.2) < 1e-8, dihedral
            assert abs(dihedral.type.scnb - 2.0) < 1e-8, dihedral
        if prev_idx is not None and prev_idx == (
            dihedral.atom1.idx,
            dihedral.atom2.idx,
            dihedral.atom3.idx,
            dihedral.atom4.idx,
        ):
            ffdata[f"FF_{type}Torsions_periodicity"][-1].append(dihedral.type.per)
            ffdata[f"FF_{type}Torsions_phase"][-1].append(dihedral.type.phase)
            ffdata[f"FF_{type}Torsions_k"][-1].append(dihedral.type.phi_k)
        else:
            prev_idx = (
                dihedral.atom1.idx,
                dihedral.atom2.idx,
                dihedral.atom3.idx,
                dihedral.atom4.idx,
            )
            ffdata[f"FF_{type}Torsions_atomidx"].append(
                (
                    dihedral.atom1.idx,
                    dihedral.atom2.idx,
                    dihedral.atom3.idx,
                    dihedral.atom4.idx,
                )
            )
            ffdata[f"FF_{type}Torsions_periodicity"].append([dihedral.type.per])
            ffdata[f"FF_{type}Torsions_phase"].append([dihedral.type.phase])
            ffdata[f"FF_{type}Torsions_k"].append([dihedral.type.phi_k])

    return ffdata
