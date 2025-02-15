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
from itertools import combinations

import parmed as pmd
from parmed.amber import AmberParm

from pxdock.common import get_logger
from pxdock.common.amber_forcefield import (
    AMBERFF_ATOMTYPES,
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

        pdb4amber.main(command[1:])
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
    run_cmd(working_dir, command, "tleap", "", logger)
    parm = pmd.load_file(prmtop, inpcrd)

    return parm


def parm_to_ffdata(parm: AmberParm) -> dict[str, list]:
    """
    Parse required forcefield parameters from an AmberParm object.
    Only restricted types of forcefield terms are allowed.

    Args:
        parm (AmberParm): The input AmberParm object.

    Returns:
        dict: A dictionary containing forcefield data.
        list: The xyz coordinates of the atoms.
    """
    ffdata = {}
    assert not parm.adjusts
    assert not parm.impropers
    assert not parm.rb_torsions
    assert all([not atom.children for atom in parm.atoms])
    assert all([i == atom.idx for i, atom in enumerate(parm.atoms)])
    ffdata["partial_charges"] = [atom.charge for atom in parm.atoms]
    ffdata["atomic_numbers"] = [atom.atomic_number for atom in parm.atoms]
    ffdata["FF_vdW_sigma"] = [float(atom.sigma) for atom in parm.atoms]
    ffdata["FF_vdW_epsilon"] = [float(atom.epsilon) for atom in parm.atoms]
    ffdata["FF_vdW_atomidx"] = [[atom.idx] for atom in parm.atoms]
    ffdata["FF_vdW_paraidx"] = [
        AMBERFF_ATOMTYPES.index(atom.atom_type.name) + 1000 for atom in parm.atoms
    ]
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

    xyz = [[atom.xx, atom.xy, atom.xz] for atom in parm.atoms]
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
    nonbonded12 = set([tuple(sorted([p[0], p[1]])) for p in ffdata["FF_Bonds_atomidx"]])
    nonbonded13 = set(
        [tuple(sorted([p[0], p[2]])) for p in ffdata["FF_Angles_atomidx"]]
    )
    nonbonded14 = (
        set([tuple(sorted([p[0], p[3]])) for p in ffdata["FF_ProperTorsions_atomidx"]])
        - nonbonded12
        - nonbonded13
    )
    nonbondedall = set(
        map(
            tuple,
            map(
                sorted,
                combinations([atomidx[0] for atomidx in ffdata["FF_vdW_atomidx"]], 2),
            ),
        )
    )

    nonbondedall -= nonbonded14 | nonbonded13 | nonbonded12
    nonbondedall, nonbonded14 = list(nonbondedall), list(nonbonded14)
    if sort_idx:
        nonbondedall, nonbonded14 = sorted(nonbondedall), sorted(nonbonded14)

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
        [angle.atom1.idx, angle.atom2.idx, angle.atom3.idx] for angle in parm.angles
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
                [
                    dihedral.atom1.idx,
                    dihedral.atom2.idx,
                    dihedral.atom3.idx,
                    dihedral.atom4.idx,
                ]
            )
            ffdata[f"FF_{type}Torsions_periodicity"].append([dihedral.type.per])
            ffdata[f"FF_{type}Torsions_phase"].append([dihedral.type.phase])
            ffdata[f"FF_{type}Torsions_k"].append([dihedral.type.phi_k])

    return ffdata
