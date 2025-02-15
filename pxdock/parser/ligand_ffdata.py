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

########################################################################
#
# This module collects the function of generating force field parameters
# for target systems. Related tools such as antechamber might need to
# be installed priorly.
#
########################################################################

import os
import shutil
import traceback
import uuid
from typing import Dict, Iterable, List

import parmed as pmd
from byteff.mol import Molecule

from pxdock.common import get_logger, kWorkDir, my_random_string
from pxdock.common.utilities import convert_value_to_list, run_command_and_check
from pxdock.parser.ffdata import parm_to_ffdata

logger = get_logger(__name__)

DEFAULT_FORCEFIELD = "gaff2"
DEFAULT_CHARGE_MODEL = "am1bcc"

TLEAP_TEMP = """
verbosity 1
source leaprc.gaff2
mods = loadamberparams {mol_name}_gaff2.frcmod
ligand = loadmol2 MOL_user_gaff2.mol2
check ligand
saveamberparm ligand {mol_name}_gaff2.prmtop {mol_name}_gaff2.inpcrd
saveoff ligand {mol_name}_gaff2.lib
quit
"""

def call_antechamber(
    sdf_file: str,
    natoms: int,
    tot_formal_charge: int,
    *,
    charge_model: str,
    sqm_opt: bool,
):
    """raise exception if any error is detected"""

    env = {"AMBER_BCCDAT": f"BCCPARM.DAT.{charge_model}"}

    assert os.path.exists(sdf_file)

    sdf_file = os.path.abspath(sdf_file)

    if sqm_opt:
        sqm_params = [
            "-ek",
            "\"qm_theory='AM1', scfconv=1.d-10, grms_tol=0.0005, ndiis_attempts=700\" ",
        ]
    else:
        sqm_params = [
            "-ek",
            "\"qm_theory='AM1', scfconv=1.d-10, maxcyc=0, grms_tol=0.0005, ndiis_attempts=700\" ",
        ]

    # run antechamber
    cmd = " ".join(
        [
            "antechamber",
            "-i",
            f"{sdf_file}",
            "-fi",
            "sdf",
            "-o",
            "charged.mol2",
            "-fo",
            "mol2",
            "-pf",
            "yes",
            "-dr",
            "n",
            "-c",
            "bcc",
            "-nc",
            str(tot_formal_charge),
        ]
        + sqm_params
    )
    run_command_and_check(cmd, env=env)

    # Write out just charges
    cmd = "antechamber -dr n -i charged.mol2 -fi mol2 -o charges2.mol2 -fo mol2 -c wc -cf charges.txt -pf yes"
    run_command_and_check(cmd, env=env)

    # Save .mol2 file
    cmd = f"antechamber -dr no -i charges2.mol2 -fi mol2 -o MOL_user_gaff2.mol2 -fo mol2 -nc {tot_formal_charge} -m 1 -s 2 -df 2 -at gaff2 -pf n"
    run_command_and_check(cmd, env=env)
    # read charges.txt to get partial charges
    with open("charges.txt", "r") as infile:
        contents = infile.read()
    text_charges = contents.split()
    partial_charges = [float(c) for c in text_charges]

    assert len(partial_charges) == natoms
    return partial_charges


def assign_partial_charges(
    mol: Molecule,
    *,
    charge_model: str = "am1bcc",  # or path to json file
    sqm_opt: bool = True,
    conf_id: int = 0,
) -> List:
    assert isinstance(mol, Molecule)
    assert len(mol.conformers) > conf_id

    logger.debug(
        f"smiles {mol.get_smiles(isomeric=True)} formal_charges {mol.formal_charges}"
    )

    total_formal_charge = int(sum(mol.formal_charges))
    error_message = None

    try:
        mol.to_sdf("mol.sdf", conf_id=conf_id)
        partial_charges = call_antechamber(
            "mol.sdf",
            mol.natoms,
            total_formal_charge,
            charge_model=charge_model,
            sqm_opt=sqm_opt,
        )
    except Exception:
        error_message = traceback.format_exc()
        logger.warning("%s", error_message)

    if error_message is not None:
        raise RuntimeError(f"{error_message}")

    return partial_charges


def generate_parm(
    mol: Molecule,
    save_path: str,
    *,
    partial_charges_key: str = "partial_charges",  # use this in mol.conformers[conf_id].confdata
    forcefield: str = DEFAULT_FORCEFIELD,
    mol_name: str = "MOL",
    # params shared by forcefield and assign_partial_charges
    conf_id: int = 0,
    verbose: bool = False,
    # params passed to assign_partial_charges
    charge_model: str = DEFAULT_CHARGE_MODEL,
    sqm_opt: bool = True,
) -> Dict:
    """
    Generate parm file with the specified forcefield and charge model.

    Args:
        mol: the molecule to be converted to a inpcrd and prmtop files.
        save_path: where to save the parm file.
        forcefield: only support gaff2.
        charge_model: am1bcc or abcg2
        mol_name: moleculetype in parm file.
        sqm_opt: whether to use am1 optimizatoin in charge calculation.
        conf_id: which conformer to to use in charge calculation.
    Return:
        Dict[str, ...] useful results
    """
    assert forcefield == "gaff2", f"Forcefield {forcefield} not found!"

    # sanitize supported charge model
    assert charge_model in [
        "am1bcc",
        "abcg2",
    ], f"Charge model {charge_model} not implemented!"

    # sanitize other parameters
    if not os.path.exists(save_path):
        os.makedirs(save_path, exist_ok=True)

    # convert all path to abs path
    save_path = os.path.abspath(save_path)

    result = dict()

    # conformation generation
    if mol.nconfs == 0:
        mol = Molecule.from_mapped_smiles(
            mol.get_mapped_smiles(isomeric=True), nconfs=1
        )
        conf_id = 0
        logger.info("Created conformer to generate partial charges")
    else:
        assert -conf_id <= conf_id < mol.nconfs

    result["mapped_smiles_nonisomeric"] = mol.get_mapped_smiles(isomeric=False)
    result["mapped_smiles_isomeric"] = mol.get_mapped_smiles(isomeric=True)

    # stage 1, assign partial charges
    logger.info(f"Using conf {conf_id} to generate partial charges")
    pc = assign_partial_charges(
        mol,
        charge_model=charge_model,
        conf_id=conf_id,
        sqm_opt=sqm_opt,
    )

    result[partial_charges_key] = pc

    # .sdf to .mol2
    mol2_path = "MOL_user_gaff2.mol2"

    # .mol2 to .itp .inpcrd .prmtop ... files
    total_formal_charges = sum(mol.formal_charges)
    # generate .frcmod file
    run_command_and_check(
        f"parmchk2 -i {mol2_path} -f mol2 -o {mol_name}_gaff2.frcmod -s 2",
        allow_error=False,
        separate_stderr=True,
    )
    # write tleap.in file and generate .inpcrd and .prmtop files
    with open("tleap.in", "w") as wf:
        wf.write(TLEAP_TEMP.format(mol_name=mol_name))

    run_command_and_check(
        "tleap -f tleap.in", allow_error=False, separate_stderr=True
    )
    shutil.copy(f"{mol_name}_gaff2.inpcrd", save_path)
    shutil.copy(f"{mol_name}_gaff2.prmtop", save_path)
    result["inpcrd_file"] = os.path.join(save_path, f"{mol_name}_gaff2.inpcrd")
    result["prmtop_file"] = os.path.join(save_path, f"{mol_name}_gaff2.prmtop")
    if verbose:
        logger.info(f"Created {save_path} by {forcefield} successfully!")

    return result


def modify_molecule_weights(m):
    """
    Cap atomic masses such that all atoms with mass greater than chlorine
    are assigned the mass of chlorine.
    """
    from rdkit import Chem

    m = Chem.Mol(m)
    total_cut = 0.0
    for atom in m.GetAtoms():
        if atom.GetMass() > 35.45:  # 35.45 is the mass of chlorine atom
            total_cut += 35.45 - atom.GetMass()
    return total_cut


def calculate_molecular_weight(mapped_smiles_lists):
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    if isinstance(mapped_smiles_lists, str):
        mapped_smiles_lists = [mapped_smiles_lists]
    mols = [Chem.MolFromSmiles(smiles) for smiles in mapped_smiles_lists]
    total_cutoff = [modify_molecule_weights(mol) for mol in mols]
    molecular_weights = [
        Descriptors.MolWt(mol) + cut for (mol, cut) in zip(mols, total_cutoff)
    ]
    return molecular_weights


def _fix_rounding_error(partial_charges: Iterable, total_formal_charge: int) -> List:
    """normalize the partial charges so the sum equals the tot_formal_charges"""
    total_formal_charge = int(sum(total_formal_charge))

    deviate = sum(partial_charges) - total_formal_charge  # on order of 1e-3
    # step 1, distribute error to all atoms
    shift = deviate / len(partial_charges)
    partial_charges = [c - shift for c in partial_charges]
    assert (
        abs(total_formal_charge - sum(partial_charges)) < 1e-6
    ), "sum of partial charges must exactly match total formal charge"

    return partial_charges


def get_param_values(ffdata, ffparams):
    if not isinstance(ffdata, dict):
        ffdata = ffdata.as_dict()

    all_fields = set([x.split("_")[1] for x in ffdata.keys()])

    for field in all_fields:
        if field == "vdW":
            sigma = []
            epsilon = []
            for i, paraidx in enumerate(ffdata[f"FF_{field}_paraidx"]):
                p0 = ffparams.get_param(field, paraidx)
                sigma.append(p0.sigma)
                epsilon.append(p0.epsilon)
            ffdata["FF_vdW_sigma"] = sigma
            ffdata["FF_vdW_epsilon"] = epsilon

        elif field == "Bonds":
            field = field[:-1]
            ffdata["FF_Bonds_k"] = []
            ffdata["FF_Bonds_length"] = []
            for i, paraidx in enumerate(ffdata[f"FF_{field}s_paraidx"]):
                p0 = ffparams.get_param(field, paraidx)
                ffdata["FF_Bonds_k"].append(p0.k)
                ffdata["FF_Bonds_length"].append(p0.length)

        elif field == "Angles":
            field = field[:-1]
            ffdata["FF_Angles_k"] = []
            ffdata["FF_Angles_angle"] = []
            for i, paraidx in enumerate(ffdata[f"FF_{field}s_paraidx"]):
                p0 = ffparams.get_param(field, paraidx)
                ffdata["FF_Angles_k"].append(p0.k)
                ffdata["FF_Angles_angle"].append(p0.angle)

        elif field in ["ProperTorsions", "ImproperTorsions"]:
            field = field[:-1]
            ffdata[f"FF_{field}s_k"] = []
            ffdata[f"FF_{field}s_phase"] = []
            ffdata[f"FF_{field}s_periodicity"] = []

            for i, paraidx in enumerate(ffdata[f"FF_{field}s_paraidx"]):
                p0 = ffparams.get_param(field, paraidx)
                k = []
                phase = []
                periodicity = []
                for term in p0.expansion:
                    k.append(term.k)
                    phase.append(term.phase)
                    periodicity.append(term.periodicity)
                ffdata[f"FF_{field}s_k"].append(k)
                ffdata[f"FF_{field}s_phase"].append(phase)
                ffdata[f"FF_{field}s_periodicity"].append(periodicity)

    return ffdata


def get_ffdata(
    molecule: Molecule,
    forcefield: str = "gaff2",
    conf_id: int = 0,
    charge_model: str = "am1bcc",
    sqm_opt: bool = False,
    parm_save_path: str = None,
    verbose: bool = False,
):
    """
    Generate force field data (ffdata) for a given molecule using the specified force field and charge model.

    Args:
        molecule: the molecule to be converted to a ffdata
        forcefield: The force field to be used for generating the parameters. currently only support gaff2.
        parm_save_path: where to save the ligand parm files, .prmtop and .inpcrd files
        - if None, the parm file will be removed.
        charge_model: am1bcc or abcg2 or path to bcc json file
        sqm_opt: whether to use am1 optimizatoin in charge calculation.
    Return:
        Dict[str, ...] useful results
    """

    # get forcefield params
    save_path = parm_save_path or os.path.join(kWorkDir, my_random_string(10))
    working_directory = os.path.join(kWorkDir, my_random_string(10))
    os.makedirs(working_directory, exist_ok=True)
    old_directory = os.getcwd()
    os.chdir(working_directory)

    results = generate_parm(
        molecule,
        save_path=save_path,
        forcefield=forcefield,
        conf_id=conf_id,
        verbose=verbose,
        charge_model=charge_model,
        sqm_opt=sqm_opt,
    )

    os.chdir(old_directory)
    # genetate_parm generate standard .frcmod and .prmtop files and partial charges, fix the round error now
    round_partial_charges = _fix_rounding_error(
        results["partial_charges"], molecule.formal_charges
    )

    if os.path.exists(working_directory):
        shutil.rmtree(os.path.abspath(working_directory))

    # get ffdata and ff_params from amber prmtop and inpcrd files
    parm = pmd.load_file(results["prmtop_file"], results["inpcrd_file"])
    ffdata, _ = parm_to_ffdata(parm)

    # Add am1bcc partial charges
    ffdata["partial_charges"] = round_partial_charges

    # Add formal charge
    ffdata["formal_charges"] = molecule.formal_charges

    if parm_save_path is None:
        # remove parm dir
        shutil.rmtree(save_path, ignore_errors=True)

    ffdata = convert_value_to_list(ffdata)

    return ffdata
