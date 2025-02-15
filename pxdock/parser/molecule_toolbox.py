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

import math
import os
from collections import defaultdict

from pxdock.common import get_logger
from pxdock.common.utilities import TemplateRender

logger = get_logger(__name__)

# standard residue name
# REF: http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_79.html
standard_residue_name = {
    "ALA": "Alanine",
    "ARG": "Arginine",
    "ASN": "Asparagine",
    "ASP": "Aspartic acid",
    "ASH": "Aspartic acid neutral",
    "ASX": "ASP/ASN ambiguous",
    "CYS": "Cysteine",
    "CYX": "Cysteines linked via disulfide bonds",
    "CYM": "Cysteines with negative charge",
    "GLN": "Glutamine",
    "GLU": "Glutamic acid",
    "GLH": "Glutamic acid neutral",
    "GLX": "GLU/GLN ambiguous",
    "GLY": "Glycine",
    "HIS": "Histidine",
    "HIE": "Histidine with hydrogen on the epsilon nitrogen",
    "HID": "Histidine with hydrogen on the delta nitrogen",
    "HIP": "Histidine with hydrogens on both nitrogens; this is positively charged",
    "ILE": "Isoleucine",
    "LEU": "Leucine",
    "LYS": "Lysine",
    "MET": "Methionine",
    "PHE": "Phenylalanine",
    "PRO": "Proline",
    "SER": "Serine",
    "THR": "Threonine",
    "TRP": "Tryptophan",
    "TYR": "Tyrosine",
    "UNK": "Unknown",
    "VAL": "Valine",
}

supported_ions = ["CA", "CL", "NA", "MG", "K", "RB", "CS", "LI", "ZN", "II"]

water_residue = ["HOH", "SOL", "AION"]

supported_capping = ["NME", "ACE", "NMA", "NHE"]

ligand_name = ["LIG"]


def patch_PDB_record(atom_record: dict):
    """Update residue information to be accepted by AMBER tools

    Args:
        atom_record (dict): Atom record from PDB
    """
    if atom_record["residue_name"] == "NMA":
        atom_record["residue_name"] = "NME"
        if atom_record["name"] == "CA":
            atom_record["name"] = "C"
    return


def parse_PDB_ATOM_HETATM_record(raw_pdb_line: str) -> dict:
    """parse the line of raw pdb record

    Args:
        raw_pdb_line (str): a line of raw pdb record

    Returns:
        dict: parsed record
    """
    atom = {
        "record_name": str(raw_pdb_line[0:6]).strip(),
        "nr": int(raw_pdb_line[6:11]),
        "x": float(raw_pdb_line[30:38]) / 10.0,  # A to nm
        "y": float(raw_pdb_line[38:46]) / 10.0,  # A to nm
        "z": float(raw_pdb_line[46:54]) / 10.0,  # A to nm
        "name": str(raw_pdb_line[12:16]).strip().upper(),
        "altloc": str(raw_pdb_line[16:17]).strip(),
        "residue_name": str(raw_pdb_line[17:20]).strip().upper(),
        "residue_num": int(raw_pdb_line[22:26]),
        "chain_id": str(raw_pdb_line[21:22]).strip().upper(),
        "icode": str(raw_pdb_line[26:27]).strip(),
        "occupancy": float(raw_pdb_line[54:60]),
        "temp_factor": 0.0,
        "charge": 0.0,
        "charge_str": "",
    }
    try:
        atom["element"] = str(raw_pdb_line[76:78]).strip().upper()
    except Exception:
        atom["element"] = ""

    try:
        atom["charge_str"] = str(raw_pdb_line[78:80]).strip().upper()
        sign = 1.0
        if atom["charge_str"][1] == "-":
            sign = -1.0
        atom["charge"] = sign * float(atom["charge_str"][0])
    except Exception:
        atom["charge"] = 0.0
        atom["charge_str"] = ""

    # patch_PDB_record(atom)

    # if atom["residue_name"] in standard_residue_name.keys():
    #     # atom["record_name"] = "ATOM"
    #     return atom

    # atom["record_name"] = "HETATM"
    if atom["residue_name"] in water_residue:
        return atom
    if atom["name"] in supported_ions and atom["residue_name"] in supported_ions:
        atom["record_name"] = "HETATM"
        return atom
    if atom["residue_name"] in supported_capping:
        return atom
    if atom["residue_name"] in ligand_name:
        return atom

    if len(atom["residue_name"]) > 0:
        return atom

    raise ValueError(f"Invalid residue name {atom['residue_name']}")


def measure_atom_distance_pdb_record(atom1: dict, atom2: dict) -> float:
    x1 = atom1["x"]
    y1 = atom1["y"]
    z1 = atom1["z"]

    x2 = atom2["x"]
    y2 = atom2["y"]
    z2 = atom2["z"]

    dist = math.sqrt((x2 - x1) ** 2.0 + (y2 - y1) ** 2.0 + (z2 - z1) ** 2.0)
    return dist


def fix_residue_name_for_amber(residue_atoms: list) -> list:
    H_tol = 0.12
    atom0 = residue_atoms[0]
    if atom0["record_name"] != "ATOM":
        return residue_atoms
    # select alternative location of atom
    alt = defaultdict(lambda: 0.0)
    for atom in residue_atoms:
        if atom["altloc"] != " " and atom["altloc"] != "":
            if atom["altloc"] in alt.keys():
                alt[atom["altloc"]] = min(alt[atom["altloc"]], atom["occupancy"])
            else:
                alt[atom["altloc"]] = atom["occupancy"]

    if len(alt.keys()) > 0:
        select_key = list(alt.keys())[0]
        max_val = alt[select_key]
        for key, val in alt.items():
            if val > max_val:
                max_val = val
                select_key = key
        new_residue_atoms = []
        for atom in residue_atoms:
            if (
                atom["altloc"] == " "
                or atom["altloc"] == ""
                or atom["altloc"] == select_key
            ):
                new_residue_atoms.append(atom)
    else:
        new_residue_atoms = residue_atoms

    if new_residue_atoms[0]["residue_name"] == "HIS":
        HID = 0
        HIE = 0
        for atom in new_residue_atoms:
            if "ND1" in atom["name"]:
                for atom_H in new_residue_atoms:
                    if atom_H["element"] == "H":
                        dist = measure_atom_distance_pdb_record(atom, atom_H)
                        if dist < H_tol:
                            HID += 1
            if "NE2" in atom["name"]:
                for atom_H in new_residue_atoms:
                    if atom_H["element"] == "H":
                        dist = measure_atom_distance_pdb_record(atom, atom_H)
                        if dist < H_tol:
                            HIE += 1
        new_name = "HIS"
        if HID == 1 and HIE == 0:
            new_name = "HID"
        elif HID == 0 and HIE == 1:
            new_name = "HIE"
        elif HID == 1 and HIE == 1:
            new_name = "HIP"
        else:
            raise ValueError(
                f"Can not recongnize HIS, residue number: {new_residue_atoms[0]['residue_num']}"
            )

        for atom in new_residue_atoms:
            atom["residue_name"] = new_name

    if new_residue_atoms[0]["residue_name"] == "CYS":
        CYM = False
        for atom in new_residue_atoms:
            if "SG" in atom["name"]:
                if atom["charge"] < 0:
                    CYM = True
                    break

                CYM = True
                for atom_H in new_residue_atoms:
                    if atom_H["element"] == "H":
                        dist = measure_atom_distance_pdb_record(atom, atom_H)
                        if dist < 0.15:
                            CYM = False
                            break
                break
        if CYM:
            for atom in new_residue_atoms:
                atom["residue_name"] = "CYM"

    if new_residue_atoms[0]["residue_name"] == "GLU":
        GLH = False
        for atom in new_residue_atoms:
            if "OE2" in atom["name"]:
                for atom_H in new_residue_atoms:
                    if atom_H["element"] == "H":
                        dist = measure_atom_distance_pdb_record(atom, atom_H)
                        if dist < H_tol:
                            GLH = True
                            break
        if GLH:
            for atom in new_residue_atoms:
                atom["residue_name"] = "GLH"

    if new_residue_atoms[0]["residue_name"] == "ASP":
        ASH = 0
        for atom in new_residue_atoms:
            if "OD2" in atom["name"]:
                for atom_H in new_residue_atoms:
                    if atom_H["element"] == "H":
                        dist = measure_atom_distance_pdb_record(atom, atom_H)
                        if dist < H_tol:
                            ASH += 1

        new_name = "ASP"
        if ASH == 0:
            new_name = "ASP"
        elif ASH == 1:
            new_name = "ASH"
        else:
            raise ValueError(
                f"Can not recongnize ASP, residue number: {new_residue_atoms[0]['residue_num']}"
            )

        for atom in new_residue_atoms:
            atom["residue_name"] = new_name

    return new_residue_atoms


def check_CYX_SS_bond(residues: list, connections: list) -> None:
    for residue_A in residues:
        if (
            residue_A[0]["residue_name"] != "CYS"
            and residue_A[0]["residue_name"] != "CYM"
            and residue_A[0]["residue_name"] != "CYX"
        ):
            continue
        for atom1 in residue_A:
            if atom1["name"] == "SG":
                break
        if atom1["name"] != "SG":
            continue

        SS_bond = False
        for connection in connections:
            if atom1["nr"] == connection["nr"][0]:
                atom2_nr = connection["nr"][1]
                for residue_B in residues:
                    for atom_B in residue_B:
                        if atom_B["nr"] == atom2_nr and atom_B["name"] == "SG":
                            SS_bond = True
                            break
                    if SS_bond:
                        break
                if SS_bond:
                    break

        if SS_bond == False:
            continue
        for atom in residue_A:
            atom["residue_name"] = "CYX"
        for atom in residue_B:
            atom["residue_name"] = "CYX"

    return


def generate_sanitized_protein_pdb(
    src_pdb: str,
    output_dir: str,
    output_pdb_filename: str,
    remove_hydrogen: bool = True,
):
    """Generate a sanitized pdb and remove SOL, IONS, Hydrogen etc.

    Args:
        src_pdb (str): source pdb file
        output_dir (str): dir to output pdb
        output_pdb_filename (str): sanitized pdb filename
        remove_hydrogen (bool): remove hydrogen in pdb
    """

    def add_atom_to_chain(residue_atoms: list, chain: list):
        for residue_atom in residue_atoms:
            if (
                residue_atom["record_name"] != "ATOM"
                and residue_atom["residue_name"] not in supported_capping
            ):
                chain.append(residue_atom)
            elif remove_hydrogen and residue_atom["element"] == "H":
                continue
            else:
                chain.append(residue_atom)
        return

    with open(src_pdb, "r", encoding="utf-8") as fin:
        src_pdb = fin.readlines()

    pdb_length = len(src_pdb)
    protein_chains = []
    chain = []
    old_connections = []
    pre_residue_id = None
    residue_atoms = []
    residues = []

    for idx, line in enumerate(src_pdb):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_record = parse_PDB_ATOM_HETATM_record(line)
            if (
                atom_record["residue_name"] in supported_ions
                or atom_record["residue_name"] in water_residue
            ):
                continue

            residue_id = (atom_record["residue_num"], atom_record["icode"])
            if pre_residue_id != residue_id and len(residue_atoms) > 0:
                residue_atoms = fix_residue_name_for_amber(residue_atoms)
                residues.append(residue_atoms)
                add_atom_to_chain(residue_atoms, chain)

            if len(chain) > 0:
                if atom_record["chain_id"] != chain[-1]["chain_id"] or (
                    pre_residue_id != None
                    and abs(atom_record["residue_num"] - pre_residue_id[0]) > 1
                ):
                    protein_chains.append(chain)
                    chain = []

            if pre_residue_id != residue_id:
                pre_residue_id = residue_id
                residue_atoms = []

            residue_atoms.append(atom_record)

        if line.startswith("TER") or line.startswith("END") or idx == pdb_length - 1:
            if len(residue_atoms) > 0:
                residue_atoms = fix_residue_name_for_amber(residue_atoms)
                residues.append(residue_atoms)
                add_atom_to_chain(residue_atoms, chain)
                pre_residue_id = None
                residue_atoms = []

            if len(chain) > 0:
                protein_chains.append(chain)
                chain = []

            continue

        if line.startswith("CONECT"):
            connect_record = [line[0:6]]
            sub_len = 5
            for i in range(6, len(line), sub_len):
                if i + sub_len > len(line):
                    break
                connect_record.append(line[i : i + sub_len])
            connect_nr = [int(item) for item in connect_record[1:]]
            connect = {"record_name": "CONECT", "nr": connect_nr}
            old_connections.append(connect)

    # check CYX S-S bond
    check_CYX_SS_bond(residues, old_connections)

    residue_renum = {}
    for idx, residue in enumerate(residues):
        for atom in residue:
            residue_renum[idx + 1] = atom["residue_num"]
            atom["residue_num"] = idx + 1
            atom["icode"] = ""

    nr_old_to_new = {}

    nr = 1
    for chain in protein_chains:
        for record in chain:
            old_nr = record["nr"]
            nr_old_to_new[old_nr] = nr
            record["nr"] = nr
            nr += 1

        ter_record = {
            "record_name": "TER",
            "name": "",
            "nr": nr,
            "residue_name": chain[-1]["residue_name"],
            "chain_id": chain[-1]["chain_id"],
            "residue_num": chain[-1]["residue_num"],
            "icode": chain[-1]["icode"],
        }
        chain.append(ter_record)
        nr += 1

    new_connections = []
    for connection in old_connections:
        nr_all_in_map = True
        for nr in connection["nr"]:
            if nr not in nr_old_to_new.keys():
                nr_all_in_map = False
                break
        if nr_all_in_map:
            for idx, old_nr in enumerate(connection["nr"]):
                connection["nr"][idx] = nr_old_to_new[old_nr]
            new_connections.append(connection)

    sanitized_pdb_record = []
    for chain in protein_chains:
        for record in chain:
            sanitized_pdb_record.append(record)
    for connection in new_connections:
        sanitized_pdb_record.append(connection)

    rend_dict = {"records": sanitized_pdb_record}
    tempate_dir = os.path.join(os.path.dirname(__file__), "../data/templates")
    render = TemplateRender(tempate_dir, "sanitized_protein.pdb")
    render.rend_template(
        rend_dict=rend_dict, output_dir=output_dir, filename=output_pdb_filename
    )
    return
