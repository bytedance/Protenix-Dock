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

import json
import os
import traceback

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from pxdock.common import get_logger
from pxdock.geometry.geometry import LigandGeometry
from pxdock.parser.ligand import LigandParser, generate_candidate_confs

logger = get_logger(__name__)


def parse_geometry(text: str) -> dict:
    ligand_data = json.loads(text)
    geo = LigandGeometry(ligand_data["bond_index"], ligand_data["is_rotatable"])
    ligand_data["geometry"] = {
        "num_atoms": geo.num_atoms.tolist(),  # numpy.int64 -> int
        "permute_index": geo.permute_index,
        # The followings are consumed in function `LigandGeometry.torsion6_to_xyz(...)`
        "frag_traverse_levels": geo.frag_traverse_levels,
        "edge_to_bond": geo.frag_tree.edge_to_bond,
        "frag_split_index": geo.frag_split_index,
    }
    return ligand_data


def prepare_ligand(
    ligand_file: str, out_dir: str=None, include_geometry: bool=True, debug: bool=False, **kwargs
) -> list[dict]:
    if not os.path.exists(ligand_file):
        raise RuntimeError(f"ligand file {ligand_file} not exists")
    if ligand_file.endswith(".sdf"):
        df = prepare_ligand_from_sdf(ligand_file, debug=debug, **kwargs)
    elif ligand_file.endswith(".csv"):
        df = prepare_ligand_from_pose(ligand_file, debug=debug, **kwargs)
    else:
        raise RuntimeError(
            f"The ligand_file should be in sdf or csv format, but got {ligand_file}"
        )

    # find failed rows
    df_failed = df[df["ligand_data"].isnull()]
    df = df[df["ligand_data"].notnull()]
    df = df.drop(columns=["ligand_data_error"])
    if include_geometry:
        df["ligand_data"] = df["ligand_data"].apply(parse_geometry)

    if out_dir is None or out_dir == "":
        out_json_basename = os.path.join(
            os.path.dirname(ligand_file),
            f"{os.path.splitext(os.path.basename(ligand_file))[0]}-prepared-ligand",
        )
    else:
        out_json_basename = os.path.join(
            out_dir,
            f"{os.path.splitext(os.path.basename(ligand_file))[0]}-prepared-ligand",
        )
    out_dir = os.path.abspath(os.path.dirname(out_json_basename))
    os.makedirs(out_dir, exist_ok=True)
    if len(df_failed) > 0:
        logger.warning(f"Failed to parse {len(df_failed)} ligands")

    ligand_data_size = len(df["ligand_data"].tolist())
    if ligand_data_size <= 0:
        raise RuntimeError(f"prepare_ligand_from_sdf failed. Check the log file.")
    prepared_jsons = []
    num_conf = kwargs.get("num_conf", -1)
    for idx, df_ligand_data in enumerate(df["ligand_data"].tolist()):
        out_json = f"{out_json_basename}-{idx}.json"
        if isinstance(df_ligand_data, str):
            ligand_json_data = json.loads(df_ligand_data)
        elif isinstance(df_ligand_data, dict):
            ligand_json_data = df_ligand_data
        else:
            raise RuntimeError(
                f"df_ligand_data should be `str`, but got {type(df_ligand_data)}"
            )
        if num_conf > 0:
            ligand_json_data["xyz"] = generate_candidate_confs(
                ligand_json_data["mapped_smiles"],
                num_candidate_confs=num_conf,
                ffopt=True,
            )
        with open(out_json, "w") as f:
            json.dump(ligand_json_data, f)
        prepared_jsons.append(out_json)
    return prepared_jsons


def prepare_ligand_from_sdf(sdf_path: str, debug: bool=False, **kwargs) -> pd.DataFrame:
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    datas = []
    for i, rkmol in tqdm(enumerate(suppl)):
        # bytedock ligand_data
        try:
            # name
            properties_dict = rkmol.GetPropsAsDict()
            if "ligand" in properties_dict:
                ligand = rkmol.GetProp("ligand")
            else:
                ligand = rkmol.GetProp("_Name")
            data = {"ligand": ligand}
            ligand_data = LigandParser(rkmol, **kwargs).get_data()
            data["mapped_smiles"] = ligand_data["mapped_smiles"]
            data["ligand_data"] = json.dumps(ligand_data)
            data["ligand_data_error"] = None
        except Exception as exc:
            tb = traceback.format_exc()
            logger.warning(f"Error in parsing ligand: {str(exc)} \n {tb}")
            data["mapped_smiles"] = None
            data["ligand_data"] = None
            data["ligand_data_error"] = f"{str(exc)} \n {tb}"

        datas.append(data)
        if debug and i == 5:
            break

    df = pd.DataFrame.from_dict(datas)

    # reorganize columns
    lead_columns = ["ligand", "mapped_smiles", "ligand_data", "ligand_data_error"]
    columns = lead_columns + [c for c in df.columns if c not in lead_columns]
    df = df[columns]

    return df


def prepare_ligand_from_pose(ligand_fpath: str, debug=False, **kwargs) -> pd.DataFrame:
    df = pd.read_csv(ligand_fpath, index_col=False)

    assert "ligand" in df.columns, f"ligand column is missing in {ligand_fpath}"
    assert "pose" in df.columns, f"pose column is missing in {ligand_fpath}"

    ligand_data_list = []
    mapped_smiles_list = []
    ligand_error_list = []
    ligand_pdbqt_list = df["pose"].tolist()

    # for ilig in range (1):
    for pdbqt in tqdm(ligand_pdbqt_list):
        try:
            if pdbqt is None:
                raise ValueError(f"pdbqt is None.")
            ligand_data = LigandParser(pdbqt, **kwargs).get_data()
            mapped_smiles_list.append(ligand_data["mapped_smiles"])
            ligand_data_list.append(json.dumps(ligand_data))
            ligand_error_list.append(None)
        except Exception as exc:
            tb = traceback.format_exc()
            logger.warning(f"Error in parsing ligand: {str(exc)} \n {tb}")
            ligand_data_list.append(None)
            mapped_smiles_list.append(None)
            ligand_error_list.append(f"{str(exc)} \n {tb}")

    df["mapped_smiles"] = mapped_smiles_list
    df["ligand_data"] = ligand_data_list
    df["ligand_data_error"] = ligand_error_list
    df = df.drop(columns=["pose"])

    # reorganize columns
    lead_columns = ["ligand", "mapped_smiles", "ligand_data", "ligand_data_error"]
    columns = lead_columns + [c for c in df.columns if c not in lead_columns]
    df = df[columns]

    return df
