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

import inspect
import json
import os
import tempfile
import time
import traceback

import pandas as pd
from meeko import MoleculePreparation
from rdkit import Chem
from tqdm import tqdm

from pxdock.common import get_logger

logger = get_logger(__name__)


"""
    Executes a molecular docking simulation using AutoDock Vina.

    Parameters:
        input_sdf_path (str): Path to the ligand file in SDF format, which contains the 3D structure 
                        of small molecules to be docked. If the file contains multiple structures,
                        they will be docked one by one.
        receptor_pdbqt (str): Path to the receptor file in PDBQT format. 
        box_center (list or tuple): Three-dimensional coordinates (x, y, z) specifying the center 
                                    of the docking box.
        box_size (list or tuple): Three-dimensional size (width, height, depth) of the docking box.
        **task_configs (dict): Configuration options for the docking simulation. It allows for overriding
                                the default settings specified in 'default_task_configs'.

    Returns:
        conforms (list): [conform_{i} for i in range(number of mols in input_sdf_path)].
        scores (list): Corresponding docking scores for each mol in input_sdf_path.
        poses (list): The pdbqt string of each docked ligand pose.
        pose_i_scores (list of list): Docking scores for each pose of each mol.
        dock_seconds (list)
        errors (list)

"""

default_task_configs = {
    "score_func": "vina",
    "seed": 12345,
    "mode": "dock",
    "exhaustiveness": 32,
    "save_pose": True,
    "cpu": 8,
    "n_poses": 20,
    "search_depth": 0,  # 0: default heuristic
    "energy_range": 100,
}


class VinaDock(object):
    def __init__(self, pdbqt_list, prot_pdbqt, nconfs=-1, is_service=False):
        self.is_service = is_service
        self.prot_pdbqt = prot_pdbqt
        self.lig_pdbqt = pdbqt_list

    def set_box(self, pocket_center, box_size, buffer=0):
        self.pocket_center = pocket_center
        self.box_size = [b + buffer for b in box_size]
        logger.info(
            "pocket_center: %s"
            % (",".join([str(round(a, 3)) for a in self.pocket_center]))
        )
        logger.info(
            "box_size: %s" % (",".join([str(round(a, 3)) for a in self.box_size]))
        )

    def dock(
        self,
        score_func="vina",
        seed=0,
        mode="dock",
        exhaustiveness=32,
        search_depth=0,
        save_pose=True,
        cpu=8,
        n_poses=5,
        energy_range=3,
    ):
        if score_func != "vina":
            raise ValueError("VinaDock only support vina score function.")

        scores = []
        poses = []
        pose_i_scores = []
        conforms = []
        errors = []
        dock_seconds = []

        from vina import Vina

        v = Vina(sf_name=score_func, cpu=cpu, seed=seed)
        v.set_receptor(self.prot_pdbqt)
        v.compute_vina_maps(center=self.pocket_center, box_size=self.box_size)

        dock_params = inspect.signature(v.dock).parameters

        count = 0
        pbar = tqdm(self.lig_pdbqt)
        for lig in pbar:
            pbar.set_description(f"Vina: start to work on {count+1}-th molecule.")
            score = None
            pose = None
            pose_score = None
            error = None
            cur_time = time.time()

            try:
                if lig is None:
                    raise ValueError("input_pdbqt is None.")

                v.set_ligand_from_string(lig)

                if mode == "score_only":
                    score = v.score()[0]
                elif mode == "minimize":
                    score = v.optimize()[0]
                elif mode == "dock":
                    if "search_depth" in dock_params:
                        v.dock(
                            exhaustiveness=exhaustiveness,
                            n_poses=n_poses,
                            search_depth=search_depth,
                        )  # vina scm version
                    else:
                        v.dock(
                            exhaustiveness=exhaustiveness,
                            n_poses=n_poses,
                        )  # vina official version
                    score = v.energies(n_poses=1)[0][0]
                else:
                    raise ValueError("dock(): unknown dock mode %s" % (mode))

                if save_pose:
                    if mode == "score_only":
                        pose = None
                    elif mode == "minimize":
                        tmp = tempfile.NamedTemporaryFile()
                        with open(tmp.name, "w") as f:
                            v.write_pose(tmp.name, overwrite=True)
                        with open(tmp.name, "r") as f:
                            pose = f.read()
                    elif mode == "dock":
                        pose = v.poses(n_poses=n_poses, energy_range=energy_range)
                        pose_score = v.energies(
                            n_poses=n_poses, energy_range=energy_range
                        )
                        pose_score = [ps[0] for ps in pose_score]
                    else:
                        raise ValueError("dock(): unknown dock mode %s" % (mode))
            except Exception as exc:
                tb = traceback.format_exc()
                error = f"{exc}\n {tb}"
                score = None
                logger.info(f"failed to vina dock, exception: {error}")

            scores.append(score)
            poses.append(pose)
            pose_i_scores.append(pose_score)
            conforms.append(f"conform_{count}")
            errors.append(error)
            dock_seconds.append(round(time.time() - cur_time, 3))

            count = count + 1

        return conforms, scores, poses, pose_i_scores, dock_seconds, errors


def prepare_ligand_pdbqt(sdf_path, nconfs=-1):
    lig_properties = []
    lig_pdbqts = []
    errors = []

    count = 0
    st = time.time()
    for rdm in Chem.ForwardSDMolSupplier(sdf_path, removeHs=False):
        if nconfs > 0 and count >= nconfs:
            break
        properties = None
        error = None
        pdbqt = None

        try:
            properties = rdm.GetPropsAsDict()

            if "ligand" not in properties:
                # Check if the "ligand" property is not in the properties dictionary
                # If not, use the "_Name" property as the ligand name
                properties["ligand"] = rdm.GetProp("_Name")

            # convert to pdbqt
            preparator = MoleculePreparation(add_index_map=True, rigid_macrocycles=True)
            preparator.prepare(rdm)
            pdbqt = preparator.write_pdbqt_string()

        except Exception as e:
            pdbqt = None
            error = e
            logger.info(
                f"failed to prepare ligand pdbqt for the {count}-th molecule, exception: {error}"
            )

        lig_pdbqts.append(pdbqt)
        lig_properties.append(properties)
        errors.append(error)
        count = count + 1

    tt = time.time() - st
    logger.info(
        f"meeko preparation: prepared {count} molecules in {tt:.3f} seconds (avg: {tt/count:.3f} seconds per mol)."
    )

    return lig_properties, lig_pdbqts, errors


def run_vina_docking(
    input_sdf_path,
    receptor_pdbqt,
    box_center,
    box_size,
    output_csv_path=None,
    nconfs=-1,
    **task_configs,
):
    if output_csv_path is not None:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)

    # setup task configs
    docking_task_configs = dict(default_task_configs)
    for k in docking_task_configs:
        if k in task_configs:
            docking_task_configs[k] = task_configs[k]

    assert (
        docking_task_configs["score_func"] == "vina"
    ), f"Currently only support score_func=vina."

    logger.info(f"run_vina_docking: {docking_task_configs}")
    logger.info("sdf_path is {}".format(input_sdf_path))

    # Prepare ligand PDBQT files and read ligand properties from the SDF file
    lig_properties, lig_pdbqts, prepare_errors = prepare_ligand_pdbqt(
        input_sdf_path, nconfs=nconfs
    )

    dock = VinaDock(
        pdbqt_list=lig_pdbqts,
        prot_pdbqt=receptor_pdbqt,
        nconfs=nconfs,
        is_service=False,
    )
    dock.set_box(box_center, box_size)

    conforms, scores, poses, pose_i_scores, dock_seconds, dock_errors = dock.dock(
        **docking_task_configs
    )

    # merge error message
    errors = []
    for e1, e2 in zip(prepare_errors, dock_errors):
        if e1 is not None:
            errors.append(f"Meeko preparation error: {e1}")
            continue
        if e2 is not None:
            errors.append(f"Vina dock error: {e2}")
            continue
        errors.append(None)

    # organize the Vina docking results into a DataFrame
    df = pd.DataFrame(lig_properties)  # collect ligand properties saved in sdf_path

    df["score"] = scores
    df["pose"] = poses
    df["pose_i_score"] = pose_i_scores
    df["vina_dock_seconds"] = dock_seconds
    df["vina_error"] = errors
    df["sdf_conform_id"] = conforms

    def to_json(x):
        if x is not None:
            x = json.dumps(x)
        return x

    df["pose_i_score"] = df["pose_i_score"].apply(to_json)

    # reorder the columns in the DataFrame
    specified_order = [
        "ligand",
        "sdf_conform_id",
        "score",
        "pose",
        "pose_i_score",
        "vina_dock_seconds",
        "vina_error",
    ]
    new_order = specified_order + [
        col for col in df.columns if col not in specified_order
    ]
    df = df[new_order]
    if output_csv_path is not None:
        df.to_csv(output_csv_path)
        logger.info(f"Vina docking results are saved to {output_csv_path}")

    return df
