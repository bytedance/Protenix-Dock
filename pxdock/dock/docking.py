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

import logging
import multiprocessing
import os
from typing import List

from pxdock.common import kWorkDir, my_random_string
from pxdock.dock.utils import kScoreConfigFile
from pxdock.engine import ReusableEngine
from pxdock.pipeline.prepare_ligand import prepare_ligand
from pxdock.pipeline.prepare_receptor import prepare_receptor
from pxdock.pipeline.get_pocket_info import compute_pocket_box
from pxdock.vina.prepare_receptor import prepare_receptor_pdbqt
from pxdock.vina.vina_docking import run_vina_docking


class ProtenixDock(object):
    def __init__(
        self,
        receptor_pdb_file: str,
        cofactor_fpath: str = None,
        sf_file: str = None,
        nthreads: int = 0,
    ):
        self._receptor_pdb = receptor_pdb_file
        self._cofactor_fpath = cofactor_fpath
        self._sf_file = kScoreConfigFile if not sf_file else sf_file
        self._nthreads = nthreads

    def set_box(
        self, pocket_center: List[float], box_size: List[float], buffer: float = 0
    ):
        self.pocket_center = pocket_center
        self.box_size = [b + buffer for b in box_size]
        logging.info(
            "pocket_center: %s"
            % (",".join([str(round(a, 3)) for a in self.pocket_center]))
        )
        logging.info(
            "box_size: %s" % (",".join([str(round(a, 3)) for a in self.box_size]))
        )
        receptor_json = prepare_receptor(
            self._receptor_pdb,
            pocket_center,
            self.box_size,
            cofactor_fpath=self._cofactor_fpath,
        )
        if not receptor_json:
            raise RuntimeError(f"prepare_receptor failed {self._receptor_pdb}")
        self._receptor_pdbqt = prepare_receptor_pdbqt(self._receptor_pdb, kWorkDir)
        self._engine = ReusableEngine(receptor_json, self._sf_file, self._nthreads)
        self._engine.set_box(*self.pocket_center, *self.box_size)

    def autobox(self, ligand_path: str, buffer: float = 5.0):
        box_center, box_size = compute_pocket_box(ligand_path, buffer=buffer)
        self.set_box(box_center, box_size)

    def generate_cache_maps(
        self,
        spacing: float = 0.175,
        cache_max_size: float = 0.0,
        output_dir: str = None,
    ) -> None:
        self._engine.generate_maps(spacing, cache_max_size)
        if output_dir is not None:
            self._engine.dump_maps(output_dir)

    def restore_cache_maps(self, output_dir: str = None) -> None:
        self._engine.load_maps(output_dir)

    def drop_cache_maps(self) -> None:
        self._engine.drop_maps()

    def run_pose_opt(
        self,
        ligand_sdf_file: str,
        out_dir: str = None,
        penalty: float = 1e6,
        vina_nposes: int = 50,
    ):
        ligand_jsons = self._prepare_json(
            ligand_sdf_file, mode="pose_opt", vina_nposes=vina_nposes
        )

        ligand_index_file = os.path.join(
            kWorkDir, f"{my_random_string()}_ligands.index"
        )
        with open(ligand_index_file, "w") as f:
            for ligand_json in ligand_jsons:
                ligand_abs_path = os.path.abspath(ligand_json)
                f.write(ligand_abs_path)
                f.write("\n")
        if out_dir is None:
            out_dir = os.path.join(kWorkDir, f"{my_random_string()}_pose_opt_out")
        self._engine.optimize(ligand_index_file, out_dir, slope=penalty)

        return out_dir

    def run_docking(
        self,
        ligand_sdf_file: str,
        out_dir: str = None,
        penalty: float = 1e6,
        seed: int = 101,
        num_walker: int = 256,
        relax_nsteps: int = 40,
        mc_prune_energy_threshold: float = -1,
    ):
        ligand_jsons = self._prepare_json(
            ligand_sdf_file,
            mode="mc",
            num_conf=num_walker,  # if mc_prune_energy_threshold > 0 else -1
        )

        ligand_index_file = os.path.join(
            kWorkDir, f"{my_random_string()}_ligands.index"
        )
        with open(ligand_index_file, "w") as f:
            for ligand_json in ligand_jsons:
                ligand_abs_path = os.path.abspath(ligand_json)
                f.write(ligand_abs_path)
                f.write("\n")
        if out_dir is None:
            out_dir = os.path.join(kWorkDir, f"{my_random_string()}_docking_out")
        self._engine.search(
            ligand_index_file,
            out_dir,
            seed,
            num_walker,
            relax_nsteps=relax_nsteps,
            mc_mmenergy_threshold=mc_prune_energy_threshold,
            slope=penalty,
        )

        return out_dir

    def _prepare_json(
        self,
        ligand_sdf: str,
        mode: str = "pose_opt",
        vina_nposes: int = 50,
        num_conf: int = -1,
    ):
        if not os.path.exists(ligand_sdf):
            raise RuntimeError(f"input ligand_sdf is invalid: {ligand_sdf}")
        if mode == "mc":
            return self._process_ligand(ligand_sdf, num_conf=num_conf)
        if mode != "pose_opt":
            raise RuntimeError(f"invalid mode: {mode}")
        kVinaOpts = {
            "cpu": int(multiprocessing.cpu_count()),
            "seed": 509,
            "exhaustiveness": 32,
            "search_depth": 0,
            "energy_range": 1000,
        }
        tmp_csv_path = os.path.join(kWorkDir, f"{my_random_string()}_vina_docking.csv")
        _ = run_vina_docking(
            input_sdf_path=ligand_sdf,
            receptor_pdbqt=self._receptor_pdbqt,
            box_center=self.pocket_center,
            box_size=self.box_size,
            n_poses=vina_nposes,
            output_csv_path=tmp_csv_path,
            **kVinaOpts,
        )
        return self._process_ligand(tmp_csv_path)

    def _process_ligand(self, ligand_file: str, num_conf: int = -1):
        ligand_jsons = prepare_ligand(ligand_file, num_conf=num_conf)
        return ligand_jsons
