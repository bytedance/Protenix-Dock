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
import unittest

from pxdock.common import get_logger, my_random_string
from pxdock.dock.utils import kScoreConfigFile, parse_pocket_config
from pxdock.engine import ReusableEngine
from pxdock.pipeline.prepare_ligand import prepare_ligand
from pxdock.pipeline.prepare_receptor import prepare_receptor

logger = get_logger(__name__)


class DataParserTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        cls.ligand_sdf = os.path.join(
            current_dir, "../examples/5s8i_ligand_prepared.sdf"
        )
        cls.receptor_pdb = os.path.join(
            current_dir, "../examples/5s8i_protein_prepared.pdb"
        )
        cls.pocket_config = os.path.join(current_dir, "../examples/5s8i_pocket.config")

    def _test_prepare_ligand(self):
        prepared_jsons = prepare_ligand(self.ligand_sdf)

        if len(prepared_jsons) == 0:
            raise RuntimeError(f"prepare_ligand failed {self.ligand_sdf}")
        for ligand_json in prepared_jsons:
            with open(ligand_json, "r") as f:
                pre_ligand_data = json.load(f)
            logger.info(f"ligand keys: {pre_ligand_data.keys()}")
            logger.info(f"ligand ffdata keys: {pre_ligand_data['ffdata'].keys()}")
            logger.info(f"ligand geometry keys: {pre_ligand_data['geometry'].keys()}")
        return prepared_jsons

    def _test_prepare_receptor(self):
        pocket_config = parse_pocket_config(self.pocket_config)
        logger.info(f"receptor_path: {self.receptor_pdb}, {self.pocket_config}")
        prepared_json = prepare_receptor(
            self.receptor_pdb,
            box_center=[
                pocket_config["center_x"],
                pocket_config["center_y"],
                pocket_config["center_z"],
            ],
            box_size=[
                pocket_config["size_x"],
                pocket_config["size_y"],
                pocket_config["size_z"],
            ],
        )
        if prepared_json is None:
            raise RuntimeError(f"prepare_receptor failed {self.receptor_pdb}")
        else:
            logger.info(f"prepare_receptor result save to: {prepared_json}")
        return prepared_json, pocket_config

    def test_bdock(self):
        prepared_ligand_json = self._test_prepare_ligand()[0]
        prepared_receptor_json, box = self._test_prepare_receptor()
        if not (
            os.path.exists(prepared_ligand_json)
            and os.path.exists(prepared_receptor_json)
        ):
            raise RuntimeError(
                f"invalid ligand or receptor josn: {prepared_ligand_json} , {prepared_receptor_json}"
            )
        ligand_abs_path = os.path.abspath(prepared_ligand_json)
        tmp_work_dir = f"/tmp/{my_random_string()}"
        os.makedirs(tmp_work_dir, exist_ok=True)
        ligand_index_file = os.path.join(
            tmp_work_dir, f"{my_random_string()}_ligands.index"
        )
        with open(ligand_index_file, "w") as f:
            f.write(ligand_abs_path)
        penalty = 100.0
        seed = 101
        exhaustiveness = 32
        relax_nsteps = 40
        nthreads = 40
        output_dir = f"{tmp_work_dir}/bdock_out"

        engine = ReusableEngine(prepared_receptor_json, kScoreConfigFile, nthreads)
        engine.set_box(
            box["center_x"],
            box["center_y"],
            box["center_z"],
            box["size_x"],
            box["size_y"],
            box["size_z"],
        )
        engine.search(
            ligand_index_file,
            output_dir,
            seed,
            exhaustiveness,
            relax_nsteps=relax_nsteps,
            slope=penalty,
        )


if __name__ == "__main__":
    unittest.main()
