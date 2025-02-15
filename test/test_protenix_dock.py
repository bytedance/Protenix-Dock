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
import unittest

import pytest
from pxdock import ProtenixDock
from pxdock.common import get_logger
from pxdock.dock.utils import parse_pocket_config

logger = get_logger(__name__)


class ProtenixDockTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        cls.ligand_sdf = os.path.join(
            current_dir, "../examples/5s8i_ligand_prepared.sdf"
        )
        cls.receptor_pdb = os.path.join(
            current_dir, "../examples/5s8i_protein_prepared.pdb"
        )
        pocket_config = parse_pocket_config(
            os.path.join(current_dir, "../examples/5s8i_pocket.config")
        )
        cls.box_center = [
            pocket_config["center_x"],
            pocket_config["center_y"],
            pocket_config["center_z"],
        ]
        cls.box_size = [
            pocket_config["size_x"],
            pocket_config["size_y"],
            pocket_config["size_z"],
        ]

    @pytest.mark.slow
    def test_pose_opt_with_cache(self):
        dock_instance = ProtenixDock(self.receptor_pdb)
        dock_instance.set_box(self.box_center, self.box_size)
        dock_instance.generate_cache_maps(spacing=0.5)
        opt_res_files = dock_instance.run_pose_opt(self.ligand_sdf)
        logger.info(f"pose optimization without cache result: {opt_res_files}")
        dock_instance.drop_cache_maps()

    @pytest.mark.slow
    def test_pose_opt_no_cache(self):
        dock_instance = ProtenixDock(self.receptor_pdb)
        dock_instance.set_box(self.box_center, self.box_size)
        opt_res_files = dock_instance.run_pose_opt(self.ligand_sdf)
        logger.info(f"pose optimization without cache result: {opt_res_files}")

    @pytest.mark.slow
    def test_mc_docking_with_cache(self):
        dock_instance = ProtenixDock(self.receptor_pdb)
        dock_instance.set_box(self.box_center, self.box_size)
        dock_instance.generate_cache_maps(spacing=0.5)
        docking_res_files = dock_instance.run_docking(self.ligand_sdf, num_walker=16)
        logger.info(f"docking with cache result: {docking_res_files}")
        dock_instance.drop_cache_maps()

    @pytest.mark.slow
    def test_mc_docking_no_cache(self):
        dock_instance = ProtenixDock(self.receptor_pdb)
        dock_instance.set_box(self.box_center, self.box_size)
        docking_res_files = dock_instance.run_docking(self.ligand_sdf, num_walker=16)
        logger.info(f"docking without cache result: {docking_res_files}")

    @pytest.mark.slow
    def test_mc_prune_docking(self):
        dock_instance = ProtenixDock(self.receptor_pdb)
        dock_instance.set_box(self.box_center, self.box_size)
        dock_instance.generate_cache_maps(spacing=0.5)
        docking_res_files = dock_instance.run_docking(
            self.ligand_sdf, mc_prune_energy_threshold=500, num_walker=16
        )
        logger.info(f"prune_docking with cache result: {docking_res_files}")
        dock_instance.drop_cache_maps()


if __name__ == "__main__":
    unittest.main()
