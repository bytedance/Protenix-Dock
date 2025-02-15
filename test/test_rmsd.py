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
from typing import List

from pxdock.common import get_logger, my_random_string

logger = get_logger(__name__)


class RMSDTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        cls.json_output = os.path.join(
            current_dir, "../examples/5sak_docking_result.json"
        )
        cls.local_ref_sdf = os.path.join(current_dir, "../examples/5sak_ligand_ref.sdf")

    def calculate_pose_rmsd(
        self,
        mapped_smiles: str,
        xyz: List[List[float]],
        local_ref_file: str,
    ) -> float:
        from pxdock.common.rmsd_calculator import (
            calculate_rmsd,
            write_ligand_to_sdf,
        )

        predict_sdf_file = os.path.abspath(f"{my_random_string(15)}.sdf")
        write_ligand_to_sdf(mapped_smiles, [xyz], predict_sdf_file)
        value_list = calculate_rmsd(predict_sdf_file, local_ref_file)
        assert len(value_list) == 1, "Only one pose is expected!"
        os.remove(predict_sdf_file)
        return value_list[0]

    def test_calculate_rmsd(self):
        with open(self.json_output, "r") as f:
            json_out_data = json.load(f)
        best_id = json_out_data["best_pose"]["index"]
        best_pose = json_out_data["poses"][best_id]
        rmsd = self.calculate_pose_rmsd(
            json_out_data["mapped_smiles"],
            best_pose["ligand"]["xyz"],
            self.local_ref_sdf,
        )
        logger.info(f"rmsd result is {rmsd}")


if __name__ == "__main__":
    unittest.main()
