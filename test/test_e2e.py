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
import shutil
import subprocess
import tempfile
import unittest

import pytest
from pxdock import ProtenixDock


@pytest.mark.e2e
class EndToEndTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        cls.ligand_sdf = os.path.join(
            current_dir, "../examples/5s8i_ligand_prepared.sdf"
        )
        cls.input_pdb = os.path.join(
            current_dir, "../examples/5s8i_protein_prepared.pdb"
        )
        cls.work_dir = tempfile.mkdtemp(prefix="pxdock_e2e_")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work_dir, ignore_errors=True)

    def test_full_pipeline(self):
        if not shutil.which("pdb2pqr30"):
            pytest.skip("pdb2pqr30 not found — install pdb2pqr to run e2e tests")
        if not shutil.which("tleap"):
            pytest.skip("tleap not found — install AmberTools to run e2e tests")

        # pdb2pqr step: re-prepare the PDB (tests the binary works)
        prepared_pdb = os.path.join(self.work_dir, "receptor_prepared.pdb")
        output_pqr = os.path.join(self.work_dir, "receptor.pqr")
        result = subprocess.run(
            [
                "pdb2pqr30",
                "--ff=AMBER",
                "--ffout=AMBER",
                "--pdb-output",
                prepared_pdb,
                self.input_pdb,
                output_pqr,
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertEqual(
            result.returncode,
            0,
            f"pdb2pqr30 failed:\nstdout: {result.stdout}\nstderr: {result.stderr}",
        )
        self.assertTrue(
            os.path.exists(prepared_pdb), "pdb2pqr30 did not produce output PDB"
        )

        # autobox + docking
        dock_instance = ProtenixDock(prepared_pdb)
        dock_instance.autobox(self.ligand_sdf)
        docking_res_files = dock_instance.run_docking(
            self.ligand_sdf, num_walker=16
        )
        self.assertIsNotNone(docking_res_files, "run_docking returned None")

        # score_complex
        scores = dock_instance.score_complex(self.ligand_sdf)
        self.assertIsInstance(scores, list)
        self.assertGreater(len(scores), 0)
        for entry in scores:
            self.assertIn("ligand_file", entry)
            self.assertIn("pose_id", entry)
            self.assertIn("pscore", entry)
            self.assertIn("bscore", entry)
            self.assertIsInstance(entry["pose_id"], int)
            self.assertIsInstance(entry["pscore"], float)


if __name__ == "__main__":
    unittest.main()
