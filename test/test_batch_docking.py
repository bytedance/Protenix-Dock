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
import tempfile
import unittest

import pytest
from pxdock import screen_ligands, dock_against_pockets
from pxdock.common import get_logger
from pxdock.dock.utils import parse_pocket_config

logger = get_logger(__name__)


class BatchDockingTestBase(unittest.TestCase):
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


class ScreenLigandsTest(BatchDockingTestBase):

    @pytest.mark.slow
    def test_screen_single_ligand(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results = screen_ligands(
                receptor_pdb=self.receptor_pdb,
                pocket_center=self.box_center,
                box_size=self.box_size,
                ligand_files=[self.ligand_sdf],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["ligand"], self.ligand_sdf)
            self.assertIsNotNone(results[0]["out_dir"])
            self.assertTrue(os.path.isdir(results[0]["out_dir"]))

    @pytest.mark.slow
    def test_screen_multiple_ligands_reuses_receptor(self):
        """Dock the same ligand file twice to verify receptor reuse works."""
        with tempfile.TemporaryDirectory() as tmpdir:
            results = screen_ligands(
                receptor_pdb=self.receptor_pdb,
                pocket_center=self.box_center,
                box_size=self.box_size,
                ligand_files=[self.ligand_sdf, self.ligand_sdf],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 2)
            for entry in results:
                self.assertIsNotNone(entry["out_dir"])
                self.assertTrue(os.path.isdir(entry["out_dir"]))
            # Output dirs should be different
            self.assertNotEqual(results[0]["out_dir"], results[1]["out_dir"])

    @pytest.mark.slow
    def test_screen_invalid_ligand_continues(self):
        """Invalid ligand paths should be caught and not crash the batch."""
        with tempfile.TemporaryDirectory() as tmpdir:
            results = screen_ligands(
                receptor_pdb=self.receptor_pdb,
                pocket_center=self.box_center,
                box_size=self.box_size,
                ligand_files=["/nonexistent/fake.sdf", self.ligand_sdf],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 2)
            self.assertIsNone(results[0]["out_dir"])
            self.assertIn("error", results[0])
            self.assertIsNotNone(results[1]["out_dir"])


class DockAgainstPocketsTest(BatchDockingTestBase):

    @pytest.mark.slow
    def test_single_pocket_tuple(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results = dock_against_pockets(
                ligand_sdf=self.ligand_sdf,
                pockets=[
                    (self.receptor_pdb, self.box_center, self.box_size),
                ],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 1)
            self.assertIn("pocket_0", results)
            self.assertIsNotNone(results["pocket_0"])
            self.assertTrue(os.path.isdir(results["pocket_0"]))

    @pytest.mark.slow
    def test_single_pocket_dict(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results = dock_against_pockets(
                ligand_sdf=self.ligand_sdf,
                pockets=[
                    {
                        "receptor": self.receptor_pdb,
                        "center": self.box_center,
                        "box_size": self.box_size,
                        "name": "active_site",
                    },
                ],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 1)
            self.assertIn("active_site", results)
            self.assertIsNotNone(results["active_site"])

    @pytest.mark.slow
    def test_multiple_pockets_same_receptor(self):
        """Dock against two different boxes on the same receptor."""
        shifted_center = [
            self.box_center[0] + 5,
            self.box_center[1],
            self.box_center[2],
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            results = dock_against_pockets(
                ligand_sdf=self.ligand_sdf,
                pockets=[
                    {
                        "receptor": self.receptor_pdb,
                        "center": self.box_center,
                        "box_size": self.box_size,
                        "name": "site_A",
                    },
                    {
                        "receptor": self.receptor_pdb,
                        "center": shifted_center,
                        "box_size": self.box_size,
                        "name": "site_B",
                    },
                ],
                cache_spacing=0.5,
                num_walker=16,
                out_dir=tmpdir,
            )
            self.assertEqual(len(results), 2)
            self.assertIn("site_A", results)
            self.assertIn("site_B", results)
            for name, out in results.items():
                self.assertIsNotNone(out, f"Pocket {name} failed")


if __name__ == "__main__":
    unittest.main()
