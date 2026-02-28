#!/usr/bin/env python
"""Profile the Protenix-Dock pipeline with cProfile.

Usage:
    python test/profile_docking.py          # run profiling, save stats
    snakeviz /tmp/pxdock_profile.prof       # visualize in browser
"""

import cProfile
import os
import pstats

from pxdock import ProtenixDock
from pxdock.dock.utils import parse_pocket_config

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIGAND_SDF = os.path.join(SCRIPT_DIR, "../examples/5s8i_ligand_prepared.sdf")
RECEPTOR_PDB = os.path.join(SCRIPT_DIR, "../examples/5s8i_protein_prepared.pdb")
POCKET_CFG = os.path.join(SCRIPT_DIR, "../examples/5s8i_pocket.config")
PROFILE_OUT = os.getenv("PXDOCK_PROFILE_OUT", "/tmp/pxdock_profile.prof")


def run_pipeline():
    pocket = parse_pocket_config(POCKET_CFG)
    box_center = (pocket["center_x"], pocket["center_y"], pocket["center_z"])
    box_size = (pocket["size_x"], pocket["size_y"], pocket["size_z"])

    dock = ProtenixDock(RECEPTOR_PDB)
    dock.set_box(box_center, box_size)
    dock.run_docking(LIGAND_SDF, num_walker=4)


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    run_pipeline()
    profiler.disable()

    profiler.dump_stats(PROFILE_OUT)
    print(f"\nProfile saved to {PROFILE_OUT}")
    print("Run: snakeviz", PROFILE_OUT)
    print("\nTop 30 cumulative:")
    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    stats.print_stats(30)
