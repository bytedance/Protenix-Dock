#!/usr/bin/env python
"""Benchmark the C++ docking engine across multiple systems.

Usage:
    python test/bench_engine.py
"""
import os
import shutil
import time

from pxdock import ProtenixDock
from pxdock.dock.utils import parse_pocket_config

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.join(SCRIPT_DIR, "..", "examples")
OUT_DIR = "/tmp/pxdock_bench_engine"

PGO_TRAINING = os.environ.get("PXDOCK_PGO_TRAINING") == "1"
N_RUNS = 1 if PGO_TRAINING else 3
NUM_WALKERS = 4

SYSTEMS = [
    "5s8i",
    "8ay3",
    "7bmi",
    "5sb2",
    "7t0d",
]


def get_paths(system_id):
    return (
        os.path.join(EXAMPLES_DIR, f"{system_id}_protein_prepared.pdb"),
        os.path.join(EXAMPLES_DIR, f"{system_id}_ligand_prepared.sdf"),
        os.path.join(EXAMPLES_DIR, f"{system_id}_pocket.config"),
    )


def run_system(system_id):
    receptor_pdb, ligand_sdf, pocket_cfg = get_paths(system_id)
    pocket = parse_pocket_config(pocket_cfg)
    box_center = (pocket["center_x"], pocket["center_y"], pocket["center_z"])
    box_size = (pocket["size_x"], pocket["size_y"], pocket["size_z"])

    dock = ProtenixDock(receptor_pdb)

    t0 = time.perf_counter()
    dock.set_box(box_center, box_size)
    t_prep = time.perf_counter() - t0

    t0 = time.perf_counter()
    dock.run_docking(ligand_sdf, out_dir=OUT_DIR, num_walker=NUM_WALKERS)
    t_dock = time.perf_counter() - t0

    return t_prep, t_dock


def main():
    if PGO_TRAINING:
        print("PGO training mode: single iteration for profile data")
        for sid in SYSTEMS:
            if os.path.exists(get_paths(sid)[0]):
                run_system(sid)
        print("PGO training complete.")
        return

    # Warm up with first system
    print(f"Warming up ({SYSTEMS[0]})...")
    if os.path.exists(OUT_DIR):
        shutil.rmtree(OUT_DIR)
    run_system(SYSTEMS[0])

    # Benchmark each system
    all_results = {}
    for sid in SYSTEMS:
        receptor_pdb = get_paths(sid)[0]
        if not os.path.exists(receptor_pdb):
            print(f"\n  {sid}: SKIPPED (files not found)")
            continue

        prep_times = []
        dock_times = []
        for run_i in range(N_RUNS):
            if os.path.exists(OUT_DIR):
                shutil.rmtree(OUT_DIR)
            t_prep, t_dock = run_system(sid)
            prep_times.append(t_prep)
            dock_times.append(t_dock)

        avg_prep = sum(prep_times) / N_RUNS
        avg_dock = sum(dock_times) / N_RUNS
        min_dock = min(dock_times)
        all_results[sid] = (avg_prep, avg_dock, min_dock)
        print(f"  {sid:6s}  prep={avg_prep:.3f}s  dock_avg={avg_dock:.3f}s  dock_min={min_dock:.3f}s")

    # Summary
    total_dock_avg = sum(r[1] for r in all_results.values())
    total_dock_min = sum(r[2] for r in all_results.values())
    print(f"\n{'='*60}")
    print(f"Results ({N_RUNS} runs, {NUM_WALKERS} walkers, {len(all_results)} systems):")
    print(f"  {'system':8s} {'prep':>8s} {'dock_avg':>10s} {'dock_min':>10s}")
    print(f"  {'-'*8} {'-'*8} {'-'*10} {'-'*10}")
    for sid, (p, d, m) in all_results.items():
        print(f"  {sid:8s} {p:8.3f} {d:10.3f} {m:10.3f}")
    print(f"  {'-'*8} {'-'*8} {'-'*10} {'-'*10}")
    print(f"  {'TOTAL':8s} {'':8s} {total_dock_avg:10.3f} {total_dock_min:10.3f}")


if __name__ == "__main__":
    main()
