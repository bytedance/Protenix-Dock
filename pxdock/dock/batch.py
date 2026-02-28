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
import os
from typing import Dict, List, Optional, Tuple, Union

from pxdock.dock.docking import ProtenixDock

logger = logging.getLogger(__name__)


def screen_ligands(
    receptor_pdb: str,
    pocket_center: List[float],
    box_size: List[float],
    ligand_files: List[str],
    cofactor_fpath: str = None,
    sf_file: str = None,
    nthreads: int = 0,
    cache_spacing: float = 0.175,
    seed: int = 101,
    num_walker: int = 256,
    relax_nsteps: int = 40,
    mc_prune_energy_threshold: float = -1,
    include_affinity: bool = False,
    out_dir: str = None,
) -> List[dict]:
    """Screen multiple ligands against a single receptor pocket.

    Loads the receptor once, generates cache maps once, then docks each
    ligand file reusing the same ProtenixDock instance.

    Args:
        receptor_pdb: Path to receptor PDB file.
        pocket_center: [x, y, z] center of the docking box.
        box_size: [sx, sy, sz] size of the docking box.
        ligand_files: List of ligand SDF file paths. Each file may
            contain one or more molecules.
        cofactor_fpath: Optional cofactor file path.
        sf_file: Optional scoring function config file.
        nthreads: Number of engine threads (0 = auto).
        cache_spacing: Grid spacing for receptor cache maps.
        seed: Random seed for Monte Carlo search.
        num_walker: Number of MC walkers per ligand.
        relax_nsteps: Relaxation steps per MC walker.
        mc_prune_energy_threshold: Prune threshold (-1 = disabled).
        include_affinity: Include binding affinity (bscore) in output.
        out_dir: Base output directory. Per-ligand results are stored
            in subdirectories named after the ligand file.

    Returns:
        List of dicts, one per ligand file, each with keys "ligand"
        (input path) and "out_dir" (output path, or None on failure).
    """
    dock = ProtenixDock(
        receptor_pdb,
        cofactor_fpath=cofactor_fpath,
        sf_file=sf_file,
        nthreads=nthreads,
    )
    dock.set_box(pocket_center, box_size)
    dock.generate_cache_maps(spacing=cache_spacing)

    if out_dir is None:
        from pxdock.common import kWorkDir, my_random_string
        out_dir = os.path.join(kWorkDir, f"{my_random_string()}_screen_out")
    os.makedirs(out_dir, exist_ok=True)

    results = []
    for idx, lig_path in enumerate(ligand_files):
        lig_name = os.path.splitext(os.path.basename(lig_path))[0]
        lig_out = os.path.join(out_dir, f"{idx}_{lig_name}")
        logger.info(f"Docking ligand [{idx}]: {lig_path}")
        try:
            dock.run_docking(
                lig_path,
                out_dir=lig_out,
                seed=seed,
                num_walker=num_walker,
                relax_nsteps=relax_nsteps,
                mc_prune_energy_threshold=mc_prune_energy_threshold,
                include_affinity=include_affinity,
            )
            results.append({"ligand": lig_path, "out_dir": lig_out})
        except Exception as e:
            logger.error(f"Failed to dock {lig_path}: {e}")
            results.append({"ligand": lig_path, "out_dir": None, "error": str(e)})

    dock.drop_cache_maps()
    return results


PocketDef = Tuple[str, List[float], List[float]]


def dock_against_pockets(
    ligand_sdf: str,
    pockets: List[Union[PocketDef, Dict[str, object]]],
    cofactor_fpath: str = None,
    sf_file: str = None,
    nthreads: int = 0,
    cache_spacing: float = 0.175,
    seed: int = 101,
    num_walker: int = 256,
    relax_nsteps: int = 40,
    mc_prune_energy_threshold: float = -1,
    include_affinity: bool = False,
    out_dir: str = None,
) -> Dict[str, str]:
    """Dock one ligand against multiple receptor pockets.

    Each pocket can be a different receptor PDB or a different box on
    the same receptor. A fresh ProtenixDock instance is created for
    each pocket definition.

    Args:
        ligand_sdf: Path to ligand SDF file.
        pockets: List of pocket definitions. Each element is either:
            - A tuple of (receptor_pdb, center, box_size)
            - A dict with keys "receptor", "center", "box_size" and
              optionally "name" (used for the output subdirectory).
        cofactor_fpath: Optional cofactor file path.
        sf_file: Optional scoring function config file.
        nthreads: Number of engine threads (0 = auto).
        cache_spacing: Grid spacing for receptor cache maps.
        seed: Random seed for Monte Carlo search.
        num_walker: Number of MC walkers.
        relax_nsteps: Relaxation steps per MC walker.
        mc_prune_energy_threshold: Prune threshold (-1 = disabled).
        include_affinity: Include binding affinity (bscore) in output.
        out_dir: Base output directory. Per-pocket results are stored
            in subdirectories.

    Returns:
        Dict mapping pocket name/index to its output directory.
    """
    if out_dir is None:
        from pxdock.common import kWorkDir, my_random_string
        out_dir = os.path.join(kWorkDir, f"{my_random_string()}_multipocket_out")
    os.makedirs(out_dir, exist_ok=True)

    results = {}
    for idx, pocket in enumerate(pockets):
        if isinstance(pocket, dict):
            rec_pdb = pocket["receptor"]
            center = pocket["center"]
            bsize = pocket["box_size"]
            name = pocket.get("name", f"pocket_{idx}")
        else:
            rec_pdb, center, bsize = pocket
            name = f"pocket_{idx}"

        pocket_out = os.path.join(out_dir, name)
        logger.info(f"Docking against pocket {name}: receptor={rec_pdb}, center={center}")

        try:
            dock = ProtenixDock(
                rec_pdb,
                cofactor_fpath=cofactor_fpath,
                sf_file=sf_file,
                nthreads=nthreads,
            )
            dock.set_box(center, bsize)
            dock.generate_cache_maps(spacing=cache_spacing)
            dock.run_docking(
                ligand_sdf,
                out_dir=pocket_out,
                seed=seed,
                num_walker=num_walker,
                relax_nsteps=relax_nsteps,
                mc_prune_energy_threshold=mc_prune_energy_threshold,
                include_affinity=include_affinity,
            )
            dock.drop_cache_maps()
            results[name] = pocket_out
        except Exception as e:
            logger.error(f"Failed to dock against pocket {name}: {e}")
            results[name] = None

    return results
