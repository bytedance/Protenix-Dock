# Protenix-Dock

This repository hosts the source code for our work "Protenix-Dock: An accurate and trainable end-to-end protein-ligand docking framework using empirical scoring functions".  

For more information about the implementation and the performance of Protenix-Dock, see our [technical report](ProtenixDock_Technical_Report.pdf).

🔍 Protenix-Dock is a classical protein-ligand docking method designed for rigid docking tasks. For our deep learning complex structure prediction model, see [Protenix](https://github.com/bytedance/Protenix).

## Features

✨ Advanced docking conformation sampling.

✨ Accurate and interpretable scoring functions incorporating force field and empirical terms.

✨ Independent scoring functions for geometry minimization, pose selection and affinity ranking.

✨ Easy-to-use Python API and command-line tools.

### Work in progress

🚧 Affinity-ranking score checkpoint and screening power evaluation result.

🚧 Traninig code.

## Installation

### 1. Create a conda environment:

To minimize environment setup cost, it is recommended to create an Conda environment.

```bash
git clone https://github.com/bytedance/Protenix-Dock.git
cd protenix-dock

sudo apt-get update && sudo apt-get install -y libxrender1 libxext6
conda env create -f environment.yml
```

### 2. Install the Python package:

For better compatibility between packages, it is safe to install Protenix-Dock from
source.

```bash
python3 setup.py install
```

If your CPU is equiped with AVX2 instructions, you can build a faster one.
```bash
export PXDOCK_ENABLE_AVX2=1
python3 setup.py install
```

### 3. Install command-line tools (Optional):

If receptors & ligands are already prepared and only docking/optimizatioin/evaluation
is required, you can install command-lines tools from source.

```bash
pushd engine

mkdir build
cd build

destdir=~/pxdock
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$destdir \
    -DBDOCK_AVX2=OFF  # If your CPU supports AVX2, turn on it for better speed
make -j8 install

confdir=$destdir/conf
mkdir $confdir
cp ../../pxdock/data/pscore-v7_and_bscore-fake.json $confdir

popd
```

## Docking

### Receptor Preparation

PDB files downloaded directly from the Protein Data Bank typically lack hydrogen atoms and proper protonation state assignments. Protenix-Dock requires preprocessed receptor files where residues like histidine are resolved to their protonated forms (HID/HIE/HIP).

The simplest approach is to add hydrogens with `reduce` (included in AmberTools):

```bash
# Download structure from RCSB
curl -s "https://files.rcsb.org/download/5S8I.pdb.gz" | gunzip > raw.pdb

# Extract protein ATOM records only (remove water, ligands, etc.)
grep "^ATOM" raw.pdb > protein_only.pdb
echo "END" >> protein_only.pdb

# Add hydrogens (required for HIS protonation detection)
reduce -BUILD -Quiet protein_only.pdb | grep "^ATOM\|^END" > receptor_prepared.pdb
```

Alternatively, use `pdb2pqr` for more control over protonation:

```bash
pdb2pqr30 --ff=AMBER --ffout=AMBER --pdb-output receptor_prepared.pdb receptor.pdb receptor.pqr
```

See the [pdb2pqr documentation](https://pdb2pqr.readthedocs.io/) for more options (e.g. `--with-ph` for custom pH).

### Ligand Preparation

Extract the co-crystallized ligand from the PDB and convert to SDF format using RDKit:

```python
from rdkit import Chem

# Read the raw PDB and extract HETATM records for the target ligand
# Replace "ZXG" with your ligand's 3-letter residue code (from RCSB)
with open("raw.pdb") as f:
    lig_lines = [l for l in f if l.startswith("HETATM") and l[17:20].strip() == "ZXG"]
lig_lines.append("END\n")

mol = Chem.MolFromPDBBlock("".join(lig_lines), removeHs=True, sanitize=True)
writer = Chem.SDWriter("ligand_prepared.sdf")
writer.write(mol)
writer.close()
```

You can find the ligand residue code on the [RCSB PDB](https://www.rcsb.org/) page for your structure under "Small Molecules".

### Pocket Definition

Define the docking search box as a center point and dimensions (in Angstroms). There are three ways:

**Option A**: Compute automatically from a reference ligand pose:
```python
from pxdock.pipeline.get_pocket_info import compute_pocket_box
center, size = compute_pocket_box("ligand_prepared.sdf", buffer=5.0)
```

**Option B**: Use autobox (does this internally):
```python
dock_instance.autobox("ligand_prepared.sdf", buffer=5.0)
```

**Option C**: Write a pocket config file manually:
```
center_x = -21.702
center_y = 13.751
center_z = 27.892
size_x = 12.062
size_y = 9.510
size_z = 12.016
```

Load it with:
```python
from pxdock.dock.utils import parse_pocket_config
pocket = parse_pocket_config("pocket.config")
box_center = (pocket["center_x"], pocket["center_y"], pocket["center_z"])
box_size = (pocket["size_x"], pocket["size_y"], pocket["size_z"])
```

### Batch Fetch from RCSB

To download and prepare multiple systems at once, edit `scripts/fetch_examples.py` and add your PDB ID / ligand code pairs:

```python
SYSTEMS = [
    ("8AY3", "OE3"),
    ("7BMI", "U4B"),
    ("5SB2", "1K2"),
    ("7T0D", "FPP"),
]
```

Then run:
```bash
python scripts/fetch_examples.py
```

This downloads each PDB, extracts the protein (with hydrogens via `reduce`), extracts the ligand as SDF, and computes the pocket box. Output files are written to `examples/`.

### Usage

#### Run with Python (recommended):

**Autobox** (easiest &mdash; computes the search box from a reference ligand):
```python
from pxdock import ProtenixDock

dock = ProtenixDock("receptor_prepared.pdb")
dock.autobox("ligand_prepared.sdf", buffer=5.0)
out_dir = dock.run_docking("ligand_prepared.sdf", out_dir="output/", num_walker=4)
```

**Manual box definition**:
```python
from pxdock import ProtenixDock

dock = ProtenixDock("receptor_prepared.pdb")
dock.set_box(
    pocket_center=[-21.702, 13.751, 27.892],
    box_size=[12.062, 9.510, 12.016],
)
out_dir = dock.run_docking("ligand_prepared.sdf", out_dir="output/", num_walker=4)
```

**Full example from PDB ID** (end-to-end):
```python
from pxdock import ProtenixDock
from pxdock.dock.utils import parse_pocket_config

# Assuming files prepared as described above
pocket = parse_pocket_config("examples/5s8i_pocket.config")
box_center = (pocket["center_x"], pocket["center_y"], pocket["center_z"])
box_size = (pocket["size_x"], pocket["size_y"], pocket["size_z"])

dock = ProtenixDock("examples/5s8i_protein_prepared.pdb")
dock.set_box(box_center, box_size)
out_dir = dock.run_docking(
    "examples/5s8i_ligand_prepared.sdf",
    out_dir="output/",
    num_walker=4,
    seed=101,
)
# Docked poses are written as SDF files in output/
```

**Grid cache for multiple ligands** (generate once, reuse):
```python
dock = ProtenixDock("receptor_prepared.pdb")
dock.set_box(box_center, box_size)

# Generate cache maps (0.175 A spacing balances accuracy and speed)
dock.generate_cache_maps(spacing=0.175, output_dir="cache/")

# Dock multiple ligands against the same receptor
dock.run_docking("ligand_1.sdf", out_dir="out1/", num_walker=4)
dock.run_docking("ligand_2.sdf", out_dir="out2/", num_walker=4)

# Restore from disk in a new session:
# dock.restore_cache_maps("cache/")
```

**Virtual screening: N ligands against 1 receptor pocket**:
```python
from pxdock import screen_ligands

results = screen_ligands(
    receptor_pdb="receptor_prepared.pdb",
    pocket_center=[-21.702, 13.751, 27.892],
    box_size=[12.062, 9.510, 12.016],
    ligand_files=["ligand_1.sdf", "ligand_2.sdf", "ligand_3.sdf"],
    cache_spacing=0.175,
    num_walker=256,
    out_dir="screen_output/",
)
# results is a list of dicts: [{"ligand": "ligand_1.sdf", "out_dir": "screen_output/0_ligand_1"}, ...]
for entry in results:
    if entry["out_dir"]:
        print(f"{entry['ligand']} -> {entry['out_dir']}")
    else:
        print(f"{entry['ligand']} FAILED: {entry.get('error')}")
```

The receptor is loaded once and the grid cache is generated once, then shared across all ligands. Failed ligands are isolated and don't crash the batch.

**Ensemble/multi-pocket docking: 1 ligand against N pockets**:
```python
from pxdock import dock_against_pockets

results = dock_against_pockets(
    ligand_sdf="ligand_prepared.sdf",
    pockets=[
        {
            "receptor": "receptor_conf1.pdb",
            "center": [-21.7, 13.8, 27.9],
            "box_size": [12.0, 9.5, 12.0],
            "name": "conf1_active_site",
        },
        {
            "receptor": "receptor_conf2.pdb",
            "center": [-18.0, 15.0, 30.0],
            "box_size": [14.0, 14.0, 14.0],
            "name": "conf2_allosteric",
        },
    ],
    cache_spacing=0.175,
    num_walker=256,
    out_dir="multipocket_output/",
)
# results is a dict: {"conf1_active_site": "multipocket_output/conf1_active_site", ...}
for name, out_dir in results.items():
    print(f"{name} -> {out_dir}")
```

Pockets can also be passed as tuples: `(receptor_pdb, center, box_size)`. Each pocket gets a fresh `ProtenixDock` instance with its own receptor and cache.

**Score an existing complex** (no docking, just scoring):
```python
dock = ProtenixDock("receptor_prepared.pdb")
dock.autobox("docked_ligand.sdf")
scores = dock.score_complex("docked_ligand.sdf")
for s in scores:
    print(f"Pose {s['pose_id']}: pscore={s['pscore']:.3f}")
```

**Pose optimization** (Vina-seeded refinement):
```python
dock = ProtenixDock("receptor_prepared.pdb")
dock.set_box(box_center, box_size)
out_dir = dock.run_pose_opt("ligand_prepared.sdf", out_dir="output/")
```

### Run tests

```bash
# In these tests, we set the spacing to 0.5 in order to quickly complete the functional test.
cd test
# performing preare ligand, receptor and docking separately.
python3 test_data_prepare.py

# run docking or pose_opt by `ProtenixDock` class.
python3 test_protenix_dock.py

# batch docking tests (screen_ligands + dock_against_pockets).
python3 test_batch_docking.py

# calculate pose rmsd.
python3 test_rmsd.py
```

### End-to-end tests (optional)

The e2e test exercises the full pipeline: raw PDB → pdb2pqr → receptor prep → docking → scoring. It requires **AmberTools** (`tleap`, `antechamber`, `parmchk2`, `pdb4amber`), which is only installable via conda/micromamba.

#### Install AmberTools

```bash
# Using conda:
conda install -c conda-forge ambertools

# Or using micromamba:
micromamba install -c conda-forge ambertools
```

#### Run e2e tests

```bash
pytest test/test_e2e.py -m e2e -v
```

If AmberTools or pdb2pqr are not installed, the tests will auto-skip with a descriptive message.

## Contribution

Please check [Contributing](CONTRIBUTING.md) for more details. If you encounter problems using Protenix—Dock, feel free to create an issue! We also welcome pull requests from the community.

## Code of Conduct

Please check [Code of Conduct](CODE_OF_CONDUCT.md) for more details.

## Security

If you discover a potential security issue in this project, or think you may
have discovered a security issue, we ask that you notify Bytedance Security via our [security center](https://security.bytedance.com/src) or [vulnerability reporting email](sec@bytedance.com).

Please do **not** create a public GitHub issue.

## License 

The Protenix-Dock project is made available under the [GPLv3 License](./LICENSE)

Portions of the source code are based on the [Meeko](https://github.com/forlilab/Meeko) and [posebusters](https://github.com/maabuu/posebusters) project.

Portions of the SMARTS patterns used in Protenix-Dock are derived from the [ProLIF](https://github.com/chemosim-lab/ProLIF) and [OpenFF](https://github.com/openforcefield/openff-forcefields) project. 

## Contact

We welcome inquiries and collaboration opportunities for advanced applications of our framework, such as developing new features and more. Please feel free to contact us at ai4s-bio@bytedance.com.