# Protenix-Dock

This repository hosts the source code for our work "Protenix-Dock: An accurate and trainable end-to-end protein-ligand docking framework using empirical scoring functions".  

For more information about the implementation and the performance of Protenix-Dock, see our [technical report](ProtenixDock_Technical_Report.pdf).

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

### 3. Install command-line tools (Optional):

If receptors & ligands are already prepared and only docking/optimizatioin/evaluation
is required, you can install command-lines tools from source.

```bash
pushd engine

mkdir build
cd build

destdir=~/pxdock
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$destdir
make -j8 install

confdir=$destdir/conf
mkdir $confdir
cp ../../pxdock/data/pscore-v7_and_bscore-fake.json $confdir

popd
```

## Docking

### Usage


#### Run with Python (recommended):

#### Run with Python (recommended):

* Autobox (easy benchmark)
```python
from pxdock import ProtenixDock
receptor_pdb = "path/to/receptor.pdb"
true_ligand_sdf = "path/to/true_ligand.sdf"
ligand_sdf = "path/to/ligand.sdf"

dock_instance = ProtenixDock(receptor_pdb)
dock_instance.autobox(true_ligand_sdf)
out_dir = dock_instance.run_docking(ligand_sdf)
```

* Manual input: 
```python
from pxdock import ProtenixDock
receptor_pdb = "path/to/receptor.pdb"
ligand_sdf = "path/to/ligand.sdf"

box_center = [0., 0., 0.] # box center for receptor
box_size = [10., 10., 10.] # box size for receptor
dock_instance = ProtenixDock(receptor_pdb)
dock_instance.set_box(box_center, box_size)

# Optional: you can generate cache maps for receptor, and then you can load it for next docking.
# out_dir = dock_instance.generate_cache_maps(spacing=0.5)
# in next run: 
# dock_engine.load_cache_maps(out_dir)

# the docking_res_files is in json format.
docking_res_files = dock_instance.run_docking(ligand_sdf)
```

### Run tests

```bash
cd test
# performing preare ligand, receptor and docking separately.
python3 test_data_prepare.py

# run docking or pose_opt by `ProtenixDock` class.
python3 test_protenix_dock.py

# calculate pose rmsd.
python3 test_rmsd.py
```


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