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
import subprocess

from pxdock.common import get_logger
from pxdock.vina.protein import PrepProt

try:
    import AutoDockTools as ADT

    ADTPath = ADT.__path__[0]
except ImportError:
    raise f"import AutoDockTools error, please install `AutoDockTools`"

logger = get_logger(__name__)


class PrepProt_py39(PrepProt):
    def get_pdbqt(self, prot_pdbqt):
        if self.prot_pqr is None:
            self.prot_pqr = self.prot
        if not os.path.exists(prot_pdbqt):
            prepare_receptor = f"{ADTPath}/Utilities24/prepare_receptor4.py"
            subprocess.Popen(
                ["python3", prepare_receptor, "-r", self.prot_pdb, "-o", prot_pdbqt]
            ).communicate()
        if not os.path.exists(prot_pdbqt):
            raise ValueError(
                f"pdbqt path: {prot_pdbqt} does not exists! Some thing must be wrong during the preparation. Please check the log."
            )

    def removeH(self, prot_pdb_noH):
        with open(self.prot, "r") as infile:
            with open(prot_pdb_noH, "w") as outfile:
                for line in infile:
                    # Checking if it's an ATOM line and the element is not 'H'
                    if line.startswith("ATOM") and line[76:78].strip() == "H":
                        continue
                    outfile.write(line)
        self.prot = prot_pdb_noH


def prepare_receptor_pdbqt(receptor_pdb, save_dir, del_metal=True):

    pdb_name = os.path.basename(receptor_pdb)
    assert pdb_name.endswith(".pdb")
    pdb_name = pdb_name[: -len(".pdb")]

    prot_dry = os.path.join(save_dir, f"{pdb_name}_dry.pdb")
    prot_metal = os.path.join(save_dir, f"{pdb_name}_metal.pdb")
    prot_pdb = os.path.join(save_dir, f"{pdb_name}_temp.pdb")
    prot_pqr = os.path.join(save_dir, f"{pdb_name}.pqr")
    prot_pdbqt = os.path.join(save_dir, f"{pdb_name}.pdbqt")
    prot_pdb_noH = os.path.join(save_dir, f"{pdb_name}_noH.pdb")

    if os.path.exists(prot_pdbqt):
        logger.info(
            f"Will skip processing the protein since an existing pdbqt file is found: {prot_pdbqt}"
        )
        return prot_pdbqt

    b = PrepProt_py39(pdb_file=receptor_pdb)
    b.del_water(prot_dry)
    if del_metal:
        b.del_metal(prot_metal)
    try:
        b.addH(prot_pqr, prot_pdb)  # call pdb2pqr30 in it
        b.get_pdbqt(prot_pdbqt)
    except Exception:
        logger.info(f"remove Hs as {prot_pdb_noH} and retry prepare_receptor")
        b.removeH(prot_pdb_noH)
        b.addH(prot_pqr, prot_pdb)
        b.get_pdbqt(prot_pdbqt)
    return prot_pdbqt
