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

import argparse
import json
import os

from pxdock.common import get_logger
from pxdock.parser.cofactor import CofactorParser
from pxdock.parser.receptor import ReceptorParser

logger = get_logger(__name__)


def prepare_receptor(
    receptor_fpath: str,
    box_center: list,
    box_size: list,
    output_fpath: str = None,
    cofactor_fpath: str = None,
    verbose=True,
    find_pocket=True,
):
    receptor_fpath = os.path.abspath(receptor_fpath)
    try:
        receptor = ReceptorParser(
            pdb_file=receptor_fpath,
            center=box_center,
            box=box_size,
            find_pocket=find_pocket,
        )
        if cofactor_fpath is not None:
            if cofactor_fpath.endswith(".sdf") or cofactor_fpath.endswith(".pdb"):
                cofactor_files = [cofactor_fpath]
            else:
                cofactor_files = os.listdir(cofactor_fpath)
                cofactor_files = [
                    os.path.join(cofactor_fpath, cofactor_file)
                    for cofactor_file in cofactor_files
                    if cofactor_file.endswith(".sdf") or cofactor_file.endswith(".pdb")
                ]
            assert len(cofactor_files) > 0
            for cofactor_file in cofactor_files:
                cofactor = CofactorParser(ligand=cofactor_file)
                receptor.add_cofactor(cofactor.get_data())
        receptor_data = receptor.get_data()
        if output_fpath is None:
            output_fpath = os.path.join(
                os.path.dirname(receptor_fpath),
                f"{os.path.splitext(os.path.basename(receptor_fpath))[0]}-prepared-receptor.json",
            )
        else:
            output_fpath = os.path.abspath(output_fpath)
        os.makedirs(os.path.dirname(output_fpath), exist_ok=True)

        with open(output_fpath, "w") as f:
            json.dump(receptor_data, f)
        if verbose:
            logger.info(f"Saved receptor_data to {output_fpath}")
        return output_fpath
    except Exception as e:
        logger.warning(
            f"Error in parsing receptor {receptor_fpath}, cofactor {cofactor_fpath}: {e} "
        )
        return None


if __name__ == "__main__":

    _parser = argparse.ArgumentParser()
    _parser.add_argument("--receptor_fpath", type=str, required=True)
    _parser.add_argument("--output_fpath", type=str, required=True)
    _parser.add_argument("--box_center", type=str, default=None)
    _parser.add_argument("--box_size", type=str, default=None)
    _parser.add_argument("--cofactor_fpath", type=str, default=None)

    args = _parser.parse_args()
    print(args)

    args.box_center = list(map(float, args.box_center.split(",")))
    args.box_size = list(map(float, args.box_size.split(",")))

    print("Got box center:")
    print(args.box_center)
    print("Got box size:")
    print(args.box_size)

    prepare_receptor(
        receptor_fpath=args.receptor_fpath,
        output_fpath=args.output_fpath,
        box_center=args.box_center,
        box_size=args.box_size,
        cofactor_fpath=args.cofactor_fpath,
    )
