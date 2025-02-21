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

import pxdock
from pxdock.common import get_logger, kWorkDir, my_random_string

logger = get_logger(__name__)

protenix_dock_dir = os.path.dirname(pxdock.__file__)
kScoreConfigFile = os.path.join(protenix_dock_dir, "data/pscore-v7_and_bscore-fake.json")


def parse_pocket_config(pocket_config: str) -> dict:
    with open(pocket_config, "r") as file:
        lines = file.readlines()

    config = {}
    for line in lines:
        if line.strip():
            key, value = line.split("=")
            config[key.strip()] = float(value.strip())

    return config


def write_pocket_config(pocket_center: list[int] | tuple[int], box_size: list[int] | tuple[int], out_fpath: str | None = None) -> str:
    if out_fpath is None:
        out_fpath = os.path.join(kWorkDir, f"{my_random_string()}_pocket.config")

    with open(out_fpath, "w") as f:
        f.write("center_x=%f\n" % pocket_center[0])
        f.write("center_y=%f\n" % pocket_center[1])
        f.write("center_z=%f\n" % pocket_center[2])
        f.write("size_x=%f\n" % box_size[0])
        f.write("size_y=%f\n" % box_size[1])
        f.write("size_z=%f\n" % box_size[2])
    return out_fpath
