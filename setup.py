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

import multiprocessing
import os
import shutil
import subprocess

from setuptools import Extension, find_namespace_packages, setup
from setuptools.command.build_ext import build_ext

package_name = "pxdock"


class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_INSTALL_PREFIX=" + extdir,
            "-DBDOCK_PYSDK=ON"
        ]
        build_args = ["--config", "Release"]
        num_jobs = multiprocessing.cpu_count()
        build_args += ["--", f"-j{num_jobs}"]

        os.makedirs(self.build_temp, exist_ok=True)
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", ".", "--target", "install"] + build_args,
            cwd=self.build_temp,
        )
        self.copy_cmake_output(extdir)

    def copy_cmake_output(self, target_dir):
        dst_dir = os.path.join(target_dir, package_name)
        src_dir = os.path.join(target_dir, "lib")
        if os.path.exists(src_dir):
            for filename in os.listdir(src_dir):
                full_file_name = os.path.join(src_dir, filename)
                if os.path.isfile(full_file_name):
                    self.copy_file(full_file_name, dst_dir)
            shutil.rmtree(src_dir)


class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])
        sourcedir = os.path.join(os.path.dirname(__file__), "engine")
        self.sourcedir = os.path.abspath(sourcedir)


with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name=package_name,
    version="0.0.1",
    python_requires=">=3.9",
    keywords=["docking"],
    description="A high-precision Docking engine.",
    license="Apache License 2.0",
    url="",
    packages=find_namespace_packages(exclude=["test", "test.*", "examples"]),
    package_data={package_name: ["data/*", "data/amberff/*", "data/templates/*"]},
    include_package_data=True,
    ext_modules=[CMakeExtension(package_name)],
    cmdclass=dict(build_ext=CMakeBuild),
    platforms="manylinux1",
    install_requires=install_requires,
)
