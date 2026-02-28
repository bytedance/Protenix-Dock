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
import platform
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
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBDOCK_PYSDK=ON"
        ]
        # On macOS with conda, help CMake find the right Boost and Python
        import sys
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        if conda_prefix:
            cmake_args += [f"-DCMAKE_PREFIX_PATH={conda_prefix}"]
            cmake_args += [f"-DCMAKE_CXX_FLAGS=-DBOOST_LOG_DYN_LINK -DBOOST_ALL_DYN_LINK"]
            # Help find boost_python with version suffix
            py_ver = f"{sys.version_info.major}{sys.version_info.minor}"
            bp_lib = os.path.join(conda_prefix, "lib", f"libboost_python{py_ver}.dylib")
            if os.path.exists(bp_lib):
                cmake_args += [f"-DBoost_PYTHON_LIBRARY={bp_lib}"]
            # Use conda compiler if available
            cxx = os.path.join(conda_prefix, "bin", "clang++")
            if os.path.exists(cxx):
                cmake_args += [f"-DCMAKE_CXX_COMPILER={cxx}"]
        if int(os.environ.get("PXDOCK_ENABLE_AVX2", "0")) > 0:
            cmake_args += ["-DBDOCK_AVX2=ON"]
        elif int(os.environ.get("PXDOCK_ENABLE_NEON", "0")) > 0 or \
             (platform.machine() == "arm64" and
              os.environ.get("PXDOCK_DISABLE_SIMD") is None):
            cmake_args += ["-DBDOCK_NEON=ON"]
        build_args = ["--config", "Release"]
        num_jobs = multiprocessing.cpu_count()
        build_args += ["--", f"-j{num_jobs}"]

        enable_pgo = int(os.environ.get("PXDOCK_ENABLE_PGO", "0")) > 0
        pgo_dir = os.path.join(self.build_temp, "pgo_profiles")

        os.makedirs(self.build_temp, exist_ok=True)

        if enable_pgo:
            # Detect compiler for PGO flag differences
            is_clang = (platform.system() == "Darwin") or shutil.which("clang++")
            boost_defs = " -DBOOST_LOG_DYN_LINK -DBOOST_ALL_DYN_LINK" if conda_prefix else ""

            if is_clang:
                gen_flag = f"-fprofile-instr-generate={pgo_dir}/default.profraw"
                profdata = os.path.join(pgo_dir, "default.profdata")
                use_flag = f"-fprofile-instr-use={profdata}"
            else:
                gen_flag = f"-fprofile-generate={pgo_dir}"
                use_flag = f"-fprofile-use={pgo_dir}"

            # Pass 1: build with profile generation
            pgo_gen_args = cmake_args + [
                f"-DCMAKE_CXX_FLAGS_RELEASE=-O3 -ffast-math {gen_flag}{boost_defs}",
            ]
            pgo_gen_args = [a for a in pgo_gen_args if not a.startswith("-DCMAKE_CXX_FLAGS=")]
            os.makedirs(pgo_dir, exist_ok=True)
            subprocess.check_call(
                ["cmake", ext.sourcedir] + pgo_gen_args, cwd=self.build_temp
            )
            subprocess.check_call(
                ["cmake", "--build", ".", "--target", "install"] + build_args,
                cwd=self.build_temp,
            )
            self.copy_cmake_output(extdir)

            # Run training workload to generate profile data
            print("PGO: Running training workload...")
            training_script = os.path.join(
                os.path.dirname(__file__), "test", "bench_engine.py"
            )
            if os.path.exists(training_script):
                subprocess.call([sys.executable, training_script],
                                env={**os.environ, "PXDOCK_PGO_TRAINING": "1"})
            else:
                print("PGO: Warning — bench_engine.py not found, skipping training")

            # For Clang: merge profraw files into profdata
            if is_clang:
                import glob as globmod
                profraw_files = globmod.glob(os.path.join(pgo_dir, "*.profraw"))
                if profraw_files:
                    llvm_profdata = shutil.which("llvm-profdata") or shutil.which("xcrun")
                    if llvm_profdata and "xcrun" in llvm_profdata:
                        merge_cmd = ["xcrun", "llvm-profdata", "merge", "-output", profdata] + profraw_files
                    elif llvm_profdata:
                        merge_cmd = [llvm_profdata, "merge", "-output", profdata] + profraw_files
                    else:
                        merge_cmd = ["llvm-profdata", "merge", "-output", profdata] + profraw_files
                    print(f"PGO: Merging {len(profraw_files)} profraw files...")
                    subprocess.check_call(merge_cmd)
                else:
                    print("PGO: Warning — no .profraw files found, PGO pass 2 may not help")

            # Pass 2: rebuild with profile data
            for f in os.listdir(self.build_temp):
                if f != "pgo_profiles":
                    p = os.path.join(self.build_temp, f)
                    if os.path.isdir(p):
                        shutil.rmtree(p)
                    else:
                        os.remove(p)
            pgo_use_args = cmake_args + [
                f"-DCMAKE_CXX_FLAGS_RELEASE=-O3 -ffast-math {use_flag}{boost_defs}",
            ]
            pgo_use_args = [a for a in pgo_use_args if not a.startswith("-DCMAKE_CXX_FLAGS=")]
            subprocess.check_call(
                ["cmake", ext.sourcedir] + pgo_use_args, cwd=self.build_temp
            )
            subprocess.check_call(
                ["cmake", "--build", ".", "--target", "install"] + build_args,
                cwd=self.build_temp,
            )
            self.copy_cmake_output(extdir)
        else:
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
        for name in ("lib", "lib64"):
            src_dir = os.path.join(target_dir, name)
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
