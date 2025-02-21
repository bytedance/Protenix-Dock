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

import copy
import dataclasses
import glob
import json
import multiprocessing
import os
import random
import subprocess
import threading
import traceback
import uuid
from collections.abc import Mapping
from typing import Dict, Mapping, Tuple

import jinja2
import numpy as np
import pandas as pd
from func_timeout import exceptions as func_timeout_exceptions
from func_timeout import func_set_timeout
from rdkit import Chem

from .logger import get_logger

logger = get_logger(__name__)


def seed_everything(seed, deterministic=False):
    random.seed(seed)
    np.random.seed(seed)


def convert_value_to_list(obj):
    if dataclasses.is_dataclass(obj):
        return convert_value_to_list(dataclasses.asdict(obj))
    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    if isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    if isinstance(obj, (list, tuple)):
        return [convert_value_to_list(item) for item in obj]
    if isinstance(obj, Mapping):
        return {str(key): convert_value_to_list(obj[key]) for key in obj}
    if isinstance(obj, np.ndarray):
        return convert_value_to_list(obj.tolist())
    return obj


def padding(lst, idx=0):
    if isinstance(lst, list) and len(lst) > 0 and isinstance(lst[0], list):
        max_len = max(map(len, lst))
        padded = [l + [idx] * (max_len - len(l)) for l in lst]
        return padded
    else:
        return lst



def my_random_string(string_length=10):
    """Returns a random string of length string_length."""
    random = str(uuid.uuid4())  # Convert UUID format to a Python string.
    random = random.upper()  # Make all characters uppercase.
    random = random.replace("-", "")  # Remove the UUID '-'.
    return random[0:string_length]  # Return the random string


def write_to_txt(lists, files):
    with open(files, "w") as f:
        lines = [str(line) + "\n" for line in lists]
        f.writelines(lines)


def read_txt(files):
    with open(files, "r") as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
    return lines


def sdf_string_iterator(input_fpath):
    molecule = []
    with open(input_fpath, "r") as sdf_file:
        for line in sdf_file:
            molecule.append(line)
            if line.startswith("$$$$"):
                yield "".join(molecule)
                molecule = []


def save_sdf_supplier(fpath, skip_failed_mol=True, removeHs=False):
    with open(fpath, "r") as f:
        sdf_strings = f.read()

    ms = Chem.SDMolSupplier()
    ms.SetData(sdf_strings, removeHs=removeHs)

    for j, mol in enumerate(ms):
        if mol is None:
            logger.warning(f"[read sdf exception] [{j}-th mol]: None")
            if skip_failed_mol:
                continue
        yield mol


def sdf_batch_iterator(input_fpath, batch_size):
    batch = []
    for j, sdf_string in enumerate(sdf_string_iterator(input_fpath)):
        batch.append(sdf_string)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def dataframe_batch_iterator(input_fpath, batch_size):
    df = pd.read_csv(input_fpath, index_col=False)
    # split csv into batch_size
    groups = df.groupby(df.index // batch_size)
    for i, group in groups:
        yield group


def modify_molecule_weights(m):
    """
    Cap atomic masses such that all atoms with mass greater than chlorine
    are assigned the mass of chlorine.
    """
    m = Chem.Mol(m)
    total_cut = 0.0
    for atom in m.GetAtoms():
        if atom.GetMass() > 35.45:  # 35.45 is the mass of chlorine atom
            total_cut += 35.45 - atom.GetMass()
    return total_cut



def load_json(x):
    if x is None:
        return None
    if (isinstance(x, float) or isinstance(x, np.ndarray)) and np.isnan(x):
        return None
    if not isinstance(x, str):
        return x
    try:
        return json.loads(x)
    except (json.JSONDecodeError, TypeError):
        return None


def get_file_paths(directory_or_pattern):
    """
    Retrieve file paths from a specified directory or based on a file pattern.

    If the provided input is a directory, the function will return paths to all files within this directory.
    If the input is a file pattern (e.g., '*.txt'), it will return paths matching the pattern.

    Args:
    directory_or_pattern (str): A directory path or a file pattern.

    Returns:
    list: A list of file paths.
    """
    if os.path.exists(directory_or_pattern) and os.path.isdir(directory_or_pattern):
        # List all files in the specified directory
        return [
            os.path.join(directory_or_pattern, filename)
            for filename in os.listdir(directory_or_pattern)
        ]
    else:
        # Assume it's a pattern and use glob to find matching files
        return glob.glob(directory_or_pattern)


def run_command_and_check(
    cmd: str,
    *,
    allow_error: bool = False,
    separate_stderr: bool = True,
    env: Dict = None,
    redirect_stdout=True,
    user_input: str = None,
    timeout: float = None,
) -> Tuple[int, str, str]:

    cmdenv = dict()
    cmdenv.update(**os.environ)
    if env is not None:
        cmdenv.update(**env)

    logger.debug(f"running cmd {cmd} with extra env {env}")
    stdout = subprocess.PIPE if redirect_stdout else None
    stderr = subprocess.PIPE if separate_stderr else subprocess.STDOUT
    stderr = stderr if redirect_stdout else None
    user_input = user_input.encode() if isinstance(user_input, str) else None
    result = subprocess.run(
        cmd,
        stdout=stdout,
        stderr=stderr,
        env=cmdenv,
        shell=True,
        check=False,
        input=user_input,
        timeout=timeout,
    )

    stdout = result.stdout.decode("utf-8") if result.stdout is not None else None
    stderr = result.stderr.decode("utf-8") if result.stderr is not None else None

    if result.returncode != 0 and not allow_error:
        logger.error(f"cmd {cmd}")
        logger.error(f"return code {result.returncode}")
        logger.error(f"stdout {stdout}")
        logger.error(f"stderr {stderr}")
        raise RuntimeError(f'fail to run "{cmd}". ')

    return result.returncode, stdout, stderr


def run_cmd(
    working_dir, cmd, name, output_dir, logger, user_input=None, out_f=None, err_f=None
):
    if isinstance(cmd, list):
        cmd_tokens = cmd
    else:
        cmd_tokens = cmd.split()
    cmd_tokens = list(map(str, cmd_tokens))
    cmd_bin = os.path.basename(cmd_tokens[0])
    out_f_none = False
    if out_f is None:
        out_f_none = True
        out_f = open(
            os.path.join(working_dir, output_dir, f"{cmd_bin}_stdout.log"), "a"
        )
    err_f_none = False
    if err_f is None:
        err_f_none = True
        err_f = open(
            os.path.join(working_dir, output_dir, f"{cmd_bin}_stderr.log"), "a"
        )
    logger.info(
        f"for {name}, output of {cmd_tokens} written to {out_f.name} and {err_f.name}, "
        f"working dir is {os.path.abspath(working_dir)}"
    )
    out_f.write(
        "\nbytedock logging:\n"
        f"for {name}, output of {cmd_tokens}\n"
        f"written to {out_f.name} and {err_f.name},\n"
        f"working dir is {os.path.abspath(working_dir)}\n"
    )
    out_f.flush()
    ret = 0
    try:
        if user_input is not None:
            subprocess.run(
                cmd_tokens,
                cwd=working_dir,
                stdout=out_f,
                stderr=err_f,
                text=False,
                input=user_input.encode(),
                check=True,
            )
        else:
            subprocess.run(
                cmd_tokens,
                cwd=working_dir,
                stdout=out_f,
                stderr=err_f,
                text=False,
                check=True,
            )
    except Exception as grepexc:
        logger.error(
            f"Error during subprocess {cmd}, please check the logs in {output_dir}, "
            f"error info: {grepexc.returncode} {grepexc.output}"
        )
        ret = 1
    finally:
        if out_f_none:
            out_f.close()
        if err_f_none:
            err_f.close()
        return ret


class TimeoutError(Exception):
    pass


def run_with_timeout(func, args=(), kwargs={}, timeout_duration=10):

    @func_set_timeout(timeout_duration)
    def run_func(*args, **kwargs):
        return func(*args, **kwargs)

    try:
        output = run_func(*args, **kwargs)
    except func_timeout_exceptions.FunctionTimedOut as exc:
        tb = traceback.format_exc()
        logger.warning(f"[ERROR fun]: {str(exc)} \n {tb}")
        return None
    return output


def run_func_with_timeout(func, args=(), kwargs={}, timeout_duration=10):
    result_queue = multiprocessing.Queue()

    def wrapper(result_queue):
        try:
            result_queue.put(func(*args, **kwargs))
        except Exception as e:
            # If an exception occurs, put the exception and its traceback in the queue
            result_queue.put((e, traceback.format_exc()))

    process = multiprocessing.Process(target=wrapper, args=(result_queue,))
    process.start()
    process.join(timeout_duration)

    if process.is_alive():
        process.terminate()
        process.join()
        raise TimeoutError(f"Function timed out after {timeout_duration} seconds")

    result = result_queue.get()
    if (
        isinstance(result, tuple)
        and len(result) == 2
        and isinstance(result[0], Exception)
    ):
        # If the result is a tuple containing an exception, re-raise it
        exception, traceback_str = result
        raise RuntimeError(f"Exception in func: {traceback_str}") from exception
    return result


def run_with_timeout_thread(func, args=(), kwargs={}, timeout_duration=10):
    # Create an inner function to execute the target function
    def target_func():
        nonlocal result
        result = func(*args, **kwargs)

    # Initialize the result to None
    result = None
    # Create and start a thread
    thread = threading.Thread(target=target_func)
    thread.start()
    # Wait for the thread to finish or the timeout to expire
    thread.join(timeout_duration)
    if thread.is_alive():
        # If the thread is still running after the timeout, raise a TimeoutError
        raise TimeoutError("Function timed out")
    else:
        # Otherwise, return the result of the function call
        return result


class TemplateRender(object):
    def __init__(self, template_dirs, template_name):
        """init the template render"""
        self.template_dirs = template_dirs
        self.template_name = template_name
        return

    def create_template(self):
        """create tempalte object"""
        loader = jinja2.FileSystemLoader(self.template_dirs)
        env = jinja2.Environment(loader=loader)
        template = env.get_template(self.template_name)
        return template

    def rend_template(self, rend_dict, output_dir, filename):
        """rend the template and generate the file"""
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        template = self.create_template()
        output = template.render(rend_dict)
        with open(output_dir + "/" + filename, "w") as file:
            file.write(output)
        return True
