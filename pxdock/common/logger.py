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


class WarningFilter(logging.Filter):
    def filter(self, record):
        # Filter out specific warnings by their message or other criteria
        if "simtk.openmm" in record.getMessage():
            return False
        return True


def get_logger(name="", loglevel="INFO", log_file_path=None):
    root_logger = logging.getLogger()
    logger = logging.getLogger(name)
    # we only add handlers to the root logger! Let the propogation handle the rest.
    add_handlers(root_logger, loglevel, log_file_path)
    logging.getLogger("numexpr").setLevel(logging.WARNING)
    logging.getLogger("MDAnalysis.coordinates.AMBER").setLevel(logging.ERROR)
    return logger


def add_handlers(logger, loglevel, log_file_path=None):
    fmt = "%(asctime)-15s [%(pathname)s:%(lineno)d] %(levelname)s %(name)s: %(message)s"
    formatter = logging.Formatter(fmt)
    loglevel = getattr(logging, loglevel.upper(), logging.INFO)
    logger.setLevel(loglevel)

    if not logger.handlers:
        handler = logging.StreamHandler()
        logger.addHandler(handler)
    else:
        handler = logger.handlers[0]
    handler.setFormatter(formatter)
    warning_filter = WarningFilter()
    handler.addFilter(warning_filter)

    # we output to at most two streams: one stdout and one file
    if log_file_path is not None and len(logger.handlers) == 1:
        handler = logging.FileHandler(log_file_path, mode="a")
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        handler.addFilter(warning_filter)

    return logger
