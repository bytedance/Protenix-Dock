import os

from .logger import get_logger
from .utilities import my_random_string

kWorkDir = f"/opt/tiger/yodel/container/tmp/protenix-dock/{my_random_string(10)}"
os.makedirs(kWorkDir, exist_ok=True)

__all__ = ["get_logger", "my_random_string", "kWorkDir"]
