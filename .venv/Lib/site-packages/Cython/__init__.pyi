from typing import Any

from .Shadow import *


# Internal interface for IPython. Not further typed since it's not meant for end users.
def load_ipython_extension(ip: Any) -> None: ...
