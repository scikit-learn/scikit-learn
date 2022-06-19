"""Helper functions to get location of UNU.RAN source files."""

import pathlib
from typing import Union


def _unuran_dir(ret_path: bool = False) -> Union[pathlib.Path, str]:
    """Directory where root unuran/ directory lives."""
    p = pathlib.Path(__file__).parent / "unuran"
    return p if ret_path else str(p)
