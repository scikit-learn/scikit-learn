from __future__ import annotations

import contextlib
from io import StringIO
from typing import TYPE_CHECKING, Any

from polars.io.csv.functions import read_csv

with contextlib.suppress(ImportError):
    from polars.polars import read_clipboard_string as _read_clipboard_string

if TYPE_CHECKING:
    from polars import DataFrame


def read_clipboard(separator: str = "\t", **kwargs: Any) -> DataFrame:
    """
    Read text from clipboard and pass to `read_csv`.

    Useful for reading data copied from Excel or other similar spreadsheet software.

    Parameters
    ----------
    separator
        Single byte character to use as separator parsing csv from clipboard.
    kwargs
        Additional arguments passed to `read_csv`.

    See Also
    --------
    read_csv : Read a csv file into a DataFrame.
    DataFrame.write_clipboard : Write a DataFrame to the clipboard.
    """
    csv_string: str = _read_clipboard_string()
    io_string = StringIO(csv_string)
    return read_csv(source=io_string, separator=separator, **kwargs)
