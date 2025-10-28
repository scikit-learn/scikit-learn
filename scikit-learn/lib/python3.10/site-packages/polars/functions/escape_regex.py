from __future__ import annotations

import contextlib

from polars._utils.various import qualified_type_name

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr
import polars._reexport as pl


def escape_regex(s: str) -> str:
    r"""
    Escapes string regex meta characters.

    Parameters
    ----------
    s
        The string whose meta characters will be escaped.

    """
    if isinstance(s, pl.Expr):
        msg = "escape_regex function is unsupported for `Expr`, you may want use `Expr.str.escape_regex` instead"
        raise TypeError(msg)
    elif not isinstance(s, str):
        msg = f"escape_regex function supports only `str` type, got `{qualified_type_name(s)}`"
        raise TypeError(msg)

    return plr.escape_regex(s)
