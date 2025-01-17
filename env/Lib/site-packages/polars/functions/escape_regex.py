from __future__ import annotations

import contextlib

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr
import polars._reexport as pl


def escape_regex(s: str) -> str:
    r"""
    Escapes string regex meta characters.

    Parameters
    ----------
    s
        The string that all of its meta characters will be escaped.

    """
    if isinstance(s, pl.Expr):
        msg = "escape_regex function is unsupported for `Expr`, you may want use `Expr.str.escape_regex` instead"
        raise TypeError(msg)
    elif not isinstance(s, str):
        msg = f"escape_regex function supports only `str` type, got `{type(s)}`"
        raise TypeError(msg)

    return plr.escape_regex(s)
