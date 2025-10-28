from __future__ import annotations

import contextlib

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr


def set_random_seed(seed: int) -> None:
    r"""
    Set the global random seed for Polars.

    This random seed is used to determine things such as shuffle ordering.


    Parameters
    ----------
    seed
        A non-negative integer < 2\ :sup:`64` used to seed the internal global
        random number generator.
    """
    plr.set_random_seed(seed)
