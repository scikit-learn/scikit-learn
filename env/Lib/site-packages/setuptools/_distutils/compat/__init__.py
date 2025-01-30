from __future__ import annotations


def consolidate_linker_args(args: list[str]) -> list[str] | str:
    """
    Ensure the return value is a string for backward compatibility.

    Retain until at least 2025-04-31. See pypa/distutils#246
    """

    if not all(arg.startswith('-Wl,') for arg in args):
        return args
    return '-Wl,' + ','.join(arg.removeprefix('-Wl,') for arg in args)
