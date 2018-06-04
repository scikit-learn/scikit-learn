""" Dataframe optimizations """
from __future__ import absolute_import, division, print_function

from ..optimization import cull, fuse_getitem, fuse
from ..context import _globals
from .. import core

try:
    import fastparquet  # noqa: F401
except ImportError:
    fastparquet = False


def optimize(dsk, keys, **kwargs):
    from .io import dataframe_from_ctable
    if isinstance(keys, list):
        dsk, dependencies = cull(dsk, list(core.flatten(keys)))
    else:
        dsk, dependencies = cull(dsk, [keys])
    dsk = fuse_getitem(dsk, dataframe_from_ctable, 3)
    if fastparquet:
        from .io.parquet import _read_parquet_row_group
        dsk = fuse_getitem(dsk, _read_parquet_row_group, 4)
    dsk, dependencies = fuse(dsk, keys, dependencies=dependencies,
                             ave_width=_globals.get('fuse_ave_width', 0))
    dsk, _ = cull(dsk, keys)
    return dsk
