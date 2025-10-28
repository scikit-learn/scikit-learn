# SPDX-FileCopyrightText: 2021 Filipe Laíns <lains@riseup.net>
# SPDX-FileCopyrightText: 2021 Quansight, LLC
# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import functools
import importlib.resources
import os
import sys
import typing


if sys.version_info >= (3, 9):
    from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
else:
    from typing import Collection, Iterable, Iterator, Mapping, Sequence


if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    cached_property = lambda x: property(functools.lru_cache(maxsize=None)(x))  # noqa: E731


if sys.version_info >= (3, 9):
    def read_binary(package: str, resource: str) -> bytes:
        return importlib.resources.files(package).joinpath(resource).read_bytes()
else:
    read_binary = importlib.resources.read_binary


if typing.TYPE_CHECKING:
    from typing import Union

    if sys.version_info >= (3, 10):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

    Path = Union[str, os.PathLike]


__all__ = [
    'cached_property',
    'read_binary',
    'Collection',
    'Iterable',
    'Iterator',
    'Mapping',
    'Path',
    'ParamSpec',
    'Self',
    'Sequence',
]
