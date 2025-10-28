"""Serialise objects to a stable representation."""

from __future__ import annotations

import hashlib
import json
import types
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


def stable_hash(obj: Any) -> str:
    """Return a stable hash for a Python data structure.

    We can't just use the md5 of str(obj) as the order of collections
    may be random.
    """
    if isinstance(obj, dict):
        obj = sorted(map(stable_hash, obj.items()))
    if isinstance(obj, list | tuple | set | frozenset):
        obj = sorted(map(stable_hash, obj))
    elif isinstance(obj, type | types.FunctionType):
        # The default repr() of functions includes the ID, which is not ideal.
        # We use the fully qualified name instead.
        obj = f'{obj.__module__}.{obj.__qualname__}'
    return hashlib.md5(str(obj).encode(), usedforsecurity=False).hexdigest()


def stable_str(obj: Any, *, indent: int | None = None) -> str:
    """Return a stable string representation of a Python data structure.

    We can't just use ``str(obj)`` as the order of collections may be random.
    """
    return json.dumps(_stable_str_prep(obj), indent=indent)


def _stable_str_prep(obj: Any) -> dict[str, Any] | list[Any] | str:
    if isinstance(obj, dict):
        # Convert to a sorted dict
        obj = [(_stable_str_prep(k), _stable_str_prep(v)) for k, v in obj.items()]
        obj.sort()
        return dict(obj)
    if isinstance(obj, list | tuple | set | frozenset):
        # Convert to a sorted list
        return sorted(map(_stable_str_prep, obj), key=str)
    if isinstance(obj, type | types.FunctionType):
        # The default repr() of functions includes the ID, which is not ideal.
        # We use the fully qualified name instead.
        return f'{obj.__module__}.{obj.__qualname__}'
    # We can't do any better, just use the string representation
    return str(obj)
