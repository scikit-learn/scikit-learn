from __future__ import annotations

from importlib import import_module
from typing import Any

from sphinx.errors import ExtensionError


def import_object(object_name: str, /, source: str = '') -> Any:
    """Import python object by qualname."""
    obj_path = object_name.split('.')
    module_name = obj_path.pop(0)
    try:
        obj = import_module(module_name)
        for name in obj_path:
            module_name += '.' + name
            try:
                obj = getattr(obj, name)
            except AttributeError:
                obj = import_module(module_name)
    except (AttributeError, ImportError) as exc:
        if source:
            msg = f'Could not import {object_name} (needed for {source})'
            raise ExtensionError(msg, exc) from exc
        msg = f'Could not import {object_name}'
        raise ExtensionError(msg, exc) from exc
    return obj
