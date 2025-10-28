"""Helpers for cross compatibility across dependency versions."""

from collections.abc import Iterable
from importlib import resources
from typing import Callable

from docutils.nodes import Element


def findall(node: Element) -> Callable[..., Iterable[Element]]:
    """Iterate through"""
    # findall replaces traverse in docutils v0.18
    # note a difference is that findall is an iterator
    return getattr(node, "findall", node.traverse)


# TODO: >= Python 3.9, only use `resources.files` and drop `resources.read_text`
def read_text(module: resources.Package, filename: str) -> str:
    if hasattr(resources, "files"):
        return resources.files(module).joinpath(filename).read_text()
    return resources.read_text(module, filename)
