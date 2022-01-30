from importlib.abc import Loader
from types import ModuleType
from typing import Mapping, Sequence

# Signature of `builtins.__import__` should be kept identical to `importlib.__import__`
def __import__(
    name: str,
    globals: Mapping[str, object] | None = ...,
    locals: Mapping[str, object] | None = ...,
    fromlist: Sequence[str] = ...,
    level: int = ...,
) -> ModuleType: ...

# `importlib.import_module` return type should be kept the same as `builtins.__import__`
def import_module(name: str, package: str | None = ...) -> ModuleType: ...
def find_loader(name: str, path: str | None = ...) -> Loader | None: ...
def invalidate_caches() -> None: ...
def reload(module: ModuleType) -> ModuleType: ...
