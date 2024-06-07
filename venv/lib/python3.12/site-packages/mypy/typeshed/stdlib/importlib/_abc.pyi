import sys
import types
from abc import ABCMeta
from importlib.machinery import ModuleSpec

if sys.version_info >= (3, 10):
    class Loader(metaclass=ABCMeta):
        def load_module(self, fullname: str) -> types.ModuleType: ...
        if sys.version_info < (3, 12):
            def module_repr(self, module: types.ModuleType) -> str: ...

        def create_module(self, spec: ModuleSpec) -> types.ModuleType | None: ...
        # Not defined on the actual class for backwards-compatibility reasons,
        # but expected in new code.
        def exec_module(self, module: types.ModuleType) -> None: ...
