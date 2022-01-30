"""
    sphinx.ext.autodoc.mock
    ~~~~~~~~~~~~~~~~~~~~~~~

    mock for autodoc

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import contextlib
import os
import sys
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from types import ModuleType
from typing import Any, Generator, Iterator, List, Optional, Sequence, Tuple, Union

from sphinx.util import logging
from sphinx.util.inspect import safe_getattr

logger = logging.getLogger(__name__)


class _MockObject:
    """Used by autodoc_mock_imports."""

    __display_name__ = '_MockObject'
    __name__ = ''
    __sphinx_mock__ = True
    __sphinx_decorator_args__: Tuple[Any, ...] = ()

    def __new__(cls, *args: Any, **kwargs: Any) -> Any:
        if len(args) == 3 and isinstance(args[1], tuple):
            superclass = args[1][-1].__class__
            if superclass is cls:
                # subclassing MockObject
                return _make_subclass(args[0], superclass.__display_name__,
                                      superclass=superclass, attributes=args[2])

        return super().__new__(cls)

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.__qualname__ = self.__name__

    def __len__(self) -> int:
        return 0

    def __contains__(self, key: str) -> bool:
        return False

    def __iter__(self) -> Iterator:
        return iter([])

    def __mro_entries__(self, bases: Tuple) -> Tuple:
        return (self.__class__,)

    def __getitem__(self, key: Any) -> "_MockObject":
        return _make_subclass(str(key), self.__display_name__, self.__class__)()

    def __getattr__(self, key: str) -> "_MockObject":
        return _make_subclass(key, self.__display_name__, self.__class__)()

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        call = self.__class__()
        call.__sphinx_decorator_args__ = args
        return call

    def __repr__(self) -> str:
        return self.__display_name__


def _make_subclass(name: str, module: str, superclass: Any = _MockObject,
                   attributes: Any = None, decorator_args: Tuple = ()) -> Any:
    attrs = {'__module__': module,
             '__display_name__': module + '.' + name,
             '__name__': name,
             '__sphinx_decorator_args__': decorator_args}
    attrs.update(attributes or {})

    return type(name, (superclass,), attrs)


class _MockModule(ModuleType):
    """Used by autodoc_mock_imports."""
    __file__ = os.devnull
    __sphinx_mock__ = True

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.__all__: List[str] = []
        self.__path__: List[str] = []

    def __getattr__(self, name: str) -> _MockObject:
        return _make_subclass(name, self.__name__)()

    def __repr__(self) -> str:
        return self.__name__


class MockLoader(Loader):
    """A loader for mocking."""
    def __init__(self, finder: "MockFinder") -> None:
        super().__init__()
        self.finder = finder

    def create_module(self, spec: ModuleSpec) -> ModuleType:
        logger.debug('[autodoc] adding a mock module as %s!', spec.name)
        self.finder.mocked_modules.append(spec.name)
        return _MockModule(spec.name)

    def exec_module(self, module: ModuleType) -> None:
        pass  # nothing to do


class MockFinder(MetaPathFinder):
    """A finder for mocking."""

    def __init__(self, modnames: List[str]) -> None:
        super().__init__()
        self.modnames = modnames
        self.loader = MockLoader(self)
        self.mocked_modules: List[str] = []

    def find_spec(self, fullname: str, path: Optional[Sequence[Union[bytes, str]]],
                  target: ModuleType = None) -> Optional[ModuleSpec]:
        for modname in self.modnames:
            # check if fullname is (or is a descendant of) one of our targets
            if modname == fullname or fullname.startswith(modname + '.'):
                return ModuleSpec(fullname, self.loader)

        return None

    def invalidate_caches(self) -> None:
        """Invalidate mocked modules on sys.modules."""
        for modname in self.mocked_modules:
            sys.modules.pop(modname, None)


@contextlib.contextmanager
def mock(modnames: List[str]) -> Generator[None, None, None]:
    """Insert mock modules during context::

        with mock(['target.module.name']):
            # mock modules are enabled here
            ...
    """
    try:
        finder = MockFinder(modnames)
        sys.meta_path.insert(0, finder)
        yield
    finally:
        sys.meta_path.remove(finder)
        finder.invalidate_caches()


def ismock(subject: Any) -> bool:
    """Check if the object is mocked."""
    # check the object has '__sphinx_mock__' attribute
    try:
        if safe_getattr(subject, '__sphinx_mock__', None) is None:
            return False
    except AttributeError:
        return False

    # check the object is mocked module
    if isinstance(subject, _MockModule):
        return True

    try:
        # check the object is mocked object
        __mro__ = safe_getattr(type(subject), '__mro__', [])
        if len(__mro__) > 2 and __mro__[-2] is _MockObject:
            # A mocked object has a MRO that ends with (..., _MockObject, object).
            return True
    except AttributeError:
        pass

    return False


def undecorate(subject: _MockObject) -> Any:
    """Unwrap mock if *subject* is decorated by mocked object.

    If not decorated, returns given *subject* itself.
    """
    if ismock(subject) and subject.__sphinx_decorator_args__:
        return subject.__sphinx_decorator_args__[0]
    else:
        return subject
