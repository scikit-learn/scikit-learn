"""Python test discovery, setup and run of test functions."""
import enum
import fnmatch
import inspect
import itertools
import os
import sys
import types
import warnings
from collections import Counter
from collections import defaultdict
from functools import partial
from typing import Any
from typing import Callable
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import Set
from typing import Tuple
from typing import Type
from typing import TYPE_CHECKING
from typing import Union

import py

import _pytest
from _pytest import fixtures
from _pytest import nodes
from _pytest._code import filter_traceback
from _pytest._code import getfslineno
from _pytest._code.code import ExceptionInfo
from _pytest._code.code import TerminalRepr
from _pytest._io import TerminalWriter
from _pytest._io.saferepr import saferepr
from _pytest.compat import ascii_escaped
from _pytest.compat import final
from _pytest.compat import get_default_arg_names
from _pytest.compat import get_real_func
from _pytest.compat import getimfunc
from _pytest.compat import getlocation
from _pytest.compat import is_async_function
from _pytest.compat import is_generator
from _pytest.compat import NOTSET
from _pytest.compat import REGEX_TYPE
from _pytest.compat import safe_getattr
from _pytest.compat import safe_isclass
from _pytest.compat import STRING_TYPES
from _pytest.config import Config
from _pytest.config import ExitCode
from _pytest.config import hookimpl
from _pytest.config.argparsing import Parser
from _pytest.deprecated import FSCOLLECTOR_GETHOOKPROXY_ISINITPATH
from _pytest.fixtures import FuncFixtureInfo
from _pytest.main import Session
from _pytest.mark import MARK_GEN
from _pytest.mark import ParameterSet
from _pytest.mark.structures import get_unpacked_marks
from _pytest.mark.structures import Mark
from _pytest.mark.structures import MarkDecorator
from _pytest.mark.structures import normalize_mark_list
from _pytest.outcomes import fail
from _pytest.outcomes import skip
from _pytest.pathlib import import_path
from _pytest.pathlib import ImportPathMismatchError
from _pytest.pathlib import parts
from _pytest.pathlib import visit
from _pytest.warning_types import PytestCollectionWarning
from _pytest.warning_types import PytestUnhandledCoroutineWarning

if TYPE_CHECKING:
    from typing_extensions import Literal
    from _pytest.fixtures import _Scope


def pytest_addoption(parser: Parser) -> None:
    group = parser.getgroup("general")
    group.addoption(
        "--fixtures",
        "--funcargs",
        action="store_true",
        dest="showfixtures",
        default=False,
        help="show available fixtures, sorted by plugin appearance "
        "(fixtures with leading '_' are only shown with '-v')",
    )
    group.addoption(
        "--fixtures-per-test",
        action="store_true",
        dest="show_fixtures_per_test",
        default=False,
        help="show fixtures per test",
    )
    parser.addini(
        "python_files",
        type="args",
        # NOTE: default is also used in AssertionRewritingHook.
        default=["test_*.py", "*_test.py"],
        help="glob-style file patterns for Python test module discovery",
    )
    parser.addini(
        "python_classes",
        type="args",
        default=["Test"],
        help="prefixes or glob names for Python test class discovery",
    )
    parser.addini(
        "python_functions",
        type="args",
        default=["test"],
        help="prefixes or glob names for Python test function and method discovery",
    )
    parser.addini(
        "disable_test_id_escaping_and_forfeit_all_rights_to_community_support",
        type="bool",
        default=False,
        help="disable string escape non-ascii characters, might cause unwanted "
        "side effects(use at your own risk)",
    )


def pytest_cmdline_main(config: Config) -> Optional[Union[int, ExitCode]]:
    if config.option.showfixtures:
        showfixtures(config)
        return 0
    if config.option.show_fixtures_per_test:
        show_fixtures_per_test(config)
        return 0
    return None


def pytest_generate_tests(metafunc: "Metafunc") -> None:
    for marker in metafunc.definition.iter_markers(name="parametrize"):
        # TODO: Fix this type-ignore (overlapping kwargs).
        metafunc.parametrize(*marker.args, **marker.kwargs, _param_mark=marker)  # type: ignore[misc]


def pytest_configure(config: Config) -> None:
    config.addinivalue_line(
        "markers",
        "parametrize(argnames, argvalues): call a test function multiple "
        "times passing in different arguments in turn. argvalues generally "
        "needs to be a list of values if argnames specifies only one name "
        "or a list of tuples of values if argnames specifies multiple names. "
        "Example: @parametrize('arg1', [1,2]) would lead to two calls of the "
        "decorated test function, one with arg1=1 and another with arg1=2."
        "see https://docs.pytest.org/en/stable/parametrize.html for more info "
        "and examples.",
    )
    config.addinivalue_line(
        "markers",
        "usefixtures(fixturename1, fixturename2, ...): mark tests as needing "
        "all of the specified fixtures. see "
        "https://docs.pytest.org/en/stable/fixture.html#usefixtures ",
    )


def async_warn_and_skip(nodeid: str) -> None:
    msg = "async def functions are not natively supported and have been skipped.\n"
    msg += (
        "You need to install a suitable plugin for your async framework, for example:\n"
    )
    msg += "  - anyio\n"
    msg += "  - pytest-asyncio\n"
    msg += "  - pytest-tornasync\n"
    msg += "  - pytest-trio\n"
    msg += "  - pytest-twisted"
    warnings.warn(PytestUnhandledCoroutineWarning(msg.format(nodeid)))
    skip(msg="async def function and no async plugin installed (see warnings)")


@hookimpl(trylast=True)
def pytest_pyfunc_call(pyfuncitem: "Function") -> Optional[object]:
    testfunction = pyfuncitem.obj
    if is_async_function(testfunction):
        async_warn_and_skip(pyfuncitem.nodeid)
    funcargs = pyfuncitem.funcargs
    testargs = {arg: funcargs[arg] for arg in pyfuncitem._fixtureinfo.argnames}
    result = testfunction(**testargs)
    if hasattr(result, "__await__") or hasattr(result, "__aiter__"):
        async_warn_and_skip(pyfuncitem.nodeid)
    return True


def pytest_collect_file(
    path: py.path.local, parent: nodes.Collector
) -> Optional["Module"]:
    ext = path.ext
    if ext == ".py":
        if not parent.session.isinitpath(path):
            if not path_matches_patterns(
                path, parent.config.getini("python_files") + ["__init__.py"]
            ):
                return None
        ihook = parent.session.gethookproxy(path)
        module: Module = ihook.pytest_pycollect_makemodule(path=path, parent=parent)
        return module
    return None


def path_matches_patterns(path: py.path.local, patterns: Iterable[str]) -> bool:
    """Return whether path matches any of the patterns in the list of globs given."""
    return any(path.fnmatch(pattern) for pattern in patterns)


def pytest_pycollect_makemodule(path: py.path.local, parent) -> "Module":
    if path.basename == "__init__.py":
        pkg: Package = Package.from_parent(parent, fspath=path)
        return pkg
    mod: Module = Module.from_parent(parent, fspath=path)
    return mod


@hookimpl(trylast=True)
def pytest_pycollect_makeitem(collector: "PyCollector", name: str, obj: object):
    # Nothing was collected elsewhere, let's do it here.
    if safe_isclass(obj):
        if collector.istestclass(obj, name):
            return Class.from_parent(collector, name=name, obj=obj)
    elif collector.istestfunction(obj, name):
        # mock seems to store unbound methods (issue473), normalize it.
        obj = getattr(obj, "__func__", obj)
        # We need to try and unwrap the function if it's a functools.partial
        # or a functools.wrapped.
        # We mustn't if it's been wrapped with mock.patch (python 2 only).
        if not (inspect.isfunction(obj) or inspect.isfunction(get_real_func(obj))):
            filename, lineno = getfslineno(obj)
            warnings.warn_explicit(
                message=PytestCollectionWarning(
                    "cannot collect %r because it is not a function." % name
                ),
                category=None,
                filename=str(filename),
                lineno=lineno + 1,
            )
        elif getattr(obj, "__test__", True):
            if is_generator(obj):
                res = Function.from_parent(collector, name=name)
                reason = "yield tests were removed in pytest 4.0 - {name} will be ignored".format(
                    name=name
                )
                res.add_marker(MARK_GEN.xfail(run=False, reason=reason))
                res.warn(PytestCollectionWarning(reason))
            else:
                res = list(collector._genfunctions(name, obj))
            return res


class PyobjMixin:
    _ALLOW_MARKERS = True

    # Function and attributes that the mixin needs (for type-checking only).
    if TYPE_CHECKING:
        name: str = ""
        parent: Optional[nodes.Node] = None
        own_markers: List[Mark] = []

        def getparent(self, cls: Type[nodes._NodeType]) -> Optional[nodes._NodeType]:
            ...

        def listchain(self) -> List[nodes.Node]:
            ...

    @property
    def module(self):
        """Python module object this node was collected from (can be None)."""
        node = self.getparent(Module)
        return node.obj if node is not None else None

    @property
    def cls(self):
        """Python class object this node was collected from (can be None)."""
        node = self.getparent(Class)
        return node.obj if node is not None else None

    @property
    def instance(self):
        """Python instance object this node was collected from (can be None)."""
        node = self.getparent(Instance)
        return node.obj if node is not None else None

    @property
    def obj(self):
        """Underlying Python object."""
        obj = getattr(self, "_obj", None)
        if obj is None:
            self._obj = obj = self._getobj()
            # XXX evil hack
            # used to avoid Instance collector marker duplication
            if self._ALLOW_MARKERS:
                self.own_markers.extend(get_unpacked_marks(self.obj))
        return obj

    @obj.setter
    def obj(self, value):
        self._obj = value

    def _getobj(self):
        """Get the underlying Python object. May be overwritten by subclasses."""
        # TODO: Improve the type of `parent` such that assert/ignore aren't needed.
        assert self.parent is not None
        obj = self.parent.obj  # type: ignore[attr-defined]
        return getattr(obj, self.name)

    def getmodpath(self, stopatmodule: bool = True, includemodule: bool = False) -> str:
        """Return Python path relative to the containing module."""
        chain = self.listchain()
        chain.reverse()
        parts = []
        for node in chain:
            if isinstance(node, Instance):
                continue
            name = node.name
            if isinstance(node, Module):
                name = os.path.splitext(name)[0]
                if stopatmodule:
                    if includemodule:
                        parts.append(name)
                    break
            parts.append(name)
        parts.reverse()
        return ".".join(parts)

    def reportinfo(self) -> Tuple[Union[py.path.local, str], int, str]:
        # XXX caching?
        obj = self.obj
        compat_co_firstlineno = getattr(obj, "compat_co_firstlineno", None)
        if isinstance(compat_co_firstlineno, int):
            # nose compatibility
            file_path = sys.modules[obj.__module__].__file__
            if file_path.endswith(".pyc"):
                file_path = file_path[:-1]
            fspath: Union[py.path.local, str] = file_path
            lineno = compat_co_firstlineno
        else:
            fspath, lineno = getfslineno(obj)
        modpath = self.getmodpath()
        assert isinstance(lineno, int)
        return fspath, lineno, modpath


# As an optimization, these builtin attribute names are pre-ignored when
# iterating over an object during collection -- the pytest_pycollect_makeitem
# hook is not called for them.
# fmt: off
class _EmptyClass: pass  # noqa: E701
IGNORED_ATTRIBUTES = frozenset.union(  # noqa: E305
    frozenset(),
    # Module.
    dir(types.ModuleType("empty_module")),
    # Some extra module attributes the above doesn't catch.
    {"__builtins__", "__file__", "__cached__"},
    # Class.
    dir(_EmptyClass),
    # Instance.
    dir(_EmptyClass()),
)
del _EmptyClass
# fmt: on


class PyCollector(PyobjMixin, nodes.Collector):
    def funcnamefilter(self, name: str) -> bool:
        return self._matches_prefix_or_glob_option("python_functions", name)

    def isnosetest(self, obj: object) -> bool:
        """Look for the __test__ attribute, which is applied by the
        @nose.tools.istest decorator.
        """
        # We explicitly check for "is True" here to not mistakenly treat
        # classes with a custom __getattr__ returning something truthy (like a
        # function) as test classes.
        return safe_getattr(obj, "__test__", False) is True

    def classnamefilter(self, name: str) -> bool:
        return self._matches_prefix_or_glob_option("python_classes", name)

    def istestfunction(self, obj: object, name: str) -> bool:
        if self.funcnamefilter(name) or self.isnosetest(obj):
            if isinstance(obj, staticmethod):
                # staticmethods need to be unwrapped.
                obj = safe_getattr(obj, "__func__", False)
            return (
                safe_getattr(obj, "__call__", False)
                and fixtures.getfixturemarker(obj) is None
            )
        else:
            return False

    def istestclass(self, obj: object, name: str) -> bool:
        return self.classnamefilter(name) or self.isnosetest(obj)

    def _matches_prefix_or_glob_option(self, option_name: str, name: str) -> bool:
        """Check if the given name matches the prefix or glob-pattern defined
        in ini configuration."""
        for option in self.config.getini(option_name):
            if name.startswith(option):
                return True
            # Check that name looks like a glob-string before calling fnmatch
            # because this is called for every name in each collected module,
            # and fnmatch is somewhat expensive to call.
            elif ("*" in option or "?" in option or "[" in option) and fnmatch.fnmatch(
                name, option
            ):
                return True
        return False

    def collect(self) -> Iterable[Union[nodes.Item, nodes.Collector]]:
        if not getattr(self.obj, "__test__", True):
            return []

        # NB. we avoid random getattrs and peek in the __dict__ instead
        # (XXX originally introduced from a PyPy need, still true?)
        dicts = [getattr(self.obj, "__dict__", {})]
        for basecls in self.obj.__class__.__mro__:
            dicts.append(basecls.__dict__)
        seen: Set[str] = set()
        values: List[Union[nodes.Item, nodes.Collector]] = []
        ihook = self.ihook
        for dic in dicts:
            # Note: seems like the dict can change during iteration -
            # be careful not to remove the list() without consideration.
            for name, obj in list(dic.items()):
                if name in IGNORED_ATTRIBUTES:
                    continue
                if name in seen:
                    continue
                seen.add(name)
                res = ihook.pytest_pycollect_makeitem(
                    collector=self, name=name, obj=obj
                )
                if res is None:
                    continue
                elif isinstance(res, list):
                    values.extend(res)
                else:
                    values.append(res)

        def sort_key(item):
            fspath, lineno, _ = item.reportinfo()
            return (str(fspath), lineno)

        values.sort(key=sort_key)
        return values

    def _genfunctions(self, name: str, funcobj) -> Iterator["Function"]:
        modulecol = self.getparent(Module)
        assert modulecol is not None
        module = modulecol.obj
        clscol = self.getparent(Class)
        cls = clscol and clscol.obj or None
        fm = self.session._fixturemanager

        definition = FunctionDefinition.from_parent(self, name=name, callobj=funcobj)
        fixtureinfo = definition._fixtureinfo

        metafunc = Metafunc(
            definition, fixtureinfo, self.config, cls=cls, module=module
        )
        methods = []
        if hasattr(module, "pytest_generate_tests"):
            methods.append(module.pytest_generate_tests)
        if cls is not None and hasattr(cls, "pytest_generate_tests"):
            methods.append(cls().pytest_generate_tests)

        self.ihook.pytest_generate_tests.call_extra(methods, dict(metafunc=metafunc))

        if not metafunc._calls:
            yield Function.from_parent(self, name=name, fixtureinfo=fixtureinfo)
        else:
            # Add funcargs() as fixturedefs to fixtureinfo.arg2fixturedefs.
            fixtures.add_funcarg_pseudo_fixture_def(self, metafunc, fm)

            # Add_funcarg_pseudo_fixture_def may have shadowed some fixtures
            # with direct parametrization, so make sure we update what the
            # function really needs.
            fixtureinfo.prune_dependency_tree()

            for callspec in metafunc._calls:
                subname = f"{name}[{callspec.id}]"
                yield Function.from_parent(
                    self,
                    name=subname,
                    callspec=callspec,
                    callobj=funcobj,
                    fixtureinfo=fixtureinfo,
                    keywords={callspec.id: True},
                    originalname=name,
                )


class Module(nodes.File, PyCollector):
    """Collector for test classes and functions."""

    def _getobj(self):
        return self._importtestmodule()

    def collect(self) -> Iterable[Union[nodes.Item, nodes.Collector]]:
        self._inject_setup_module_fixture()
        self._inject_setup_function_fixture()
        self.session._fixturemanager.parsefactories(self)
        return super().collect()

    def _inject_setup_module_fixture(self) -> None:
        """Inject a hidden autouse, module scoped fixture into the collected module object
        that invokes setUpModule/tearDownModule if either or both are available.

        Using a fixture to invoke this methods ensures we play nicely and unsurprisingly with
        other fixtures (#517).
        """
        setup_module = _get_first_non_fixture_func(
            self.obj, ("setUpModule", "setup_module")
        )
        teardown_module = _get_first_non_fixture_func(
            self.obj, ("tearDownModule", "teardown_module")
        )

        if setup_module is None and teardown_module is None:
            return

        @fixtures.fixture(
            autouse=True,
            scope="module",
            # Use a unique name to speed up lookup.
            name=f"xunit_setup_module_fixture_{self.obj.__name__}",
        )
        def xunit_setup_module_fixture(request) -> Generator[None, None, None]:
            if setup_module is not None:
                _call_with_optional_argument(setup_module, request.module)
            yield
            if teardown_module is not None:
                _call_with_optional_argument(teardown_module, request.module)

        self.obj.__pytest_setup_module = xunit_setup_module_fixture

    def _inject_setup_function_fixture(self) -> None:
        """Inject a hidden autouse, function scoped fixture into the collected module object
        that invokes setup_function/teardown_function if either or both are available.

        Using a fixture to invoke this methods ensures we play nicely and unsurprisingly with
        other fixtures (#517).
        """
        setup_function = _get_first_non_fixture_func(self.obj, ("setup_function",))
        teardown_function = _get_first_non_fixture_func(
            self.obj, ("teardown_function",)
        )
        if setup_function is None and teardown_function is None:
            return

        @fixtures.fixture(
            autouse=True,
            scope="function",
            # Use a unique name to speed up lookup.
            name=f"xunit_setup_function_fixture_{self.obj.__name__}",
        )
        def xunit_setup_function_fixture(request) -> Generator[None, None, None]:
            if request.instance is not None:
                # in this case we are bound to an instance, so we need to let
                # setup_method handle this
                yield
                return
            if setup_function is not None:
                _call_with_optional_argument(setup_function, request.function)
            yield
            if teardown_function is not None:
                _call_with_optional_argument(teardown_function, request.function)

        self.obj.__pytest_setup_function = xunit_setup_function_fixture

    def _importtestmodule(self):
        # We assume we are only called once per module.
        importmode = self.config.getoption("--import-mode")
        try:
            mod = import_path(self.fspath, mode=importmode)
        except SyntaxError as e:
            raise self.CollectError(
                ExceptionInfo.from_current().getrepr(style="short")
            ) from e
        except ImportPathMismatchError as e:
            raise self.CollectError(
                "import file mismatch:\n"
                "imported module %r has this __file__ attribute:\n"
                "  %s\n"
                "which is not the same as the test file we want to collect:\n"
                "  %s\n"
                "HINT: remove __pycache__ / .pyc files and/or use a "
                "unique basename for your test file modules" % e.args
            ) from e
        except ImportError as e:
            exc_info = ExceptionInfo.from_current()
            if self.config.getoption("verbose") < 2:
                exc_info.traceback = exc_info.traceback.filter(filter_traceback)
            exc_repr = (
                exc_info.getrepr(style="short")
                if exc_info.traceback
                else exc_info.exconly()
            )
            formatted_tb = str(exc_repr)
            raise self.CollectError(
                "ImportError while importing test module '{fspath}'.\n"
                "Hint: make sure your test modules/packages have valid Python names.\n"
                "Traceback:\n"
                "{traceback}".format(fspath=self.fspath, traceback=formatted_tb)
            ) from e
        except skip.Exception as e:
            if e.allow_module_level:
                raise
            raise self.CollectError(
                "Using pytest.skip outside of a test is not allowed. "
                "To decorate a test function, use the @pytest.mark.skip "
                "or @pytest.mark.skipif decorators instead, and to skip a "
                "module use `pytestmark = pytest.mark.{skip,skipif}."
            ) from e
        self.config.pluginmanager.consider_module(mod)
        return mod


class Package(Module):
    def __init__(
        self,
        fspath: py.path.local,
        parent: nodes.Collector,
        # NOTE: following args are unused:
        config=None,
        session=None,
        nodeid=None,
    ) -> None:
        # NOTE: Could be just the following, but kept as-is for compat.
        # nodes.FSCollector.__init__(self, fspath, parent=parent)
        session = parent.session
        nodes.FSCollector.__init__(
            self, fspath, parent=parent, config=config, session=session, nodeid=nodeid
        )
        self.name = os.path.basename(str(fspath.dirname))

    def setup(self) -> None:
        # Not using fixtures to call setup_module here because autouse fixtures
        # from packages are not called automatically (#4085).
        setup_module = _get_first_non_fixture_func(
            self.obj, ("setUpModule", "setup_module")
        )
        if setup_module is not None:
            _call_with_optional_argument(setup_module, self.obj)

        teardown_module = _get_first_non_fixture_func(
            self.obj, ("tearDownModule", "teardown_module")
        )
        if teardown_module is not None:
            func = partial(_call_with_optional_argument, teardown_module, self.obj)
            self.addfinalizer(func)

    def gethookproxy(self, fspath: py.path.local):
        warnings.warn(FSCOLLECTOR_GETHOOKPROXY_ISINITPATH, stacklevel=2)
        return self.session.gethookproxy(fspath)

    def isinitpath(self, path: py.path.local) -> bool:
        warnings.warn(FSCOLLECTOR_GETHOOKPROXY_ISINITPATH, stacklevel=2)
        return self.session.isinitpath(path)

    def _recurse(self, direntry: "os.DirEntry[str]") -> bool:
        if direntry.name == "__pycache__":
            return False
        path = py.path.local(direntry.path)
        ihook = self.session.gethookproxy(path.dirpath())
        if ihook.pytest_ignore_collect(path=path, config=self.config):
            return False
        norecursepatterns = self.config.getini("norecursedirs")
        if any(path.check(fnmatch=pat) for pat in norecursepatterns):
            return False
        return True

    def _collectfile(
        self, path: py.path.local, handle_dupes: bool = True
    ) -> Sequence[nodes.Collector]:
        assert (
            path.isfile()
        ), "{!r} is not a file (isdir={!r}, exists={!r}, islink={!r})".format(
            path, path.isdir(), path.exists(), path.islink()
        )
        ihook = self.session.gethookproxy(path)
        if not self.session.isinitpath(path):
            if ihook.pytest_ignore_collect(path=path, config=self.config):
                return ()

        if handle_dupes:
            keepduplicates = self.config.getoption("keepduplicates")
            if not keepduplicates:
                duplicate_paths = self.config.pluginmanager._duplicatepaths
                if path in duplicate_paths:
                    return ()
                else:
                    duplicate_paths.add(path)

        return ihook.pytest_collect_file(path=path, parent=self)  # type: ignore[no-any-return]

    def collect(self) -> Iterable[Union[nodes.Item, nodes.Collector]]:
        this_path = self.fspath.dirpath()
        init_module = this_path.join("__init__.py")
        if init_module.check(file=1) and path_matches_patterns(
            init_module, self.config.getini("python_files")
        ):
            yield Module.from_parent(self, fspath=init_module)
        pkg_prefixes: Set[py.path.local] = set()
        for direntry in visit(str(this_path), recurse=self._recurse):
            path = py.path.local(direntry.path)

            # We will visit our own __init__.py file, in which case we skip it.
            if direntry.is_file():
                if direntry.name == "__init__.py" and path.dirpath() == this_path:
                    continue

            parts_ = parts(direntry.path)
            if any(
                str(pkg_prefix) in parts_ and pkg_prefix.join("__init__.py") != path
                for pkg_prefix in pkg_prefixes
            ):
                continue

            if direntry.is_file():
                yield from self._collectfile(path)
            elif not direntry.is_dir():
                # Broken symlink or invalid/missing file.
                continue
            elif path.join("__init__.py").check(file=1):
                pkg_prefixes.add(path)


def _call_with_optional_argument(func, arg) -> None:
    """Call the given function with the given argument if func accepts one argument, otherwise
    calls func without arguments."""
    arg_count = func.__code__.co_argcount
    if inspect.ismethod(func):
        arg_count -= 1
    if arg_count:
        func(arg)
    else:
        func()


def _get_first_non_fixture_func(obj: object, names: Iterable[str]):
    """Return the attribute from the given object to be used as a setup/teardown
    xunit-style function, but only if not marked as a fixture to avoid calling it twice."""
    for name in names:
        meth = getattr(obj, name, None)
        if meth is not None and fixtures.getfixturemarker(meth) is None:
            return meth


class Class(PyCollector):
    """Collector for test methods."""

    @classmethod
    def from_parent(cls, parent, *, name, obj=None):
        """The public constructor."""
        return super().from_parent(name=name, parent=parent)

    def collect(self) -> Iterable[Union[nodes.Item, nodes.Collector]]:
        if not safe_getattr(self.obj, "__test__", True):
            return []
        if hasinit(self.obj):
            assert self.parent is not None
            self.warn(
                PytestCollectionWarning(
                    "cannot collect test class %r because it has a "
                    "__init__ constructor (from: %s)"
                    % (self.obj.__name__, self.parent.nodeid)
                )
            )
            return []
        elif hasnew(self.obj):
            assert self.parent is not None
            self.warn(
                PytestCollectionWarning(
                    "cannot collect test class %r because it has a "
                    "__new__ constructor (from: %s)"
                    % (self.obj.__name__, self.parent.nodeid)
                )
            )
            return []

        self._inject_setup_class_fixture()
        self._inject_setup_method_fixture()

        return [Instance.from_parent(self, name="()")]

    def _inject_setup_class_fixture(self) -> None:
        """Inject a hidden autouse, class scoped fixture into the collected class object
        that invokes setup_class/teardown_class if either or both are available.

        Using a fixture to invoke this methods ensures we play nicely and unsurprisingly with
        other fixtures (#517).
        """
        setup_class = _get_first_non_fixture_func(self.obj, ("setup_class",))
        teardown_class = getattr(self.obj, "teardown_class", None)
        if setup_class is None and teardown_class is None:
            return

        @fixtures.fixture(
            autouse=True,
            scope="class",
            # Use a unique name to speed up lookup.
            name=f"xunit_setup_class_fixture_{self.obj.__qualname__}",
        )
        def xunit_setup_class_fixture(cls) -> Generator[None, None, None]:
            if setup_class is not None:
                func = getimfunc(setup_class)
                _call_with_optional_argument(func, self.obj)
            yield
            if teardown_class is not None:
                func = getimfunc(teardown_class)
                _call_with_optional_argument(func, self.obj)

        self.obj.__pytest_setup_class = xunit_setup_class_fixture

    def _inject_setup_method_fixture(self) -> None:
        """Inject a hidden autouse, function scoped fixture into the collected class object
        that invokes setup_method/teardown_method if either or both are available.

        Using a fixture to invoke this methods ensures we play nicely and unsurprisingly with
        other fixtures (#517).
        """
        setup_method = _get_first_non_fixture_func(self.obj, ("setup_method",))
        teardown_method = getattr(self.obj, "teardown_method", None)
        if setup_method is None and teardown_method is None:
            return

        @fixtures.fixture(
            autouse=True,
            scope="function",
            # Use a unique name to speed up lookup.
            name=f"xunit_setup_method_fixture_{self.obj.__qualname__}",
        )
        def xunit_setup_method_fixture(self, request) -> Generator[None, None, None]:
            method = request.function
            if setup_method is not None:
                func = getattr(self, "setup_method")
                _call_with_optional_argument(func, method)
            yield
            if teardown_method is not None:
                func = getattr(self, "teardown_method")
                _call_with_optional_argument(func, method)

        self.obj.__pytest_setup_method = xunit_setup_method_fixture


class Instance(PyCollector):
    _ALLOW_MARKERS = False  # hack, destroy later
    # Instances share the object with their parents in a way
    # that duplicates markers instances if not taken out
    # can be removed at node structure reorganization time.

    def _getobj(self):
        # TODO: Improve the type of `parent` such that assert/ignore aren't needed.
        assert self.parent is not None
        obj = self.parent.obj  # type: ignore[attr-defined]
        return obj()

    def collect(self) -> Iterable[Union[nodes.Item, nodes.Collector]]:
        self.session._fixturemanager.parsefactories(self)
        return super().collect()

    def newinstance(self):
        self.obj = self._getobj()
        return self.obj


def hasinit(obj: object) -> bool:
    init: object = getattr(obj, "__init__", None)
    if init:
        return init != object.__init__
    return False


def hasnew(obj: object) -> bool:
    new: object = getattr(obj, "__new__", None)
    if new:
        return new != object.__new__
    return False


@final
class CallSpec2:
    def __init__(self, metafunc: "Metafunc") -> None:
        self.metafunc = metafunc
        self.funcargs: Dict[str, object] = {}
        self._idlist: List[str] = []
        self.params: Dict[str, object] = {}
        # Used for sorting parametrized resources.
        self._arg2scopenum: Dict[str, int] = {}
        self.marks: List[Mark] = []
        self.indices: Dict[str, int] = {}

    def copy(self) -> "CallSpec2":
        cs = CallSpec2(self.metafunc)
        cs.funcargs.update(self.funcargs)
        cs.params.update(self.params)
        cs.marks.extend(self.marks)
        cs.indices.update(self.indices)
        cs._arg2scopenum.update(self._arg2scopenum)
        cs._idlist = list(self._idlist)
        return cs

    def _checkargnotcontained(self, arg: str) -> None:
        if arg in self.params or arg in self.funcargs:
            raise ValueError(f"duplicate {arg!r}")

    def getparam(self, name: str) -> object:
        try:
            return self.params[name]
        except KeyError as e:
            raise ValueError(name) from e

    @property
    def id(self) -> str:
        return "-".join(map(str, self._idlist))

    def setmulti2(
        self,
        valtypes: Mapping[str, "Literal['params', 'funcargs']"],
        argnames: Sequence[str],
        valset: Iterable[object],
        id: str,
        marks: Iterable[Union[Mark, MarkDecorator]],
        scopenum: int,
        param_index: int,
    ) -> None:
        for arg, val in zip(argnames, valset):
            self._checkargnotcontained(arg)
            valtype_for_arg = valtypes[arg]
            if valtype_for_arg == "params":
                self.params[arg] = val
            elif valtype_for_arg == "funcargs":
                self.funcargs[arg] = val
            else:  # pragma: no cover
                assert False, f"Unhandled valtype for arg: {valtype_for_arg}"
            self.indices[arg] = param_index
            self._arg2scopenum[arg] = scopenum
        self._idlist.append(id)
        self.marks.extend(normalize_mark_list(marks))


@final
class Metafunc:
    """Objects passed to the :func:`pytest_generate_tests <_pytest.hookspec.pytest_generate_tests>` hook.

    They help to inspect a test function and to generate tests according to
    test configuration or values specified in the class or module where a
    test function is defined.
    """

    def __init__(
        self,
        definition: "FunctionDefinition",
        fixtureinfo: fixtures.FuncFixtureInfo,
        config: Config,
        cls=None,
        module=None,
    ) -> None:
        #: Access to the underlying :class:`_pytest.python.FunctionDefinition`.
        self.definition = definition

        #: Access to the :class:`_pytest.config.Config` object for the test session.
        self.config = config

        #: The module object where the test function is defined in.
        self.module = module

        #: Underlying Python test function.
        self.function = definition.obj

        #: Set of fixture names required by the test function.
        self.fixturenames = fixtureinfo.names_closure

        #: Class object where the test function is defined in or ``None``.
        self.cls = cls

        self._calls: List[CallSpec2] = []
        self._arg2fixturedefs = fixtureinfo.name2fixturedefs

    def parametrize(
        self,
        argnames: Union[str, List[str], Tuple[str, ...]],
        argvalues: Iterable[Union[ParameterSet, Sequence[object], object]],
        indirect: Union[bool, Sequence[str]] = False,
        ids: Optional[
            Union[
                Iterable[Union[None, str, float, int, bool]],
                Callable[[Any], Optional[object]],
            ]
        ] = None,
        scope: "Optional[_Scope]" = None,
        *,
        _param_mark: Optional[Mark] = None,
    ) -> None:
        """Add new invocations to the underlying test function using the list
        of argvalues for the given argnames.  Parametrization is performed
        during the collection phase.  If you need to setup expensive resources
        see about setting indirect to do it rather at test setup time.

        :param argnames:
            A comma-separated string denoting one or more argument names, or
            a list/tuple of argument strings.

        :param argvalues:
            The list of argvalues determines how often a test is invoked with
            different argument values.

            If only one argname was specified argvalues is a list of values.
            If N argnames were specified, argvalues must be a list of
            N-tuples, where each tuple-element specifies a value for its
            respective argname.

        :param indirect:
            A list of arguments' names (subset of argnames) or a boolean.
            If True the list contains all names from the argnames. Each
            argvalue corresponding to an argname in this list will
            be passed as request.param to its respective argname fixture
            function so that it can perform more expensive setups during the
            setup phase of a test rather than at collection time.

        :param ids:
            Sequence of (or generator for) ids for ``argvalues``,
            or a callable to return part of the id for each argvalue.

            With sequences (and generators like ``itertools.count()``) the
            returned ids should be of type ``string``, ``int``, ``float``,
            ``bool``, or ``None``.
            They are mapped to the corresponding index in ``argvalues``.
            ``None`` means to use the auto-generated id.

            If it is a callable it will be called for each entry in
            ``argvalues``, and the return value is used as part of the
            auto-generated id for the whole set (where parts are joined with
            dashes ("-")).
            This is useful to provide more specific ids for certain items, e.g.
            dates.  Returning ``None`` will use an auto-generated id.

            If no ids are provided they will be generated automatically from
            the argvalues.

        :param scope:
            If specified it denotes the scope of the parameters.
            The scope is used for grouping tests by parameter instances.
            It will also override any fixture-function defined scope, allowing
            to set a dynamic scope using test context or configuration.
        """
        from _pytest.fixtures import scope2index

        argnames, parameters = ParameterSet._for_parametrize(
            argnames,
            argvalues,
            self.function,
            self.config,
            nodeid=self.definition.nodeid,
        )
        del argvalues

        if "request" in argnames:
            fail(
                "'request' is a reserved name and cannot be used in @pytest.mark.parametrize",
                pytrace=False,
            )

        if scope is None:
            scope = _find_parametrized_scope(argnames, self._arg2fixturedefs, indirect)

        self._validate_if_using_arg_names(argnames, indirect)

        arg_values_types = self._resolve_arg_value_types(argnames, indirect)

        # Use any already (possibly) generated ids with parametrize Marks.
        if _param_mark and _param_mark._param_ids_from:
            generated_ids = _param_mark._param_ids_from._param_ids_generated
            if generated_ids is not None:
                ids = generated_ids

        ids = self._resolve_arg_ids(
            argnames, ids, parameters, nodeid=self.definition.nodeid
        )

        # Store used (possibly generated) ids with parametrize Marks.
        if _param_mark and _param_mark._param_ids_from and generated_ids is None:
            object.__setattr__(_param_mark._param_ids_from, "_param_ids_generated", ids)

        scopenum = scope2index(
            scope, descr=f"parametrize() call in {self.function.__name__}"
        )

        # Create the new calls: if we are parametrize() multiple times (by applying the decorator
        # more than once) then we accumulate those calls generating the cartesian product
        # of all calls.
        newcalls = []
        for callspec in self._calls or [CallSpec2(self)]:
            for param_index, (param_id, param_set) in enumerate(zip(ids, parameters)):
                newcallspec = callspec.copy()
                newcallspec.setmulti2(
                    arg_values_types,
                    argnames,
                    param_set.values,
                    param_id,
                    param_set.marks,
                    scopenum,
                    param_index,
                )
                newcalls.append(newcallspec)
        self._calls = newcalls

    def _resolve_arg_ids(
        self,
        argnames: Sequence[str],
        ids: Optional[
            Union[
                Iterable[Union[None, str, float, int, bool]],
                Callable[[Any], Optional[object]],
            ]
        ],
        parameters: Sequence[ParameterSet],
        nodeid: str,
    ) -> List[str]:
        """Resolve the actual ids for the given argnames, based on the ``ids`` parameter given
        to ``parametrize``.

        :param List[str] argnames: List of argument names passed to ``parametrize()``.
        :param ids: The ids parameter of the parametrized call (see docs).
        :param List[ParameterSet] parameters: The list of parameter values, same size as ``argnames``.
        :param str str: The nodeid of the item that generated this parametrized call.
        :rtype: List[str]
        :returns: The list of ids for each argname given.
        """
        if ids is None:
            idfn = None
            ids_ = None
        elif callable(ids):
            idfn = ids
            ids_ = None
        else:
            idfn = None
            ids_ = self._validate_ids(ids, parameters, self.function.__name__)
        return idmaker(argnames, parameters, idfn, ids_, self.config, nodeid=nodeid)

    def _validate_ids(
        self,
        ids: Iterable[Union[None, str, float, int, bool]],
        parameters: Sequence[ParameterSet],
        func_name: str,
    ) -> List[Union[None, str]]:
        try:
            num_ids = len(ids)  # type: ignore[arg-type]
        except TypeError:
            try:
                iter(ids)
            except TypeError as e:
                raise TypeError("ids must be a callable or an iterable") from e
            num_ids = len(parameters)

        # num_ids == 0 is a special case: https://github.com/pytest-dev/pytest/issues/1849
        if num_ids != len(parameters) and num_ids != 0:
            msg = "In {}: {} parameter sets specified, with different number of ids: {}"
            fail(msg.format(func_name, len(parameters), num_ids), pytrace=False)

        new_ids = []
        for idx, id_value in enumerate(itertools.islice(ids, num_ids)):
            if id_value is None or isinstance(id_value, str):
                new_ids.append(id_value)
            elif isinstance(id_value, (float, int, bool)):
                new_ids.append(str(id_value))
            else:
                msg = (  # type: ignore[unreachable]
                    "In {}: ids must be list of string/float/int/bool, "
                    "found: {} (type: {!r}) at index {}"
                )
                fail(
                    msg.format(func_name, saferepr(id_value), type(id_value), idx),
                    pytrace=False,
                )
        return new_ids

    def _resolve_arg_value_types(
        self, argnames: Sequence[str], indirect: Union[bool, Sequence[str]],
    ) -> Dict[str, "Literal['params', 'funcargs']"]:
        """Resolve if each parametrized argument must be considered a
        parameter to a fixture or a "funcarg" to the function, based on the
        ``indirect`` parameter of the parametrized() call.

        :param List[str] argnames: List of argument names passed to ``parametrize()``.
        :param indirect: Same as the ``indirect`` parameter of ``parametrize()``.
        :rtype: Dict[str, str]
            A dict mapping each arg name to either:
            * "params" if the argname should be the parameter of a fixture of the same name.
            * "funcargs" if the argname should be a parameter to the parametrized test function.
        """
        if isinstance(indirect, bool):
            valtypes: Dict[str, Literal["params", "funcargs"]] = dict.fromkeys(
                argnames, "params" if indirect else "funcargs"
            )
        elif isinstance(indirect, Sequence):
            valtypes = dict.fromkeys(argnames, "funcargs")
            for arg in indirect:
                if arg not in argnames:
                    fail(
                        "In {}: indirect fixture '{}' doesn't exist".format(
                            self.function.__name__, arg
                        ),
                        pytrace=False,
                    )
                valtypes[arg] = "params"
        else:
            fail(
                "In {func}: expected Sequence or boolean for indirect, got {type}".format(
                    type=type(indirect).__name__, func=self.function.__name__
                ),
                pytrace=False,
            )
        return valtypes

    def _validate_if_using_arg_names(
        self, argnames: Sequence[str], indirect: Union[bool, Sequence[str]],
    ) -> None:
        """Check if all argnames are being used, by default values, or directly/indirectly.

        :param List[str] argnames: List of argument names passed to ``parametrize()``.
        :param indirect: Same as the ``indirect`` parameter of ``parametrize()``.
        :raises ValueError: If validation fails.
        """
        default_arg_names = set(get_default_arg_names(self.function))
        func_name = self.function.__name__
        for arg in argnames:
            if arg not in self.fixturenames:
                if arg in default_arg_names:
                    fail(
                        "In {}: function already takes an argument '{}' with a default value".format(
                            func_name, arg
                        ),
                        pytrace=False,
                    )
                else:
                    if isinstance(indirect, Sequence):
                        name = "fixture" if arg in indirect else "argument"
                    else:
                        name = "fixture" if indirect else "argument"
                    fail(
                        f"In {func_name}: function uses no {name} '{arg}'",
                        pytrace=False,
                    )


def _find_parametrized_scope(
    argnames: Sequence[str],
    arg2fixturedefs: Mapping[str, Sequence[fixtures.FixtureDef[object]]],
    indirect: Union[bool, Sequence[str]],
) -> "fixtures._Scope":
    """Find the most appropriate scope for a parametrized call based on its arguments.

    When there's at least one direct argument, always use "function" scope.

    When a test function is parametrized and all its arguments are indirect
    (e.g. fixtures), return the most narrow scope based on the fixtures used.

    Related to issue #1832, based on code posted by @Kingdread.
    """
    if isinstance(indirect, Sequence):
        all_arguments_are_fixtures = len(indirect) == len(argnames)
    else:
        all_arguments_are_fixtures = bool(indirect)

    if all_arguments_are_fixtures:
        fixturedefs = arg2fixturedefs or {}
        used_scopes = [
            fixturedef[0].scope
            for name, fixturedef in fixturedefs.items()
            if name in argnames
        ]
        if used_scopes:
            # Takes the most narrow scope from used fixtures.
            for scope in reversed(fixtures.scopes):
                if scope in used_scopes:
                    return scope

    return "function"


def _ascii_escaped_by_config(val: Union[str, bytes], config: Optional[Config]) -> str:
    if config is None:
        escape_option = False
    else:
        escape_option = config.getini(
            "disable_test_id_escaping_and_forfeit_all_rights_to_community_support"
        )
    # TODO: If escaping is turned off and the user passes bytes,
    #       will return a bytes. For now we ignore this but the
    #       code *probably* doesn't handle this case.
    return val if escape_option else ascii_escaped(val)  # type: ignore


def _idval(
    val: object,
    argname: str,
    idx: int,
    idfn: Optional[Callable[[Any], Optional[object]]],
    nodeid: Optional[str],
    config: Optional[Config],
) -> str:
    if idfn:
        try:
            generated_id = idfn(val)
            if generated_id is not None:
                val = generated_id
        except Exception as e:
            prefix = f"{nodeid}: " if nodeid is not None else ""
            msg = "error raised while trying to determine id of parameter '{}' at position {}"
            msg = prefix + msg.format(argname, idx)
            raise ValueError(msg) from e
    elif config:
        hook_id: Optional[str] = config.hook.pytest_make_parametrize_id(
            config=config, val=val, argname=argname
        )
        if hook_id:
            return hook_id

    if isinstance(val, STRING_TYPES):
        return _ascii_escaped_by_config(val, config)
    elif val is None or isinstance(val, (float, int, bool)):
        return str(val)
    elif isinstance(val, REGEX_TYPE):
        return ascii_escaped(val.pattern)
    elif val is NOTSET:
        # Fallback to default. Note that NOTSET is an enum.Enum.
        pass
    elif isinstance(val, enum.Enum):
        return str(val)
    elif isinstance(getattr(val, "__name__", None), str):
        # Name of a class, function, module, etc.
        name: str = getattr(val, "__name__")
        return name
    return str(argname) + str(idx)


def _idvalset(
    idx: int,
    parameterset: ParameterSet,
    argnames: Iterable[str],
    idfn: Optional[Callable[[Any], Optional[object]]],
    ids: Optional[List[Union[None, str]]],
    nodeid: Optional[str],
    config: Optional[Config],
) -> str:
    if parameterset.id is not None:
        return parameterset.id
    id = None if ids is None or idx >= len(ids) else ids[idx]
    if id is None:
        this_id = [
            _idval(val, argname, idx, idfn, nodeid=nodeid, config=config)
            for val, argname in zip(parameterset.values, argnames)
        ]
        return "-".join(this_id)
    else:
        return _ascii_escaped_by_config(id, config)


def idmaker(
    argnames: Iterable[str],
    parametersets: Iterable[ParameterSet],
    idfn: Optional[Callable[[Any], Optional[object]]] = None,
    ids: Optional[List[Union[None, str]]] = None,
    config: Optional[Config] = None,
    nodeid: Optional[str] = None,
) -> List[str]:
    resolved_ids = [
        _idvalset(
            valindex, parameterset, argnames, idfn, ids, config=config, nodeid=nodeid
        )
        for valindex, parameterset in enumerate(parametersets)
    ]

    # All IDs must be unique!
    unique_ids = set(resolved_ids)
    if len(unique_ids) != len(resolved_ids):

        # Record the number of occurrences of each test ID.
        test_id_counts = Counter(resolved_ids)

        # Map the test ID to its next suffix.
        test_id_suffixes: Dict[str, int] = defaultdict(int)

        # Suffix non-unique IDs to make them unique.
        for index, test_id in enumerate(resolved_ids):
            if test_id_counts[test_id] > 1:
                resolved_ids[index] = "{}{}".format(test_id, test_id_suffixes[test_id])
                test_id_suffixes[test_id] += 1

    return resolved_ids


def show_fixtures_per_test(config):
    from _pytest.main import wrap_session

    return wrap_session(config, _show_fixtures_per_test)


def _show_fixtures_per_test(config: Config, session: Session) -> None:
    import _pytest.config

    session.perform_collect()
    curdir = py.path.local()
    tw = _pytest.config.create_terminal_writer(config)
    verbose = config.getvalue("verbose")

    def get_best_relpath(func):
        loc = getlocation(func, str(curdir))
        return curdir.bestrelpath(py.path.local(loc))

    def write_fixture(fixture_def: fixtures.FixtureDef[object]) -> None:
        argname = fixture_def.argname
        if verbose <= 0 and argname.startswith("_"):
            return
        if verbose > 0:
            bestrel = get_best_relpath(fixture_def.func)
            funcargspec = f"{argname} -- {bestrel}"
        else:
            funcargspec = argname
        tw.line(funcargspec, green=True)
        fixture_doc = inspect.getdoc(fixture_def.func)
        if fixture_doc:
            write_docstring(tw, fixture_doc)
        else:
            tw.line("    no docstring available", red=True)

    def write_item(item: nodes.Item) -> None:
        # Not all items have _fixtureinfo attribute.
        info: Optional[FuncFixtureInfo] = getattr(item, "_fixtureinfo", None)
        if info is None or not info.name2fixturedefs:
            # This test item does not use any fixtures.
            return
        tw.line()
        tw.sep("-", f"fixtures used by {item.name}")
        # TODO: Fix this type ignore.
        tw.sep("-", "({})".format(get_best_relpath(item.function)))  # type: ignore[attr-defined]
        # dict key not used in loop but needed for sorting.
        for _, fixturedefs in sorted(info.name2fixturedefs.items()):
            assert fixturedefs is not None
            if not fixturedefs:
                continue
            # Last item is expected to be the one used by the test item.
            write_fixture(fixturedefs[-1])

    for session_item in session.items:
        write_item(session_item)


def showfixtures(config: Config) -> Union[int, ExitCode]:
    from _pytest.main import wrap_session

    return wrap_session(config, _showfixtures_main)


def _showfixtures_main(config: Config, session: Session) -> None:
    import _pytest.config

    session.perform_collect()
    curdir = py.path.local()
    tw = _pytest.config.create_terminal_writer(config)
    verbose = config.getvalue("verbose")

    fm = session._fixturemanager

    available = []
    seen: Set[Tuple[str, str]] = set()

    for argname, fixturedefs in fm._arg2fixturedefs.items():
        assert fixturedefs is not None
        if not fixturedefs:
            continue
        for fixturedef in fixturedefs:
            loc = getlocation(fixturedef.func, str(curdir))
            if (fixturedef.argname, loc) in seen:
                continue
            seen.add((fixturedef.argname, loc))
            available.append(
                (
                    len(fixturedef.baseid),
                    fixturedef.func.__module__,
                    curdir.bestrelpath(py.path.local(loc)),
                    fixturedef.argname,
                    fixturedef,
                )
            )

    available.sort()
    currentmodule = None
    for baseid, module, bestrel, argname, fixturedef in available:
        if currentmodule != module:
            if not module.startswith("_pytest."):
                tw.line()
                tw.sep("-", f"fixtures defined from {module}")
                currentmodule = module
        if verbose <= 0 and argname[0] == "_":
            continue
        tw.write(argname, green=True)
        if fixturedef.scope != "function":
            tw.write(" [%s scope]" % fixturedef.scope, cyan=True)
        if verbose > 0:
            tw.write(" -- %s" % bestrel, yellow=True)
        tw.write("\n")
        loc = getlocation(fixturedef.func, str(curdir))
        doc = inspect.getdoc(fixturedef.func)
        if doc:
            write_docstring(tw, doc)
        else:
            tw.line(f"    {loc}: no docstring available", red=True)
        tw.line()


def write_docstring(tw: TerminalWriter, doc: str, indent: str = "    ") -> None:
    for line in doc.split("\n"):
        tw.line(indent + line)


class Function(PyobjMixin, nodes.Item):
    """An Item responsible for setting up and executing a Python test function.

    param name:
        The full function name, including any decorations like those
        added by parametrization (``my_func[my_param]``).
    param parent:
        The parent Node.
    param config:
        The pytest Config object.
    param callspec:
        If given, this is function has been parametrized and the callspec contains
        meta information about the parametrization.
    param callobj:
        If given, the object which will be called when the Function is invoked,
        otherwise the callobj will be obtained from ``parent`` using ``originalname``.
    param keywords:
        Keywords bound to the function object for "-k" matching.
    param session:
        The pytest Session object.
    param fixtureinfo:
        Fixture information already resolved at this fixture node..
    param originalname:
        The attribute name to use for accessing the underlying function object.
        Defaults to ``name``. Set this if name is different from the original name,
        for example when it contains decorations like those added by parametrization
        (``my_func[my_param]``).
    """

    # Disable since functions handle it themselves.
    _ALLOW_MARKERS = False

    def __init__(
        self,
        name: str,
        parent,
        config: Optional[Config] = None,
        callspec: Optional[CallSpec2] = None,
        callobj=NOTSET,
        keywords=None,
        session: Optional[Session] = None,
        fixtureinfo: Optional[FuncFixtureInfo] = None,
        originalname: Optional[str] = None,
    ) -> None:
        super().__init__(name, parent, config=config, session=session)

        if callobj is not NOTSET:
            self.obj = callobj

        #: Original function name, without any decorations (for example
        #: parametrization adds a ``"[...]"`` suffix to function names), used to access
        #: the underlying function object from ``parent`` (in case ``callobj`` is not given
        #: explicitly).
        #:
        #: .. versionadded:: 3.0
        self.originalname = originalname or name

        # Note: when FunctionDefinition is introduced, we should change ``originalname``
        # to a readonly property that returns FunctionDefinition.name.

        self.keywords.update(self.obj.__dict__)
        self.own_markers.extend(get_unpacked_marks(self.obj))
        if callspec:
            self.callspec = callspec
            # this is total hostile and a mess
            # keywords are broken by design by now
            # this will be redeemed later
            for mark in callspec.marks:
                # feel free to cry, this was broken for years before
                # and keywords cant fix it per design
                self.keywords[mark.name] = mark
            self.own_markers.extend(normalize_mark_list(callspec.marks))
        if keywords:
            self.keywords.update(keywords)

        # todo: this is a hell of a hack
        # https://github.com/pytest-dev/pytest/issues/4569

        self.keywords.update(
            {
                mark.name: True
                for mark in self.iter_markers()
                if mark.name not in self.keywords
            }
        )

        if fixtureinfo is None:
            fixtureinfo = self.session._fixturemanager.getfixtureinfo(
                self, self.obj, self.cls, funcargs=True
            )
        self._fixtureinfo: FuncFixtureInfo = fixtureinfo
        self.fixturenames = fixtureinfo.names_closure
        self._initrequest()

    @classmethod
    def from_parent(cls, parent, **kw):  # todo: determine sound type limitations
        """The public constructor."""
        return super().from_parent(parent=parent, **kw)

    def _initrequest(self) -> None:
        self.funcargs: Dict[str, object] = {}
        self._request = fixtures.FixtureRequest(self, _ispytest=True)

    @property
    def function(self):
        """Underlying python 'function' object."""
        return getimfunc(self.obj)

    def _getobj(self):
        assert self.parent is not None
        return getattr(self.parent.obj, self.originalname)  # type: ignore[attr-defined]

    @property
    def _pyfuncitem(self):
        """(compatonly) for code expecting pytest-2.2 style request objects."""
        return self

    def runtest(self) -> None:
        """Execute the underlying test function."""
        self.ihook.pytest_pyfunc_call(pyfuncitem=self)

    def setup(self) -> None:
        if isinstance(self.parent, Instance):
            self.parent.newinstance()
            self.obj = self._getobj()
        self._request._fillfixtures()

    def _prunetraceback(self, excinfo: ExceptionInfo[BaseException]) -> None:
        if hasattr(self, "_obj") and not self.config.getoption("fulltrace", False):
            code = _pytest._code.Code.from_function(get_real_func(self.obj))
            path, firstlineno = code.path, code.firstlineno
            traceback = excinfo.traceback
            ntraceback = traceback.cut(path=path, firstlineno=firstlineno)
            if ntraceback == traceback:
                ntraceback = ntraceback.cut(path=path)
                if ntraceback == traceback:
                    ntraceback = ntraceback.filter(filter_traceback)
                    if not ntraceback:
                        ntraceback = traceback

            excinfo.traceback = ntraceback.filter()
            # issue364: mark all but first and last frames to
            # only show a single-line message for each frame.
            if self.config.getoption("tbstyle", "auto") == "auto":
                if len(excinfo.traceback) > 2:
                    for entry in excinfo.traceback[1:-1]:
                        entry.set_repr_style("short")

    # TODO: Type ignored -- breaks Liskov Substitution.
    def repr_failure(  # type: ignore[override]
        self, excinfo: ExceptionInfo[BaseException],
    ) -> Union[str, TerminalRepr]:
        style = self.config.getoption("tbstyle", "auto")
        if style == "auto":
            style = "long"
        return self._repr_failure_py(excinfo, style=style)


class FunctionDefinition(Function):
    """
    This class is a step gap solution until we evolve to have actual function definition nodes
    and manage to get rid of ``metafunc``.
    """

    def runtest(self) -> None:
        raise RuntimeError("function definitions are not supposed to be run as tests")

    setup = runtest
