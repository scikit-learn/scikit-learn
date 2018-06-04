""" Python test discovery, setup and run of test functions. """
from __future__ import absolute_import, division, print_function

import fnmatch
import inspect
import sys
import os
import collections
import warnings
from textwrap import dedent
from itertools import count


import py
import six
from _pytest.mark import MarkerError
from _pytest.config import hookimpl

import _pytest
import pluggy
from _pytest import fixtures
from _pytest import nodes
from _pytest import deprecated
from _pytest.compat import (
    isclass, isfunction, is_generator, ascii_escaped,
    REGEX_TYPE, STRING_TYPES, NoneType, NOTSET,
    get_real_func, getfslineno, safe_getattr,
    safe_str, getlocation, enum, get_default_arg_names
)
from _pytest.outcomes import fail
from _pytest.mark.structures import transfer_markers


# relative paths that we use to filter traceback entries from appearing to the user;
# see filter_traceback
# note: if we need to add more paths than what we have now we should probably use a list
# for better maintenance
_pluggy_dir = py.path.local(pluggy.__file__.rstrip("oc"))
# pluggy is either a package or a single module depending on the version
if _pluggy_dir.basename == '__init__.py':
    _pluggy_dir = _pluggy_dir.dirpath()
_pytest_dir = py.path.local(_pytest.__file__).dirpath()
_py_dir = py.path.local(py.__file__).dirpath()


def filter_traceback(entry):
    """Return True if a TracebackEntry instance should be removed from tracebacks:
    * dynamically generated code (no code to show up for it);
    * internal traceback from pytest or its internal libraries, py and pluggy.
    """
    # entry.path might sometimes return a str object when the entry
    # points to dynamically generated code
    # see https://bitbucket.org/pytest-dev/py/issues/71
    raw_filename = entry.frame.code.raw.co_filename
    is_generated = '<' in raw_filename and '>' in raw_filename
    if is_generated:
        return False
    # entry.path might point to an non-existing file, in which case it will
    # also return a str object. see #1133
    p = py.path.local(entry.path)
    return not p.relto(_pluggy_dir) and not p.relto(_pytest_dir) and not p.relto(_py_dir)


def pyobj_property(name):
    def get(self):
        node = self.getparent(getattr(__import__('pytest'), name))
        if node is not None:
            return node.obj
    doc = "python %s object this node was collected from (can be None)." % (
          name.lower(),)
    return property(get, None, None, doc)


def pytest_addoption(parser):
    group = parser.getgroup("general")
    group.addoption('--fixtures', '--funcargs',
                    action="store_true", dest="showfixtures", default=False,
                    help="show available fixtures, sorted by plugin appearance "
                    "(fixtures with leading '_' are only shown with '-v')")
    group.addoption(
        '--fixtures-per-test',
        action="store_true",
        dest="show_fixtures_per_test",
        default=False,
        help="show fixtures per test",
    )
    parser.addini("usefixtures", type="args", default=[],
                  help="list of default fixtures to be used with this project")
    parser.addini("python_files", type="args",
                  default=['test_*.py', '*_test.py'],
                  help="glob-style file patterns for Python test module discovery")
    parser.addini("python_classes", type="args", default=["Test", ],
                  help="prefixes or glob names for Python test class discovery")
    parser.addini("python_functions", type="args", default=["test", ],
                  help="prefixes or glob names for Python test function and "
                  "method discovery")

    group.addoption("--import-mode", default="prepend",
                    choices=["prepend", "append"], dest="importmode",
                    help="prepend/append to sys.path when importing test modules, "
                    "default is to prepend.")


def pytest_cmdline_main(config):
    if config.option.showfixtures:
        showfixtures(config)
        return 0
    if config.option.show_fixtures_per_test:
        show_fixtures_per_test(config)
        return 0


def pytest_generate_tests(metafunc):
    # those alternative spellings are common - raise a specific error to alert
    # the user
    alt_spellings = ['parameterize', 'parametrise', 'parameterise']
    for attr in alt_spellings:
        if hasattr(metafunc.function, attr):
            msg = "{0} has '{1}', spelling should be 'parametrize'"
            raise MarkerError(msg.format(metafunc.function.__name__, attr))
    try:
        markers = metafunc.function.parametrize
    except AttributeError:
        return
    for marker in markers:
        metafunc.parametrize(*marker.args, **marker.kwargs)


def pytest_configure(config):
    config.addinivalue_line("markers",
                            "parametrize(argnames, argvalues): call a test function multiple "
                            "times passing in different arguments in turn. argvalues generally "
                            "needs to be a list of values if argnames specifies only one name "
                            "or a list of tuples of values if argnames specifies multiple names. "
                            "Example: @parametrize('arg1', [1,2]) would lead to two calls of the "
                            "decorated test function, one with arg1=1 and another with arg1=2."
                            "see http://pytest.org/latest/parametrize.html for more info and "
                            "examples."
                            )
    config.addinivalue_line("markers",
                            "usefixtures(fixturename1, fixturename2, ...): mark tests as needing "
                            "all of the specified fixtures. see http://pytest.org/latest/fixture.html#usefixtures "
                            )


@hookimpl(trylast=True)
def pytest_pyfunc_call(pyfuncitem):
    testfunction = pyfuncitem.obj
    if pyfuncitem._isyieldedfunction():
        testfunction(*pyfuncitem._args)
    else:
        funcargs = pyfuncitem.funcargs
        testargs = {}
        for arg in pyfuncitem._fixtureinfo.argnames:
            testargs[arg] = funcargs[arg]
        testfunction(**testargs)
    return True


def pytest_collect_file(path, parent):
    ext = path.ext
    if ext == ".py":
        if not parent.session.isinitpath(path):
            for pat in parent.config.getini('python_files'):
                if path.fnmatch(pat):
                    break
            else:
                return
        ihook = parent.session.gethookproxy(path)
        return ihook.pytest_pycollect_makemodule(path=path, parent=parent)


def pytest_pycollect_makemodule(path, parent):
    return Module(path, parent)


@hookimpl(hookwrapper=True)
def pytest_pycollect_makeitem(collector, name, obj):
    outcome = yield
    res = outcome.get_result()
    if res is not None:
        return
    # nothing was collected elsewhere, let's do it here
    if isclass(obj):
        if collector.istestclass(obj, name):
            Class = collector._getcustomclass("Class")
            outcome.force_result(Class(name, parent=collector))
    elif collector.istestfunction(obj, name):
        # mock seems to store unbound methods (issue473), normalize it
        obj = getattr(obj, "__func__", obj)
        # We need to try and unwrap the function if it's a functools.partial
        # or a funtools.wrapped.
        # We musn't if it's been wrapped with mock.patch (python 2 only)
        if not (isfunction(obj) or isfunction(get_real_func(obj))):
            collector.warn(code="C2", message="cannot collect %r because it is not a function."
                           % name, )
        elif getattr(obj, "__test__", True):
            if is_generator(obj):
                res = Generator(name, parent=collector)
            else:
                res = list(collector._genfunctions(name, obj))
            outcome.force_result(res)


def pytest_make_parametrize_id(config, val, argname=None):
    return None


class PyobjContext(object):
    module = pyobj_property("Module")
    cls = pyobj_property("Class")
    instance = pyobj_property("Instance")


class PyobjMixin(PyobjContext):
    def obj():
        def fget(self):
            obj = getattr(self, '_obj', None)
            if obj is None:
                self._obj = obj = self._getobj()
            return obj

        def fset(self, value):
            self._obj = value

        return property(fget, fset, None, "underlying python object")

    obj = obj()

    def _getobj(self):
        return getattr(self.parent.obj, self.name)

    def getmodpath(self, stopatmodule=True, includemodule=False):
        """ return python path relative to the containing module. """
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
        s = ".".join(parts)
        return s.replace(".[", "[")

    def _getfslineno(self):
        return getfslineno(self.obj)

    def reportinfo(self):
        # XXX caching?
        obj = self.obj
        compat_co_firstlineno = getattr(obj, 'compat_co_firstlineno', None)
        if isinstance(compat_co_firstlineno, int):
            # nose compatibility
            fspath = sys.modules[obj.__module__].__file__
            if fspath.endswith(".pyc"):
                fspath = fspath[:-1]
            lineno = compat_co_firstlineno
        else:
            fspath, lineno = getfslineno(obj)
        modpath = self.getmodpath()
        assert isinstance(lineno, int)
        return fspath, lineno, modpath


class PyCollector(PyobjMixin, nodes.Collector):

    def funcnamefilter(self, name):
        return self._matches_prefix_or_glob_option('python_functions', name)

    def isnosetest(self, obj):
        """ Look for the __test__ attribute, which is applied by the
        @nose.tools.istest decorator
        """
        # We explicitly check for "is True" here to not mistakenly treat
        # classes with a custom __getattr__ returning something truthy (like a
        # function) as test classes.
        return safe_getattr(obj, '__test__', False) is True

    def classnamefilter(self, name):
        return self._matches_prefix_or_glob_option('python_classes', name)

    def istestfunction(self, obj, name):
        if self.funcnamefilter(name) or self.isnosetest(obj):
            if isinstance(obj, staticmethod):
                # static methods need to be unwrapped
                obj = safe_getattr(obj, '__func__', False)
                if obj is False:
                    # Python 2.6 wraps in a different way that we won't try to handle
                    msg = "cannot collect static method %r because " \
                          "it is not a function (always the case in Python 2.6)"
                    self.warn(
                        code="C2", message=msg % name)
                    return False
            return (
                safe_getattr(obj, "__call__", False) and fixtures.getfixturemarker(obj) is None
            )
        else:
            return False

    def istestclass(self, obj, name):
        return self.classnamefilter(name) or self.isnosetest(obj)

    def _matches_prefix_or_glob_option(self, option_name, name):
        """
        checks if the given name matches the prefix or glob-pattern defined
        in ini configuration.
        """
        for option in self.config.getini(option_name):
            if name.startswith(option):
                return True
            # check that name looks like a glob-string before calling fnmatch
            # because this is called for every name in each collected module,
            # and fnmatch is somewhat expensive to call
            elif ('*' in option or '?' in option or '[' in option) and \
                    fnmatch.fnmatch(name, option):
                return True
        return False

    def collect(self):
        if not getattr(self.obj, "__test__", True):
            return []

        # NB. we avoid random getattrs and peek in the __dict__ instead
        # (XXX originally introduced from a PyPy need, still true?)
        dicts = [getattr(self.obj, '__dict__', {})]
        for basecls in inspect.getmro(self.obj.__class__):
            dicts.append(basecls.__dict__)
        seen = {}
        values = []
        for dic in dicts:
            for name, obj in list(dic.items()):
                if name in seen:
                    continue
                seen[name] = True
                res = self._makeitem(name, obj)
                if res is None:
                    continue
                if not isinstance(res, list):
                    res = [res]
                values.extend(res)
        values.sort(key=lambda item: item.reportinfo()[:2])
        return values

    def makeitem(self, name, obj):
        warnings.warn(deprecated.COLLECTOR_MAKEITEM, stacklevel=2)
        self._makeitem(name, obj)

    def _makeitem(self, name, obj):
        # assert self.ihook.fspath == self.fspath, self
        return self.ihook.pytest_pycollect_makeitem(
            collector=self, name=name, obj=obj)

    def _genfunctions(self, name, funcobj):
        module = self.getparent(Module).obj
        clscol = self.getparent(Class)
        cls = clscol and clscol.obj or None
        transfer_markers(funcobj, cls, module)
        fm = self.session._fixturemanager
        fixtureinfo = fm.getfixtureinfo(self, funcobj, cls)
        metafunc = Metafunc(funcobj, fixtureinfo, self.config,
                            cls=cls, module=module)
        methods = []
        if hasattr(module, "pytest_generate_tests"):
            methods.append(module.pytest_generate_tests)
        if hasattr(cls, "pytest_generate_tests"):
            methods.append(cls().pytest_generate_tests)
        if methods:
            self.ihook.pytest_generate_tests.call_extra(methods,
                                                        dict(metafunc=metafunc))
        else:
            self.ihook.pytest_generate_tests(metafunc=metafunc)

        Function = self._getcustomclass("Function")
        if not metafunc._calls:
            yield Function(name, parent=self, fixtureinfo=fixtureinfo)
        else:
            # add funcargs() as fixturedefs to fixtureinfo.arg2fixturedefs
            fixtures.add_funcarg_pseudo_fixture_def(self, metafunc, fm)

            for callspec in metafunc._calls:
                subname = "%s[%s]" % (name, callspec.id)
                yield Function(name=subname, parent=self,
                               callspec=callspec, callobj=funcobj,
                               fixtureinfo=fixtureinfo,
                               keywords={callspec.id: True},
                               originalname=name,
                               )


class Module(nodes.File, PyCollector):
    """ Collector for test classes and functions. """

    def _getobj(self):
        return self._importtestmodule()

    def collect(self):
        self.session._fixturemanager.parsefactories(self)
        return super(Module, self).collect()

    def _importtestmodule(self):
        # we assume we are only called once per module
        importmode = self.config.getoption("--import-mode")
        try:
            mod = self.fspath.pyimport(ensuresyspath=importmode)
        except SyntaxError:
            raise self.CollectError(
                _pytest._code.ExceptionInfo().getrepr(style="short"))
        except self.fspath.ImportMismatchError:
            e = sys.exc_info()[1]
            raise self.CollectError(
                "import file mismatch:\n"
                "imported module %r has this __file__ attribute:\n"
                "  %s\n"
                "which is not the same as the test file we want to collect:\n"
                "  %s\n"
                "HINT: remove __pycache__ / .pyc files and/or use a "
                "unique basename for your test file modules"
                % e.args
            )
        except ImportError:
            from _pytest._code.code import ExceptionInfo
            exc_info = ExceptionInfo()
            if self.config.getoption('verbose') < 2:
                exc_info.traceback = exc_info.traceback.filter(filter_traceback)
            exc_repr = exc_info.getrepr(style='short') if exc_info.traceback else exc_info.exconly()
            formatted_tb = safe_str(exc_repr)
            raise self.CollectError(
                "ImportError while importing test module '{fspath}'.\n"
                "Hint: make sure your test modules/packages have valid Python names.\n"
                "Traceback:\n"
                "{traceback}".format(fspath=self.fspath, traceback=formatted_tb)
            )
        except _pytest.runner.Skipped as e:
            if e.allow_module_level:
                raise
            raise self.CollectError(
                "Using pytest.skip outside of a test is not allowed. "
                "To decorate a test function, use the @pytest.mark.skip "
                "or @pytest.mark.skipif decorators instead, and to skip a "
                "module use `pytestmark = pytest.mark.{skip,skipif}."
            )
        self.config.pluginmanager.consider_module(mod)
        return mod

    def setup(self):
        setup_module = _get_xunit_setup_teardown(self.obj, "setUpModule")
        if setup_module is None:
            setup_module = _get_xunit_setup_teardown(self.obj, "setup_module")
        if setup_module is not None:
            setup_module()

        teardown_module = _get_xunit_setup_teardown(self.obj, 'tearDownModule')
        if teardown_module is None:
            teardown_module = _get_xunit_setup_teardown(self.obj, 'teardown_module')
        if teardown_module is not None:
            self.addfinalizer(teardown_module)


def _get_xunit_setup_teardown(holder, attr_name, param_obj=None):
    """
    Return a callable to perform xunit-style setup or teardown if
    the function exists in the ``holder`` object.
    The ``param_obj`` parameter is the parameter which will be passed to the function
    when the callable is called without arguments, defaults to the ``holder`` object.
    Return ``None`` if a suitable callable is not found.
    """
    param_obj = param_obj if param_obj is not None else holder
    result = _get_xunit_func(holder, attr_name)
    if result is not None:
        arg_count = result.__code__.co_argcount
        if inspect.ismethod(result):
            arg_count -= 1
        if arg_count:
            return lambda: result(param_obj)
        else:
            return result


def _get_xunit_func(obj, name):
    """Return the attribute from the given object to be used as a setup/teardown
    xunit-style function, but only if not marked as a fixture to
    avoid calling it twice.
    """
    meth = getattr(obj, name, None)
    if fixtures.getfixturemarker(meth) is None:
        return meth


class Class(PyCollector):
    """ Collector for test methods. """

    def collect(self):
        if not safe_getattr(self.obj, "__test__", True):
            return []
        if hasinit(self.obj):
            self.warn("C1", "cannot collect test class %r because it has a "
                      "__init__ constructor" % self.obj.__name__)
            return []
        elif hasnew(self.obj):
            self.warn("C1", "cannot collect test class %r because it has a "
                            "__new__ constructor" % self.obj.__name__)
            return []
        return [self._getcustomclass("Instance")(name="()", parent=self)]

    def setup(self):
        setup_class = _get_xunit_func(self.obj, 'setup_class')
        if setup_class is not None:
            setup_class = getattr(setup_class, 'im_func', setup_class)
            setup_class = getattr(setup_class, '__func__', setup_class)
            setup_class(self.obj)

        fin_class = getattr(self.obj, 'teardown_class', None)
        if fin_class is not None:
            fin_class = getattr(fin_class, 'im_func', fin_class)
            fin_class = getattr(fin_class, '__func__', fin_class)
            self.addfinalizer(lambda: fin_class(self.obj))


class Instance(PyCollector):
    def _getobj(self):
        return self.parent.obj()

    def collect(self):
        self.session._fixturemanager.parsefactories(self)
        return super(Instance, self).collect()

    def newinstance(self):
        self.obj = self._getobj()
        return self.obj


class FunctionMixin(PyobjMixin):
    """ mixin for the code common to Function and Generator.
    """

    def setup(self):
        """ perform setup for this test function. """
        if hasattr(self, '_preservedparent'):
            obj = self._preservedparent
        elif isinstance(self.parent, Instance):
            obj = self.parent.newinstance()
            self.obj = self._getobj()
        else:
            obj = self.parent.obj
        if inspect.ismethod(self.obj):
            setup_name = 'setup_method'
            teardown_name = 'teardown_method'
        else:
            setup_name = 'setup_function'
            teardown_name = 'teardown_function'
        setup_func_or_method = _get_xunit_setup_teardown(obj, setup_name, param_obj=self.obj)
        if setup_func_or_method is not None:
            setup_func_or_method()
        teardown_func_or_method = _get_xunit_setup_teardown(obj, teardown_name, param_obj=self.obj)
        if teardown_func_or_method is not None:
            self.addfinalizer(teardown_func_or_method)

    def _prunetraceback(self, excinfo):
        if hasattr(self, '_obj') and not self.config.option.fulltrace:
            code = _pytest._code.Code(get_real_func(self.obj))
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
            # only show a single-line message for each frame
            if self.config.option.tbstyle == "auto":
                if len(excinfo.traceback) > 2:
                    for entry in excinfo.traceback[1:-1]:
                        entry.set_repr_style('short')

    def _repr_failure_py(self, excinfo, style="long"):
        if excinfo.errisinstance(fail.Exception):
            if not excinfo.value.pytrace:
                return py._builtin._totext(excinfo.value)
        return super(FunctionMixin, self)._repr_failure_py(excinfo,
                                                           style=style)

    def repr_failure(self, excinfo, outerr=None):
        assert outerr is None, "XXX outerr usage is deprecated"
        style = self.config.option.tbstyle
        if style == "auto":
            style = "long"
        return self._repr_failure_py(excinfo, style=style)


class Generator(FunctionMixin, PyCollector):
    def collect(self):
        # test generators are seen as collectors but they also
        # invoke setup/teardown on popular request
        # (induced by the common "test_*" naming shared with normal tests)
        from _pytest import deprecated
        self.session._setupstate.prepare(self)
        # see FunctionMixin.setup and test_setupstate_is_preserved_134
        self._preservedparent = self.parent.obj
        values = []
        seen = {}
        for i, x in enumerate(self.obj()):
            name, call, args = self.getcallargs(x)
            if not callable(call):
                raise TypeError("%r yielded non callable test %r" % (self.obj, call,))
            if name is None:
                name = "[%d]" % i
            else:
                name = "['%s']" % name
            if name in seen:
                raise ValueError("%r generated tests with non-unique name %r" % (self, name))
            seen[name] = True
            values.append(self.Function(name, self, args=args, callobj=call))
        self.warn('C1', deprecated.YIELD_TESTS)
        return values

    def getcallargs(self, obj):
        if not isinstance(obj, (tuple, list)):
            obj = (obj,)
        # explicit naming
        if isinstance(obj[0], six.string_types):
            name = obj[0]
            obj = obj[1:]
        else:
            name = None
        call, args = obj[0], obj[1:]
        return name, call, args


def hasinit(obj):
    init = getattr(obj, '__init__', None)
    if init:
        return init != object.__init__


def hasnew(obj):
    new = getattr(obj, '__new__', None)
    if new:
        return new != object.__new__


class CallSpec2(object):
    def __init__(self, metafunc):
        self.metafunc = metafunc
        self.funcargs = {}
        self._idlist = []
        self.params = {}
        self._globalid = NOTSET
        self._globalid_args = set()
        self._globalparam = NOTSET
        self._arg2scopenum = {}  # used for sorting parametrized resources
        self.marks = []
        self.indices = {}

    def copy(self, metafunc):
        cs = CallSpec2(self.metafunc)
        cs.funcargs.update(self.funcargs)
        cs.params.update(self.params)
        cs.marks.extend(self.marks)
        cs.indices.update(self.indices)
        cs._arg2scopenum.update(self._arg2scopenum)
        cs._idlist = list(self._idlist)
        cs._globalid = self._globalid
        cs._globalid_args = self._globalid_args
        cs._globalparam = self._globalparam
        return cs

    def _checkargnotcontained(self, arg):
        if arg in self.params or arg in self.funcargs:
            raise ValueError("duplicate %r" % (arg,))

    def getparam(self, name):
        try:
            return self.params[name]
        except KeyError:
            if self._globalparam is NOTSET:
                raise ValueError(name)
            return self._globalparam

    @property
    def id(self):
        return "-".join(map(str, filter(None, self._idlist)))

    def setmulti2(self, valtypes, argnames, valset, id, marks, scopenum,
                  param_index):
        for arg, val in zip(argnames, valset):
            self._checkargnotcontained(arg)
            valtype_for_arg = valtypes[arg]
            getattr(self, valtype_for_arg)[arg] = val
            self.indices[arg] = param_index
            self._arg2scopenum[arg] = scopenum
        self._idlist.append(id)
        self.marks.extend(marks)

    def setall(self, funcargs, id, param):
        for x in funcargs:
            self._checkargnotcontained(x)
        self.funcargs.update(funcargs)
        if id is not NOTSET:
            self._idlist.append(id)
        if param is not NOTSET:
            assert self._globalparam is NOTSET
            self._globalparam = param
        for arg in funcargs:
            self._arg2scopenum[arg] = fixtures.scopenum_function


class Metafunc(fixtures.FuncargnamesCompatAttr):
    """
    Metafunc objects are passed to the :func:`pytest_generate_tests <_pytest.hookspec.pytest_generate_tests>` hook.
    They help to inspect a test function and to generate tests according to
    test configuration or values specified in the class or module where a
    test function is defined.
    """

    def __init__(self, function, fixtureinfo, config, cls=None, module=None):
        #: access to the :class:`_pytest.config.Config` object for the test session
        self.config = config

        #: the module object where the test function is defined in.
        self.module = module

        #: underlying python test function
        self.function = function

        #: set of fixture names required by the test function
        self.fixturenames = fixtureinfo.names_closure

        #: class object where the test function is defined in or ``None``.
        self.cls = cls

        self._calls = []
        self._ids = set()
        self._arg2fixturedefs = fixtureinfo.name2fixturedefs

    def parametrize(self, argnames, argvalues, indirect=False, ids=None,
                    scope=None):
        """ Add new invocations to the underlying test function using the list
        of argvalues for the given argnames.  Parametrization is performed
        during the collection phase.  If you need to setup expensive resources
        see about setting indirect to do it rather at test setup time.

        :arg argnames: a comma-separated string denoting one or more argument
                       names, or a list/tuple of argument strings.

        :arg argvalues: The list of argvalues determines how often a
            test is invoked with different argument values.  If only one
            argname was specified argvalues is a list of values.  If N
            argnames were specified, argvalues must be a list of N-tuples,
            where each tuple-element specifies a value for its respective
            argname.

        :arg indirect: The list of argnames or boolean. A list of arguments'
            names (subset of argnames). If True the list contains all names from
            the argnames. Each argvalue corresponding to an argname in this list will
            be passed as request.param to its respective argname fixture
            function so that it can perform more expensive setups during the
            setup phase of a test rather than at collection time.

        :arg ids: list of string ids, or a callable.
            If strings, each is corresponding to the argvalues so that they are
            part of the test id. If None is given as id of specific test, the
            automatically generated id for that argument will be used.
            If callable, it should take one argument (a single argvalue) and return
            a string or return None. If None, the automatically generated id for that
            argument will be used.
            If no ids are provided they will be generated automatically from
            the argvalues.

        :arg scope: if specified it denotes the scope of the parameters.
            The scope is used for grouping tests by parameter instances.
            It will also override any fixture-function defined scope, allowing
            to set a dynamic scope using test context or configuration.
        """
        from _pytest.fixtures import scope2index
        from _pytest.mark import ParameterSet
        from py.io import saferepr

        argnames, parameters = ParameterSet._for_parametrize(
            argnames, argvalues, self.function, self.config)
        del argvalues
        default_arg_names = set(get_default_arg_names(self.function))

        if scope is None:
            scope = _find_parametrized_scope(argnames, self._arg2fixturedefs, indirect)

        scopenum = scope2index(scope, descr='call to {0}'.format(self.parametrize))
        valtypes = {}
        for arg in argnames:
            if arg not in self.fixturenames:
                if arg in default_arg_names:
                    raise ValueError("%r already takes an argument %r with a default value" % (self.function, arg))
                else:
                    if isinstance(indirect, (tuple, list)):
                        name = 'fixture' if arg in indirect else 'argument'
                    else:
                        name = 'fixture' if indirect else 'argument'
                    raise ValueError(
                        "%r uses no %s %r" % (
                            self.function, name, arg))

        if indirect is True:
            valtypes = dict.fromkeys(argnames, "params")
        elif indirect is False:
            valtypes = dict.fromkeys(argnames, "funcargs")
        elif isinstance(indirect, (tuple, list)):
            valtypes = dict.fromkeys(argnames, "funcargs")
            for arg in indirect:
                if arg not in argnames:
                    raise ValueError("indirect given to %r: fixture %r doesn't exist" % (
                                     self.function, arg))
                valtypes[arg] = "params"
        idfn = None
        if callable(ids):
            idfn = ids
            ids = None
        if ids:
            if len(ids) != len(parameters):
                raise ValueError('%d tests specified with %d ids' % (
                                 len(parameters), len(ids)))
            for id_value in ids:
                if id_value is not None and not isinstance(id_value, six.string_types):
                    msg = 'ids must be list of strings, found: %s (type: %s)'
                    raise ValueError(msg % (saferepr(id_value), type(id_value).__name__))
        ids = idmaker(argnames, parameters, idfn, ids, self.config)
        newcalls = []
        for callspec in self._calls or [CallSpec2(self)]:
            elements = zip(ids, parameters, count())
            for a_id, param, param_index in elements:
                if len(param.values) != len(argnames):
                    raise ValueError(
                        'In "parametrize" the number of values ({0}) must be '
                        'equal to the number of names ({1})'.format(
                            param.values, argnames))
                newcallspec = callspec.copy(self)
                newcallspec.setmulti2(valtypes, argnames, param.values, a_id,
                                      param.marks, scopenum, param_index)
                newcalls.append(newcallspec)
        self._calls = newcalls

    def addcall(self, funcargs=None, id=NOTSET, param=NOTSET):
        """ Add a new call to the underlying test function during the collection phase of a test run.

        .. deprecated:: 3.3

            Use :meth:`parametrize` instead.

        Note that request.addcall() is called during the test collection phase prior and
        independently to actual test execution.  You should only use addcall()
        if you need to specify multiple arguments of a test function.

        :arg funcargs: argument keyword dictionary used when invoking
            the test function.

        :arg id: used for reporting and identification purposes.  If you
            don't supply an `id` an automatic unique id will be generated.

        :arg param: a parameter which will be exposed to a later fixture function
            invocation through the ``request.param`` attribute.
        """
        if self.config:
            self.config.warn('C1', message=deprecated.METAFUNC_ADD_CALL, fslocation=None)
        assert funcargs is None or isinstance(funcargs, dict)
        if funcargs is not None:
            for name in funcargs:
                if name not in self.fixturenames:
                    fail("funcarg %r not used in this function." % name)
        else:
            funcargs = {}
        if id is None:
            raise ValueError("id=None not allowed")
        if id is NOTSET:
            id = len(self._calls)
        id = str(id)
        if id in self._ids:
            raise ValueError("duplicate id %r" % id)
        self._ids.add(id)

        cs = CallSpec2(self)
        cs.setall(funcargs, id, param)
        self._calls.append(cs)


def _find_parametrized_scope(argnames, arg2fixturedefs, indirect):
    """Find the most appropriate scope for a parametrized call based on its arguments.

    When there's at least one direct argument, always use "function" scope.

    When a test function is parametrized and all its arguments are indirect
    (e.g. fixtures), return the most narrow scope based on the fixtures used.

    Related to issue #1832, based on code posted by @Kingdread.
    """
    from _pytest.fixtures import scopes
    indirect_as_list = isinstance(indirect, (list, tuple))
    all_arguments_are_fixtures = indirect is True or \
        indirect_as_list and len(indirect) == argnames
    if all_arguments_are_fixtures:
        fixturedefs = arg2fixturedefs or {}
        used_scopes = [fixturedef[0].scope for name, fixturedef in fixturedefs.items()]
        if used_scopes:
            # Takes the most narrow scope from used fixtures
            for scope in reversed(scopes):
                if scope in used_scopes:
                    return scope

    return 'function'


def _idval(val, argname, idx, idfn, config=None):
    if idfn:
        s = None
        try:
            s = idfn(val)
        except Exception:
            # See issue https://github.com/pytest-dev/pytest/issues/2169
            import warnings
            msg = "Raised while trying to determine id of parameter %s at position %d." % (argname, idx)
            msg += '\nUpdate your code as this will raise an error in pytest-4.0.'
            warnings.warn(msg, DeprecationWarning)
        if s:
            return ascii_escaped(s)

    if config:
        hook_id = config.hook.pytest_make_parametrize_id(
            config=config, val=val, argname=argname)
        if hook_id:
            return hook_id

    if isinstance(val, STRING_TYPES):
        return ascii_escaped(val)
    elif isinstance(val, (float, int, bool, NoneType)):
        return str(val)
    elif isinstance(val, REGEX_TYPE):
        return ascii_escaped(val.pattern)
    elif enum is not None and isinstance(val, enum.Enum):
        return str(val)
    elif (isclass(val) or isfunction(val)) and hasattr(val, '__name__'):
        return val.__name__
    return str(argname) + str(idx)


def _idvalset(idx, parameterset, argnames, idfn, ids, config=None):
    if parameterset.id is not None:
        return parameterset.id
    if ids is None or (idx >= len(ids) or ids[idx] is None):
        this_id = [_idval(val, argname, idx, idfn, config)
                   for val, argname in zip(parameterset.values, argnames)]
        return "-".join(this_id)
    else:
        return ascii_escaped(ids[idx])


def idmaker(argnames, parametersets, idfn=None, ids=None, config=None):
    ids = [_idvalset(valindex, parameterset, argnames, idfn, ids, config)
           for valindex, parameterset in enumerate(parametersets)]
    if len(set(ids)) != len(ids):
        # The ids are not unique
        duplicates = [testid for testid in ids if ids.count(testid) > 1]
        counters = collections.defaultdict(lambda: 0)
        for index, testid in enumerate(ids):
            if testid in duplicates:
                ids[index] = testid + str(counters[testid])
                counters[testid] += 1
    return ids


def show_fixtures_per_test(config):
    from _pytest.main import wrap_session
    return wrap_session(config, _show_fixtures_per_test)


def _show_fixtures_per_test(config, session):
    import _pytest.config
    session.perform_collect()
    curdir = py.path.local()
    tw = _pytest.config.create_terminal_writer(config)
    verbose = config.getvalue("verbose")

    def get_best_relpath(func):
        loc = getlocation(func, curdir)
        return curdir.bestrelpath(loc)

    def write_fixture(fixture_def):
        argname = fixture_def.argname
        if verbose <= 0 and argname.startswith("_"):
            return
        if verbose > 0:
            bestrel = get_best_relpath(fixture_def.func)
            funcargspec = "{0} -- {1}".format(argname, bestrel)
        else:
            funcargspec = argname
        tw.line(funcargspec, green=True)
        fixture_doc = fixture_def.func.__doc__
        if fixture_doc:
            write_docstring(tw, fixture_doc)
        else:
            tw.line('    no docstring available', red=True)

    def write_item(item):
        try:
            info = item._fixtureinfo
        except AttributeError:
            # doctests items have no _fixtureinfo attribute
            return
        if not info.name2fixturedefs:
            # this test item does not use any fixtures
            return
        tw.line()
        tw.sep('-', 'fixtures used by {0}'.format(item.name))
        tw.sep('-', '({0})'.format(get_best_relpath(item.function)))
        # dict key not used in loop but needed for sorting
        for _, fixturedefs in sorted(info.name2fixturedefs.items()):
            assert fixturedefs is not None
            if not fixturedefs:
                continue
            # last item is expected to be the one used by the test item
            write_fixture(fixturedefs[-1])

    for session_item in session.items:
        write_item(session_item)


def showfixtures(config):
    from _pytest.main import wrap_session
    return wrap_session(config, _showfixtures_main)


def _showfixtures_main(config, session):
    import _pytest.config
    session.perform_collect()
    curdir = py.path.local()
    tw = _pytest.config.create_terminal_writer(config)
    verbose = config.getvalue("verbose")

    fm = session._fixturemanager

    available = []
    seen = set()

    for argname, fixturedefs in fm._arg2fixturedefs.items():
        assert fixturedefs is not None
        if not fixturedefs:
            continue
        for fixturedef in fixturedefs:
            loc = getlocation(fixturedef.func, curdir)
            if (fixturedef.argname, loc) in seen:
                continue
            seen.add((fixturedef.argname, loc))
            available.append((len(fixturedef.baseid),
                              fixturedef.func.__module__,
                              curdir.bestrelpath(loc),
                              fixturedef.argname, fixturedef))

    available.sort()
    currentmodule = None
    for baseid, module, bestrel, argname, fixturedef in available:
        if currentmodule != module:
            if not module.startswith("_pytest."):
                tw.line()
                tw.sep("-", "fixtures defined from %s" % (module,))
                currentmodule = module
        if verbose <= 0 and argname[0] == "_":
            continue
        if verbose > 0:
            funcargspec = "%s -- %s" % (argname, bestrel,)
        else:
            funcargspec = argname
        tw.line(funcargspec, green=True)
        loc = getlocation(fixturedef.func, curdir)
        doc = fixturedef.func.__doc__ or ""
        if doc:
            write_docstring(tw, doc)
        else:
            tw.line("    %s: no docstring available" % (loc,),
                    red=True)


def write_docstring(tw, doc):
    INDENT = "    "
    doc = doc.rstrip()
    if "\n" in doc:
        firstline, rest = doc.split("\n", 1)
    else:
        firstline, rest = doc, ""

    if firstline.strip():
        tw.line(INDENT + firstline.strip())

    if rest:
        for line in dedent(rest).split("\n"):
            tw.write(INDENT + line + "\n")


class Function(FunctionMixin, nodes.Item, fixtures.FuncargnamesCompatAttr):
    """ a Function Item is responsible for setting up and executing a
    Python test function.
    """
    _genid = None

    def __init__(self, name, parent, args=None, config=None,
                 callspec=None, callobj=NOTSET, keywords=None, session=None,
                 fixtureinfo=None, originalname=None):
        super(Function, self).__init__(name, parent, config=config,
                                       session=session)
        self._args = args
        if callobj is not NOTSET:
            self.obj = callobj

        self.keywords.update(self.obj.__dict__)
        if callspec:
            self.callspec = callspec
            # this is total hostile and a mess
            # keywords are broken by design by now
            # this will be redeemed later
            for mark in callspec.marks:
                # feel free to cry, this was broken for years before
                # and keywords cant fix it per design
                self.keywords[mark.name] = mark
        if keywords:
            self.keywords.update(keywords)

        if fixtureinfo is None:
            fixtureinfo = self.session._fixturemanager.getfixtureinfo(
                self.parent, self.obj, self.cls,
                funcargs=not self._isyieldedfunction())
        self._fixtureinfo = fixtureinfo
        self.fixturenames = fixtureinfo.names_closure
        self._initrequest()

        #: original function name, without any decorations (for example
        #: parametrization adds a ``"[...]"`` suffix to function names).
        #:
        #: .. versionadded:: 3.0
        self.originalname = originalname

    def _initrequest(self):
        self.funcargs = {}
        if self._isyieldedfunction():
            assert not hasattr(self, "callspec"), (
                "yielded functions (deprecated) cannot have funcargs")
        else:
            if hasattr(self, "callspec"):
                callspec = self.callspec
                assert not callspec.funcargs
                self._genid = callspec.id
                if hasattr(callspec, "param"):
                    self.param = callspec.param
        self._request = fixtures.FixtureRequest(self)

    @property
    def function(self):
        "underlying python 'function' object"
        return getattr(self.obj, 'im_func', self.obj)

    def _getobj(self):
        name = self.name
        i = name.find("[")  # parametrization
        if i != -1:
            name = name[:i]
        return getattr(self.parent.obj, name)

    @property
    def _pyfuncitem(self):
        "(compatonly) for code expecting pytest-2.2 style request objects"
        return self

    def _isyieldedfunction(self):
        return getattr(self, "_args", None) is not None

    def runtest(self):
        """ execute the underlying test function. """
        self.ihook.pytest_pyfunc_call(pyfuncitem=self)

    def setup(self):
        super(Function, self).setup()
        fixtures.fillfixtures(self)
