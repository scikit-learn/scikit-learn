""" core implementation of testing process: init, session, runtest loop. """
from __future__ import absolute_import, division, print_function

import contextlib
import functools
import os
import pkgutil
import six
import sys

import _pytest
from _pytest import nodes
import _pytest._code
import py

from _pytest.config import directory_arg, UsageError, hookimpl
from _pytest.outcomes import exit
from _pytest.runner import collect_one_node


# exitcodes for the command line
EXIT_OK = 0
EXIT_TESTSFAILED = 1
EXIT_INTERRUPTED = 2
EXIT_INTERNALERROR = 3
EXIT_USAGEERROR = 4
EXIT_NOTESTSCOLLECTED = 5


def pytest_addoption(parser):
    parser.addini("norecursedirs", "directory patterns to avoid for recursion",
                  type="args", default=['.*', 'build', 'dist', 'CVS', '_darcs', '{arch}', '*.egg', 'venv'])
    parser.addini("testpaths", "directories to search for tests when no files or directories are given in the "
                               "command line.",
                  type="args", default=[])
    # parser.addini("dirpatterns",
    #    "patterns specifying possible locations of test files",
    #    type="linelist", default=["**/test_*.txt",
    #            "**/test_*.py", "**/*_test.py"]
    # )
    group = parser.getgroup("general", "running and selection options")
    group._addoption('-x', '--exitfirst', action="store_const",
                     dest="maxfail", const=1,
                     help="exit instantly on first error or failed test."),
    group._addoption('--maxfail', metavar="num",
                     action="store", type=int, dest="maxfail", default=0,
                     help="exit after first num failures or errors.")
    group._addoption('--strict', action="store_true",
                     help="marks not registered in configuration file raise errors.")
    group._addoption("-c", metavar="file", type=str, dest="inifilename",
                     help="load configuration from `file` instead of trying to locate one of the implicit "
                          "configuration files.")
    group._addoption("--continue-on-collection-errors", action="store_true",
                     default=False, dest="continue_on_collection_errors",
                     help="Force test execution even if collection errors occur.")
    group._addoption("--rootdir", action="store",
                     dest="rootdir",
                     help="Define root directory for tests. Can be relative path: 'root_dir', './root_dir', "
                          "'root_dir/another_dir/'; absolute path: '/home/user/root_dir'; path with variables: "
                          "'$HOME/root_dir'.")

    group = parser.getgroup("collect", "collection")
    group.addoption('--collectonly', '--collect-only', action="store_true",
                    help="only collect tests, don't execute them."),
    group.addoption('--pyargs', action="store_true",
                    help="try to interpret all arguments as python packages.")
    group.addoption("--ignore", action="append", metavar="path",
                    help="ignore path during collection (multi-allowed).")
    group.addoption("--deselect", action="append", metavar="nodeid_prefix",
                    help="deselect item during collection (multi-allowed).")
    # when changing this to --conf-cut-dir, config.py Conftest.setinitial
    # needs upgrading as well
    group.addoption('--confcutdir', dest="confcutdir", default=None,
                    metavar="dir", type=functools.partial(directory_arg, optname="--confcutdir"),
                    help="only load conftest.py's relative to specified dir.")
    group.addoption('--noconftest', action="store_true",
                    dest="noconftest", default=False,
                    help="Don't load any conftest.py files.")
    group.addoption('--keepduplicates', '--keep-duplicates', action="store_true",
                    dest="keepduplicates", default=False,
                    help="Keep duplicate tests.")
    group.addoption('--collect-in-virtualenv', action='store_true',
                    dest='collect_in_virtualenv', default=False,
                    help="Don't ignore tests in a local virtualenv directory")

    group = parser.getgroup("debugconfig",
                            "test session debugging and configuration")
    group.addoption('--basetemp', dest="basetemp", default=None, metavar="dir",
                    help="base temporary directory for this test run.")


def pytest_configure(config):
    __import__('pytest').config = config  # compatibiltiy


def wrap_session(config, doit):
    """Skeleton command line program"""
    session = Session(config)
    session.exitstatus = EXIT_OK
    initstate = 0
    try:
        try:
            config._do_configure()
            initstate = 1
            config.hook.pytest_sessionstart(session=session)
            initstate = 2
            session.exitstatus = doit(config, session) or 0
        except UsageError:
            raise
        except Failed:
            session.exitstatus = EXIT_TESTSFAILED
        except KeyboardInterrupt:
            excinfo = _pytest._code.ExceptionInfo()
            if initstate < 2 and isinstance(excinfo.value, exit.Exception):
                sys.stderr.write('{0}: {1}\n'.format(
                    excinfo.typename, excinfo.value.msg))
            config.hook.pytest_keyboard_interrupt(excinfo=excinfo)
            session.exitstatus = EXIT_INTERRUPTED
        except:  # noqa
            excinfo = _pytest._code.ExceptionInfo()
            config.notify_exception(excinfo, config.option)
            session.exitstatus = EXIT_INTERNALERROR
            if excinfo.errisinstance(SystemExit):
                sys.stderr.write("mainloop: caught Spurious SystemExit!\n")

    finally:
        excinfo = None  # Explicitly break reference cycle.
        session.startdir.chdir()
        if initstate >= 2:
            config.hook.pytest_sessionfinish(
                session=session,
                exitstatus=session.exitstatus)
        config._ensure_unconfigure()
    return session.exitstatus


def pytest_cmdline_main(config):
    return wrap_session(config, _main)


def _main(config, session):
    """ default command line protocol for initialization, session,
    running tests and reporting. """
    config.hook.pytest_collection(session=session)
    config.hook.pytest_runtestloop(session=session)

    if session.testsfailed:
        return EXIT_TESTSFAILED
    elif session.testscollected == 0:
        return EXIT_NOTESTSCOLLECTED


def pytest_collection(session):
    return session.perform_collect()


def pytest_runtestloop(session):
    if (session.testsfailed and
            not session.config.option.continue_on_collection_errors):
        raise session.Interrupted(
            "%d errors during collection" % session.testsfailed)

    if session.config.option.collectonly:
        return True

    for i, item in enumerate(session.items):
        nextitem = session.items[i + 1] if i + 1 < len(session.items) else None
        item.config.hook.pytest_runtest_protocol(item=item, nextitem=nextitem)
        if session.shouldfail:
            raise session.Failed(session.shouldfail)
        if session.shouldstop:
            raise session.Interrupted(session.shouldstop)
    return True


def _in_venv(path):
    """Attempts to detect if ``path`` is the root of a Virtual Environment by
    checking for the existence of the appropriate activate script"""
    bindir = path.join('Scripts' if sys.platform.startswith('win') else 'bin')
    if not bindir.isdir():
        return False
    activates = ('activate', 'activate.csh', 'activate.fish',
                 'Activate', 'Activate.bat', 'Activate.ps1')
    return any([fname.basename in activates for fname in bindir.listdir()])


def pytest_ignore_collect(path, config):
    ignore_paths = config._getconftest_pathlist("collect_ignore", path=path.dirpath())
    ignore_paths = ignore_paths or []
    excludeopt = config.getoption("ignore")
    if excludeopt:
        ignore_paths.extend([py.path.local(x) for x in excludeopt])

    if py.path.local(path) in ignore_paths:
        return True

    allow_in_venv = config.getoption("collect_in_virtualenv")
    if _in_venv(path) and not allow_in_venv:
        return True

    # Skip duplicate paths.
    keepduplicates = config.getoption("keepduplicates")
    duplicate_paths = config.pluginmanager._duplicatepaths
    if not keepduplicates:
        if path in duplicate_paths:
            return True
        else:
            duplicate_paths.add(path)

    return False


def pytest_collection_modifyitems(items, config):
    deselect_prefixes = tuple(config.getoption("deselect") or [])
    if not deselect_prefixes:
        return

    remaining = []
    deselected = []
    for colitem in items:
        if colitem.nodeid.startswith(deselect_prefixes):
            deselected.append(colitem)
        else:
            remaining.append(colitem)

    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = remaining


@contextlib.contextmanager
def _patched_find_module():
    """Patch bug in pkgutil.ImpImporter.find_module

    When using pkgutil.find_loader on python<3.4 it removes symlinks
    from the path due to a call to os.path.realpath. This is not consistent
    with actually doing the import (in these versions, pkgutil and __import__
    did not share the same underlying code). This can break conftest
    discovery for pytest where symlinks are involved.

    The only supported python<3.4 by pytest is python 2.7.
    """
    if six.PY2:  # python 3.4+ uses importlib instead
        def find_module_patched(self, fullname, path=None):
            # Note: we ignore 'path' argument since it is only used via meta_path
            subname = fullname.split(".")[-1]
            if subname != fullname and self.path is None:
                return None
            if self.path is None:
                path = None
            else:
                # original: path = [os.path.realpath(self.path)]
                path = [self.path]
            try:
                file, filename, etc = pkgutil.imp.find_module(subname,
                                                              path)
            except ImportError:
                return None
            return pkgutil.ImpLoader(fullname, file, filename, etc)

        old_find_module = pkgutil.ImpImporter.find_module
        pkgutil.ImpImporter.find_module = find_module_patched
        try:
            yield
        finally:
            pkgutil.ImpImporter.find_module = old_find_module
    else:
        yield


class FSHookProxy(object):
    def __init__(self, fspath, pm, remove_mods):
        self.fspath = fspath
        self.pm = pm
        self.remove_mods = remove_mods

    def __getattr__(self, name):
        x = self.pm.subset_hook_caller(name, remove_plugins=self.remove_mods)
        self.__dict__[name] = x
        return x


class NoMatch(Exception):
    """ raised if matching cannot locate a matching names. """


class Interrupted(KeyboardInterrupt):
    """ signals an interrupted test run. """
    __module__ = 'builtins'  # for py3


class Failed(Exception):
    """ signals an stop as failed test run. """


class Session(nodes.FSCollector):
    Interrupted = Interrupted
    Failed = Failed

    def __init__(self, config):
        nodes.FSCollector.__init__(
            self, config.rootdir, parent=None,
            config=config, session=self, nodeid="")
        self.testsfailed = 0
        self.testscollected = 0
        self.shouldstop = False
        self.shouldfail = False
        self.trace = config.trace.root.get("collection")
        self._norecursepatterns = config.getini("norecursedirs")
        self.startdir = py.path.local()

        self.config.pluginmanager.register(self, name="session")

    @hookimpl(tryfirst=True)
    def pytest_collectstart(self):
        if self.shouldfail:
            raise self.Failed(self.shouldfail)
        if self.shouldstop:
            raise self.Interrupted(self.shouldstop)

    @hookimpl(tryfirst=True)
    def pytest_runtest_logreport(self, report):
        if report.failed and not hasattr(report, 'wasxfail'):
            self.testsfailed += 1
            maxfail = self.config.getvalue("maxfail")
            if maxfail and self.testsfailed >= maxfail:
                self.shouldfail = "stopping after %d failures" % (
                    self.testsfailed)
    pytest_collectreport = pytest_runtest_logreport

    def isinitpath(self, path):
        return path in self._initialpaths

    def gethookproxy(self, fspath):
        # check if we have the common case of running
        # hooks with all conftest.py files
        pm = self.config.pluginmanager
        my_conftestmodules = pm._getconftestmodules(fspath)
        remove_mods = pm._conftest_plugins.difference(my_conftestmodules)
        if remove_mods:
            # one or more conftests are not in use at this fspath
            proxy = FSHookProxy(fspath, pm, remove_mods)
        else:
            # all plugis are active for this fspath
            proxy = self.config.hook
        return proxy

    def perform_collect(self, args=None, genitems=True):
        hook = self.config.hook
        try:
            items = self._perform_collect(args, genitems)
            self.config.pluginmanager.check_pending()
            hook.pytest_collection_modifyitems(session=self,
                                               config=self.config, items=items)
        finally:
            hook.pytest_collection_finish(session=self)
        self.testscollected = len(items)
        return items

    def _perform_collect(self, args, genitems):
        if args is None:
            args = self.config.args
        self.trace("perform_collect", self, args)
        self.trace.root.indent += 1
        self._notfound = []
        self._initialpaths = set()
        self._initialparts = []
        self.items = items = []
        for arg in args:
            parts = self._parsearg(arg)
            self._initialparts.append(parts)
            self._initialpaths.add(parts[0])
        rep = collect_one_node(self)
        self.ihook.pytest_collectreport(report=rep)
        self.trace.root.indent -= 1
        if self._notfound:
            errors = []
            for arg, exc in self._notfound:
                line = "(no name %r in any of %r)" % (arg, exc.args[0])
                errors.append("not found: %s\n%s" % (arg, line))
                # XXX: test this
            raise UsageError(*errors)
        if not genitems:
            return rep.result
        else:
            if rep.passed:
                for node in rep.result:
                    self.items.extend(self.genitems(node))
            return items

    def collect(self):
        for parts in self._initialparts:
            arg = "::".join(map(str, parts))
            self.trace("processing argument", arg)
            self.trace.root.indent += 1
            try:
                for x in self._collect(arg):
                    yield x
            except NoMatch:
                # we are inside a make_report hook so
                # we cannot directly pass through the exception
                self._notfound.append((arg, sys.exc_info()[1]))

            self.trace.root.indent -= 1

    def _collect(self, arg):
        names = self._parsearg(arg)
        path = names.pop(0)
        if path.check(dir=1):
            assert not names, "invalid arg %r" % (arg,)
            for path in path.visit(fil=lambda x: x.check(file=1),
                                   rec=self._recurse, bf=True, sort=True):
                for x in self._collectfile(path):
                    yield x
        else:
            assert path.check(file=1)
            for x in self.matchnodes(self._collectfile(path), names):
                yield x

    def _collectfile(self, path):
        ihook = self.gethookproxy(path)
        if not self.isinitpath(path):
            if ihook.pytest_ignore_collect(path=path, config=self.config):
                return ()
        return ihook.pytest_collect_file(path=path, parent=self)

    def _recurse(self, path):
        ihook = self.gethookproxy(path.dirpath())
        if ihook.pytest_ignore_collect(path=path, config=self.config):
            return
        for pat in self._norecursepatterns:
            if path.check(fnmatch=pat):
                return False
        ihook = self.gethookproxy(path)
        ihook.pytest_collect_directory(path=path, parent=self)
        return True

    def _tryconvertpyarg(self, x):
        """Convert a dotted module name to path.

        """

        try:
            with _patched_find_module():
                loader = pkgutil.find_loader(x)
        except ImportError:
            return x
        if loader is None:
            return x
        # This method is sometimes invoked when AssertionRewritingHook, which
        # does not define a get_filename method, is already in place:
        try:
            with _patched_find_module():
                path = loader.get_filename(x)
        except AttributeError:
            # Retrieve path from AssertionRewritingHook:
            path = loader.modules[x][0].co_filename
        if loader.is_package(x):
            path = os.path.dirname(path)
        return path

    def _parsearg(self, arg):
        """ return (fspath, names) tuple after checking the file exists. """
        parts = str(arg).split("::")
        if self.config.option.pyargs:
            parts[0] = self._tryconvertpyarg(parts[0])
        relpath = parts[0].replace("/", os.sep)
        path = self.config.invocation_dir.join(relpath, abs=True)
        if not path.check():
            if self.config.option.pyargs:
                raise UsageError(
                    "file or package not found: " + arg +
                    " (missing __init__.py?)")
            else:
                raise UsageError("file not found: " + arg)
        parts[0] = path
        return parts

    def matchnodes(self, matching, names):
        self.trace("matchnodes", matching, names)
        self.trace.root.indent += 1
        nodes = self._matchnodes(matching, names)
        num = len(nodes)
        self.trace("matchnodes finished -> ", num, "nodes")
        self.trace.root.indent -= 1
        if num == 0:
            raise NoMatch(matching, names[:1])
        return nodes

    def _matchnodes(self, matching, names):
        if not matching or not names:
            return matching
        name = names[0]
        assert name
        nextnames = names[1:]
        resultnodes = []
        for node in matching:
            if isinstance(node, nodes.Item):
                if not names:
                    resultnodes.append(node)
                continue
            assert isinstance(node, nodes.Collector)
            rep = collect_one_node(node)
            if rep.passed:
                has_matched = False
                for x in rep.result:
                    # TODO: remove parametrized workaround once collection structure contains parametrization
                    if x.name == name or x.name.split("[")[0] == name:
                        resultnodes.extend(self.matchnodes([x], nextnames))
                        has_matched = True
                # XXX accept IDs that don't have "()" for class instances
                if not has_matched and len(rep.result) == 1 and x.name == "()":
                    nextnames.insert(0, name)
                    resultnodes.extend(self.matchnodes([x], nextnames))
            else:
                # report collection failures here to avoid failing to run some test
                # specified in the command line because the module could not be
                # imported (#134)
                node.ihook.pytest_collectreport(report=rep)
        return resultnodes

    def genitems(self, node):
        self.trace("genitems", node)
        if isinstance(node, nodes.Item):
            node.ihook.pytest_itemcollected(item=node)
            yield node
        else:
            assert isinstance(node, nodes.Collector)
            rep = collect_one_node(node)
            if rep.passed:
                for subnode in rep.result:
                    for x in self.genitems(subnode):
                        yield x
            node.ihook.pytest_collectreport(report=rep)
