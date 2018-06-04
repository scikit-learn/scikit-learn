"""(disabled by default) support for testing pytest and pytest plugins."""
from __future__ import absolute_import, division, print_function

import codecs
import gc
import os
import platform
import re
import subprocess
import six
import sys
import time
import traceback
from fnmatch import fnmatch

from weakref import WeakKeyDictionary

from _pytest.capture import MultiCapture, SysCapture
from _pytest._code import Source
import py
import pytest
from _pytest.main import Session, EXIT_OK
from _pytest.assertion.rewrite import AssertionRewritingHook


PYTEST_FULLPATH = os.path.abspath(pytest.__file__.rstrip("oc")).replace("$py.class", ".py")


IGNORE_PAM = [  # filenames added when obtaining details about the current user
     u'/var/lib/sss/mc/passwd'
]


def pytest_addoption(parser):
    parser.addoption('--lsof',
                     action="store_true", dest="lsof", default=False,
                     help=("run FD checks if lsof is available"))

    parser.addoption('--runpytest', default="inprocess", dest="runpytest",
                     choices=("inprocess", "subprocess"),
                     help=("run pytest sub runs in tests using an 'inprocess' "
                           "or 'subprocess' (python -m main) method"))


def pytest_configure(config):
    if config.getvalue("lsof"):
        checker = LsofFdLeakChecker()
        if checker.matching_platform():
            config.pluginmanager.register(checker)


class LsofFdLeakChecker(object):
    def get_open_files(self):
        out = self._exec_lsof()
        open_files = self._parse_lsof_output(out)
        return open_files

    def _exec_lsof(self):
        pid = os.getpid()
        return py.process.cmdexec("lsof -Ffn0 -p %d" % pid)

    def _parse_lsof_output(self, out):
        def isopen(line):
            return line.startswith('f') and ("deleted" not in line and
                                             'mem' not in line and "txt" not in line and 'cwd' not in line)

        open_files = []

        for line in out.split("\n"):
            if isopen(line):
                fields = line.split('\0')
                fd = fields[0][1:]
                filename = fields[1][1:]
                if filename in IGNORE_PAM:
                    continue
                if filename.startswith('/'):
                    open_files.append((fd, filename))

        return open_files

    def matching_platform(self):
        try:
            py.process.cmdexec("lsof -v")
        except (py.process.cmdexec.Error, UnicodeDecodeError):
            # cmdexec may raise UnicodeDecodeError on Windows systems with
            # locale other than English:
            # https://bitbucket.org/pytest-dev/py/issues/66
            return False
        else:
            return True

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_protocol(self, item):
        lines1 = self.get_open_files()
        yield
        if hasattr(sys, "pypy_version_info"):
            gc.collect()
        lines2 = self.get_open_files()

        new_fds = set([t[0] for t in lines2]) - set([t[0] for t in lines1])
        leaked_files = [t for t in lines2 if t[0] in new_fds]
        if leaked_files:
            error = []
            error.append("***** %s FD leakage detected" % len(leaked_files))
            error.extend([str(f) for f in leaked_files])
            error.append("*** Before:")
            error.extend([str(f) for f in lines1])
            error.append("*** After:")
            error.extend([str(f) for f in lines2])
            error.append(error[0])
            error.append("*** function %s:%s: %s " % item.location)
            error.append("See issue #2366")
            item.warn('', "\n".join(error))


# XXX copied from execnet's conftest.py - needs to be merged
winpymap = {
    'python2.7': r'C:\Python27\python.exe',
    'python3.4': r'C:\Python34\python.exe',
    'python3.5': r'C:\Python35\python.exe',
    'python3.6': r'C:\Python36\python.exe',
}


def getexecutable(name, cache={}):
    try:
        return cache[name]
    except KeyError:
        executable = py.path.local.sysfind(name)
        if executable:
            import subprocess
            popen = subprocess.Popen([str(executable), "--version"],
                                     universal_newlines=True, stderr=subprocess.PIPE)
            out, err = popen.communicate()
            if name == "jython":
                if not err or "2.5" not in err:
                    executable = None
                if "2.5.2" in err:
                    executable = None  # http://bugs.jython.org/issue1790
            elif popen.returncode != 0:
                # handle pyenv's 127
                executable = None
        cache[name] = executable
        return executable


@pytest.fixture(params=['python2.7', 'python3.4', 'pypy', 'pypy3'])
def anypython(request):
    name = request.param
    executable = getexecutable(name)
    if executable is None:
        if sys.platform == "win32":
            executable = winpymap.get(name, None)
            if executable:
                executable = py.path.local(executable)
                if executable.check():
                    return executable
        pytest.skip("no suitable %s found" % (name,))
    return executable

# used at least by pytest-xdist plugin


@pytest.fixture
def _pytest(request):
    """Return a helper which offers a gethookrecorder(hook) method which
    returns a HookRecorder instance which helps to make assertions about called
    hooks.

    """
    return PytestArg(request)


class PytestArg(object):
    def __init__(self, request):
        self.request = request

    def gethookrecorder(self, hook):
        hookrecorder = HookRecorder(hook._pm)
        self.request.addfinalizer(hookrecorder.finish_recording)
        return hookrecorder


def get_public_names(values):
    """Only return names from iterator values without a leading underscore."""
    return [x for x in values if x[0] != "_"]


class ParsedCall(object):
    def __init__(self, name, kwargs):
        self.__dict__.update(kwargs)
        self._name = name

    def __repr__(self):
        d = self.__dict__.copy()
        del d['_name']
        return "<ParsedCall %r(**%r)>" % (self._name, d)


class HookRecorder(object):
    """Record all hooks called in a plugin manager.

    This wraps all the hook calls in the plugin manager, recording each call
    before propagating the normal calls.

    """

    def __init__(self, pluginmanager):
        self._pluginmanager = pluginmanager
        self.calls = []

        def before(hook_name, hook_impls, kwargs):
            self.calls.append(ParsedCall(hook_name, kwargs))

        def after(outcome, hook_name, hook_impls, kwargs):
            pass

        self._undo_wrapping = pluginmanager.add_hookcall_monitoring(before, after)

    def finish_recording(self):
        self._undo_wrapping()

    def getcalls(self, names):
        if isinstance(names, str):
            names = names.split()
        return [call for call in self.calls if call._name in names]

    def assert_contains(self, entries):
        __tracebackhide__ = True
        i = 0
        entries = list(entries)
        backlocals = sys._getframe(1).f_locals
        while entries:
            name, check = entries.pop(0)
            for ind, call in enumerate(self.calls[i:]):
                if call._name == name:
                    print("NAMEMATCH", name, call)
                    if eval(check, backlocals, call.__dict__):
                        print("CHECKERMATCH", repr(check), "->", call)
                    else:
                        print("NOCHECKERMATCH", repr(check), "-", call)
                        continue
                    i += ind + 1
                    break
                print("NONAMEMATCH", name, "with", call)
            else:
                pytest.fail("could not find %r check %r" % (name, check))

    def popcall(self, name):
        __tracebackhide__ = True
        for i, call in enumerate(self.calls):
            if call._name == name:
                del self.calls[i]
                return call
        lines = ["could not find call %r, in:" % (name,)]
        lines.extend(["  %s" % str(x) for x in self.calls])
        pytest.fail("\n".join(lines))

    def getcall(self, name):
        values = self.getcalls(name)
        assert len(values) == 1, (name, values)
        return values[0]

    # functionality for test reports

    def getreports(self,
                   names="pytest_runtest_logreport pytest_collectreport"):
        return [x.report for x in self.getcalls(names)]

    def matchreport(self, inamepart="",
                    names="pytest_runtest_logreport pytest_collectreport", when=None):
        """return a testreport whose dotted import path matches"""
        values = []
        for rep in self.getreports(names=names):
            try:
                if not when and rep.when != "call" and rep.passed:
                    # setup/teardown passing reports - let's ignore those
                    continue
            except AttributeError:
                pass
            if when and getattr(rep, 'when', None) != when:
                continue
            if not inamepart or inamepart in rep.nodeid.split("::"):
                values.append(rep)
        if not values:
            raise ValueError("could not find test report matching %r: "
                             "no test reports at all!" % (inamepart,))
        if len(values) > 1:
            raise ValueError(
                "found 2 or more testreports matching %r: %s" % (inamepart, values))
        return values[0]

    def getfailures(self,
                    names='pytest_runtest_logreport pytest_collectreport'):
        return [rep for rep in self.getreports(names) if rep.failed]

    def getfailedcollections(self):
        return self.getfailures('pytest_collectreport')

    def listoutcomes(self):
        passed = []
        skipped = []
        failed = []
        for rep in self.getreports(
                "pytest_collectreport pytest_runtest_logreport"):
            if rep.passed:
                if getattr(rep, "when", None) == "call":
                    passed.append(rep)
            elif rep.skipped:
                skipped.append(rep)
            elif rep.failed:
                failed.append(rep)
        return passed, skipped, failed

    def countoutcomes(self):
        return [len(x) for x in self.listoutcomes()]

    def assertoutcome(self, passed=0, skipped=0, failed=0):
        realpassed, realskipped, realfailed = self.listoutcomes()
        assert passed == len(realpassed)
        assert skipped == len(realskipped)
        assert failed == len(realfailed)

    def clear(self):
        self.calls[:] = []


@pytest.fixture
def linecomp(request):
    return LineComp()


@pytest.fixture(name='LineMatcher')
def LineMatcher_fixture(request):
    return LineMatcher


@pytest.fixture
def testdir(request, tmpdir_factory):
    return Testdir(request, tmpdir_factory)


rex_outcome = re.compile(r"(\d+) ([\w-]+)")


class RunResult(object):
    """The result of running a command.

    Attributes:

    :ret: the return value
    :outlines: list of lines captured from stdout
    :errlines: list of lines captures from stderr
    :stdout: :py:class:`LineMatcher` of stdout, use ``stdout.str()`` to
       reconstruct stdout or the commonly used ``stdout.fnmatch_lines()``
       method
    :stderr: :py:class:`LineMatcher` of stderr
    :duration: duration in seconds

    """

    def __init__(self, ret, outlines, errlines, duration):
        self.ret = ret
        self.outlines = outlines
        self.errlines = errlines
        self.stdout = LineMatcher(outlines)
        self.stderr = LineMatcher(errlines)
        self.duration = duration

    def parseoutcomes(self):
        """Return a dictionary of outcomestring->num from parsing the terminal
        output that the test process produced.

        """
        for line in reversed(self.outlines):
            if 'seconds' in line:
                outcomes = rex_outcome.findall(line)
                if outcomes:
                    d = {}
                    for num, cat in outcomes:
                        d[cat] = int(num)
                    return d
        raise ValueError("Pytest terminal report not found")

    def assert_outcomes(self, passed=0, skipped=0, failed=0, error=0):
        """Assert that the specified outcomes appear with the respective
        numbers (0 means it didn't occur) in the text output from a test run.

        """
        d = self.parseoutcomes()
        obtained = {
            'passed': d.get('passed', 0),
            'skipped': d.get('skipped', 0),
            'failed': d.get('failed', 0),
            'error': d.get('error', 0),
        }
        assert obtained == dict(passed=passed, skipped=skipped, failed=failed, error=error)


class CwdSnapshot(object):
    def __init__(self):
        self.__saved = os.getcwd()

    def restore(self):
        os.chdir(self.__saved)


class SysModulesSnapshot(object):
    def __init__(self, preserve=None):
        self.__preserve = preserve
        self.__saved = dict(sys.modules)

    def restore(self):
        if self.__preserve:
            self.__saved.update(
                (k, m) for k, m in sys.modules.items() if self.__preserve(k))
        sys.modules.clear()
        sys.modules.update(self.__saved)


class SysPathsSnapshot(object):
    def __init__(self):
        self.__saved = list(sys.path), list(sys.meta_path)

    def restore(self):
        sys.path[:], sys.meta_path[:] = self.__saved


class Testdir(object):
    """Temporary test directory with tools to test/run pytest itself.

    This is based on the ``tmpdir`` fixture but provides a number of methods
    which aid with testing pytest itself.  Unless :py:meth:`chdir` is used all
    methods will use :py:attr:`tmpdir` as their current working directory.

    Attributes:

    :tmpdir: The :py:class:`py.path.local` instance of the temporary directory.

    :plugins: A list of plugins to use with :py:meth:`parseconfig` and
       :py:meth:`runpytest`.  Initially this is an empty list but plugins can
       be added to the list.  The type of items to add to the list depends on
       the method using them so refer to them for details.

    """

    def __init__(self, request, tmpdir_factory):
        self.request = request
        self._mod_collections = WeakKeyDictionary()
        name = request.function.__name__
        self.tmpdir = tmpdir_factory.mktemp(name, numbered=True)
        self.plugins = []
        self._cwd_snapshot = CwdSnapshot()
        self._sys_path_snapshot = SysPathsSnapshot()
        self._sys_modules_snapshot = self.__take_sys_modules_snapshot()
        self.chdir()
        self.request.addfinalizer(self.finalize)
        method = self.request.config.getoption("--runpytest")
        if method == "inprocess":
            self._runpytest_method = self.runpytest_inprocess
        elif method == "subprocess":
            self._runpytest_method = self.runpytest_subprocess

    def __repr__(self):
        return "<Testdir %r>" % (self.tmpdir,)

    def finalize(self):
        """Clean up global state artifacts.

        Some methods modify the global interpreter state and this tries to
        clean this up.  It does not remove the temporary directory however so
        it can be looked at after the test run has finished.

        """
        self._sys_modules_snapshot.restore()
        self._sys_path_snapshot.restore()
        self._cwd_snapshot.restore()

    def __take_sys_modules_snapshot(self):
        # some zope modules used by twisted-related tests keep internal state
        # and can't be deleted; we had some trouble in the past with
        # `zope.interface` for example
        def preserve_module(name):
            return name.startswith("zope")
        return SysModulesSnapshot(preserve=preserve_module)

    def make_hook_recorder(self, pluginmanager):
        """Create a new :py:class:`HookRecorder` for a PluginManager."""
        assert not hasattr(pluginmanager, "reprec")
        pluginmanager.reprec = reprec = HookRecorder(pluginmanager)
        self.request.addfinalizer(reprec.finish_recording)
        return reprec

    def chdir(self):
        """Cd into the temporary directory.

        This is done automatically upon instantiation.

        """
        self.tmpdir.chdir()

    def _makefile(self, ext, args, kwargs, encoding='utf-8'):
        items = list(kwargs.items())

        def to_text(s):
            return s.decode(encoding) if isinstance(s, bytes) else six.text_type(s)

        if args:
            source = u"\n".join(to_text(x) for x in args)
            basename = self.request.function.__name__
            items.insert(0, (basename, source))

        ret = None
        for basename, value in items:
            p = self.tmpdir.join(basename).new(ext=ext)
            p.dirpath().ensure_dir()
            source = Source(value)
            source = u"\n".join(to_text(line) for line in source.lines)
            p.write(source.strip().encode(encoding), "wb")
            if ret is None:
                ret = p
        return ret

    def makefile(self, ext, *args, **kwargs):
        """Create a new file in the testdir.

        ext: The extension the file should use, including the dot, e.g. `.py`.

        args: All args will be treated as strings and joined using newlines.
           The result will be written as contents to the file.  The name of the
           file will be based on the test function requesting this fixture.
           E.g. "testdir.makefile('.txt', 'line1', 'line2')"

        kwargs: Each keyword is the name of a file, while the value of it will
           be written as contents of the file.
           E.g. "testdir.makefile('.ini', pytest='[pytest]\naddopts=-rs\n')"

        """
        return self._makefile(ext, args, kwargs)

    def makeconftest(self, source):
        """Write a contest.py file with 'source' as contents."""
        return self.makepyfile(conftest=source)

    def makeini(self, source):
        """Write a tox.ini file with 'source' as contents."""
        return self.makefile('.ini', tox=source)

    def getinicfg(self, source):
        """Return the pytest section from the tox.ini config file."""
        p = self.makeini(source)
        return py.iniconfig.IniConfig(p)['pytest']

    def makepyfile(self, *args, **kwargs):
        """Shortcut for .makefile() with a .py extension."""
        return self._makefile('.py', args, kwargs)

    def maketxtfile(self, *args, **kwargs):
        """Shortcut for .makefile() with a .txt extension."""
        return self._makefile('.txt', args, kwargs)

    def syspathinsert(self, path=None):
        """Prepend a directory to sys.path, defaults to :py:attr:`tmpdir`.

        This is undone automatically when this object dies at the end of each
        test.

        """
        if path is None:
            path = self.tmpdir
        sys.path.insert(0, str(path))
        # a call to syspathinsert() usually means that the caller wants to
        # import some dynamically created files, thus with python3 we
        # invalidate its import caches
        self._possibly_invalidate_import_caches()

    def _possibly_invalidate_import_caches(self):
        # invalidate caches if we can (py33 and above)
        try:
            import importlib
        except ImportError:
            pass
        else:
            if hasattr(importlib, "invalidate_caches"):
                importlib.invalidate_caches()

    def mkdir(self, name):
        """Create a new (sub)directory."""
        return self.tmpdir.mkdir(name)

    def mkpydir(self, name):
        """Create a new python package.

        This creates a (sub)directory with an empty ``__init__.py`` file so it
        gets recognised as a python package.

        """
        p = self.mkdir(name)
        p.ensure("__init__.py")
        return p

    Session = Session

    def getnode(self, config, arg):
        """Return the collection node of a file.

        :param config: :py:class:`_pytest.config.Config` instance, see
           :py:meth:`parseconfig` and :py:meth:`parseconfigure` to create the
           configuration

        :param arg: a :py:class:`py.path.local` instance of the file

        """
        session = Session(config)
        assert '::' not in str(arg)
        p = py.path.local(arg)
        config.hook.pytest_sessionstart(session=session)
        res = session.perform_collect([str(p)], genitems=False)[0]
        config.hook.pytest_sessionfinish(session=session, exitstatus=EXIT_OK)
        return res

    def getpathnode(self, path):
        """Return the collection node of a file.

        This is like :py:meth:`getnode` but uses :py:meth:`parseconfigure` to
        create the (configured) pytest Config instance.

        :param path: a :py:class:`py.path.local` instance of the file

        """
        config = self.parseconfigure(path)
        session = Session(config)
        x = session.fspath.bestrelpath(path)
        config.hook.pytest_sessionstart(session=session)
        res = session.perform_collect([x], genitems=False)[0]
        config.hook.pytest_sessionfinish(session=session, exitstatus=EXIT_OK)
        return res

    def genitems(self, colitems):
        """Generate all test items from a collection node.

        This recurses into the collection node and returns a list of all the
        test items contained within.

        """
        session = colitems[0].session
        result = []
        for colitem in colitems:
            result.extend(session.genitems(colitem))
        return result

    def runitem(self, source):
        """Run the "test_func" Item.

        The calling test instance (class containing the test method) must
        provide a ``.getrunner()`` method which should return a runner which
        can run the test protocol for a single item, e.g.
        :py:func:`_pytest.runner.runtestprotocol`.

        """
        # used from runner functional tests
        item = self.getitem(source)
        # the test class where we are called from wants to provide the runner
        testclassinstance = self.request.instance
        runner = testclassinstance.getrunner()
        return runner(item)

    def inline_runsource(self, source, *cmdlineargs):
        """Run a test module in process using ``pytest.main()``.

        This run writes "source" into a temporary file and runs
        ``pytest.main()`` on it, returning a :py:class:`HookRecorder` instance
        for the result.

        :param source: the source code of the test module

        :param cmdlineargs: any extra command line arguments to use

        :return: :py:class:`HookRecorder` instance of the result

        """
        p = self.makepyfile(source)
        values = list(cmdlineargs) + [p]
        return self.inline_run(*values)

    def inline_genitems(self, *args):
        """Run ``pytest.main(['--collectonly'])`` in-process.

        Runs the :py:func:`pytest.main` function to run all of pytest inside
        the test process itself like :py:meth:`inline_run`, but returns a
        tuple of the collected items and a :py:class:`HookRecorder` instance.

        """
        rec = self.inline_run("--collect-only", *args)
        items = [x.item for x in rec.getcalls("pytest_itemcollected")]
        return items, rec

    def inline_run(self, *args, **kwargs):
        """Run ``pytest.main()`` in-process, returning a HookRecorder.

        Runs the :py:func:`pytest.main` function to run all of pytest inside
        the test process itself.  This means it can return a
        :py:class:`HookRecorder` instance which gives more detailed results
        from that run than can be done by matching stdout/stderr from
        :py:meth:`runpytest`.

        :param args: command line arguments to pass to :py:func:`pytest.main`

        :param plugin: (keyword-only) extra plugin instances the
           ``pytest.main()`` instance should use

        :return: a :py:class:`HookRecorder` instance

        """
        finalizers = []
        try:
            # When running py.test inline any plugins active in the main test
            # process are already imported.  So this disables the warning which
            # will trigger to say they can no longer be rewritten, which is
            # fine as they have already been rewritten.
            orig_warn = AssertionRewritingHook._warn_already_imported

            def revert_warn_already_imported():
                AssertionRewritingHook._warn_already_imported = orig_warn
            finalizers.append(revert_warn_already_imported)
            AssertionRewritingHook._warn_already_imported = lambda *a: None

            # Any sys.module or sys.path changes done while running py.test
            # inline should be reverted after the test run completes to avoid
            # clashing with later inline tests run within the same pytest test,
            # e.g. just because they use matching test module names.
            finalizers.append(self.__take_sys_modules_snapshot().restore)
            finalizers.append(SysPathsSnapshot().restore)

            # Important note:
            # - our tests should not leave any other references/registrations
            #   laying around other than possibly loaded test modules
            #   referenced from sys.modules, as nothing will clean those up
            #   automatically

            rec = []

            class Collect(object):
                def pytest_configure(x, config):
                    rec.append(self.make_hook_recorder(config.pluginmanager))

            plugins = kwargs.get("plugins") or []
            plugins.append(Collect())
            ret = pytest.main(list(args), plugins=plugins)
            if len(rec) == 1:
                reprec = rec.pop()
            else:
                class reprec(object):
                    pass
            reprec.ret = ret

            # typically we reraise keyboard interrupts from the child run
            # because it's our user requesting interruption of the testing
            if ret == 2 and not kwargs.get("no_reraise_ctrlc"):
                calls = reprec.getcalls("pytest_keyboard_interrupt")
                if calls and calls[-1].excinfo.type == KeyboardInterrupt:
                    raise KeyboardInterrupt()
            return reprec
        finally:
            for finalizer in finalizers:
                finalizer()

    def runpytest_inprocess(self, *args, **kwargs):
        """Return result of running pytest in-process, providing a similar
        interface to what self.runpytest() provides.

        """
        if kwargs.get("syspathinsert"):
            self.syspathinsert()
        now = time.time()
        capture = MultiCapture(Capture=SysCapture)
        capture.start_capturing()
        try:
            try:
                reprec = self.inline_run(*args, **kwargs)
            except SystemExit as e:

                class reprec(object):
                    ret = e.args[0]

            except Exception:
                traceback.print_exc()

                class reprec(object):
                    ret = 3
        finally:
            out, err = capture.readouterr()
            capture.stop_capturing()
            sys.stdout.write(out)
            sys.stderr.write(err)

        res = RunResult(reprec.ret,
                        out.split("\n"), err.split("\n"),
                        time.time() - now)
        res.reprec = reprec
        return res

    def runpytest(self, *args, **kwargs):
        """Run pytest inline or in a subprocess, depending on the command line
        option "--runpytest" and return a :py:class:`RunResult`.

        """
        args = self._ensure_basetemp(args)
        return self._runpytest_method(*args, **kwargs)

    def _ensure_basetemp(self, args):
        args = [str(x) for x in args]
        for x in args:
            if str(x).startswith('--basetemp'):
                # print("basedtemp exists: %s" %(args,))
                break
        else:
            args.append("--basetemp=%s" % self.tmpdir.dirpath('basetemp'))
            # print("added basetemp: %s" %(args,))
        return args

    def parseconfig(self, *args):
        """Return a new pytest Config instance from given commandline args.

        This invokes the pytest bootstrapping code in _pytest.config to create
        a new :py:class:`_pytest.core.PluginManager` and call the
        pytest_cmdline_parse hook to create a new
        :py:class:`_pytest.config.Config` instance.

        If :py:attr:`plugins` has been populated they should be plugin modules
        to be registered with the PluginManager.

        """
        args = self._ensure_basetemp(args)

        import _pytest.config
        config = _pytest.config._prepareconfig(args, self.plugins)
        # we don't know what the test will do with this half-setup config
        # object and thus we make sure it gets unconfigured properly in any
        # case (otherwise capturing could still be active, for example)
        self.request.addfinalizer(config._ensure_unconfigure)
        return config

    def parseconfigure(self, *args):
        """Return a new pytest configured Config instance.

        This returns a new :py:class:`_pytest.config.Config` instance like
        :py:meth:`parseconfig`, but also calls the pytest_configure hook.

        """
        config = self.parseconfig(*args)
        config._do_configure()
        self.request.addfinalizer(config._ensure_unconfigure)
        return config

    def getitem(self, source, funcname="test_func"):
        """Return the test item for a test function.

        This writes the source to a python file and runs pytest's collection on
        the resulting module, returning the test item for the requested
        function name.

        :param source: the module source

        :param funcname: the name of the test function for which to return a
            test item

        """
        items = self.getitems(source)
        for item in items:
            if item.name == funcname:
                return item
        assert 0, "%r item not found in module:\n%s\nitems: %s" % (
                  funcname, source, items)

    def getitems(self, source):
        """Return all test items collected from the module.

        This writes the source to a python file and runs pytest's collection on
        the resulting module, returning all test items contained within.

        """
        modcol = self.getmodulecol(source)
        return self.genitems([modcol])

    def getmodulecol(self, source, configargs=(), withinit=False):
        """Return the module collection node for ``source``.

        This writes ``source`` to a file using :py:meth:`makepyfile` and then
        runs the pytest collection on it, returning the collection node for the
        test module.

        :param source: the source code of the module to collect

        :param configargs: any extra arguments to pass to
            :py:meth:`parseconfigure`

        :param withinit: whether to also write an ``__init__.py`` file to the
            same directory to ensure it is a package

        """
        kw = {self.request.function.__name__: Source(source).strip()}
        path = self.makepyfile(**kw)
        if withinit:
            self.makepyfile(__init__="#")
        self.config = config = self.parseconfigure(path, *configargs)
        node = self.getnode(config, path)

        return node

    def collect_by_name(self, modcol, name):
        """Return the collection node for name from the module collection.

        This will search a module collection node for a collection node
        matching the given name.

        :param modcol: a module collection node; see :py:meth:`getmodulecol`

        :param name: the name of the node to return

        """
        if modcol not in self._mod_collections:
            self._mod_collections[modcol] = list(modcol.collect())
        for colitem in self._mod_collections[modcol]:
            if colitem.name == name:
                return colitem

    def popen(self, cmdargs, stdout, stderr, **kw):
        """Invoke subprocess.Popen.

        This calls subprocess.Popen making sure the current working directory
        is in the PYTHONPATH.

        You probably want to use :py:meth:`run` instead.

        """
        env = os.environ.copy()
        env['PYTHONPATH'] = os.pathsep.join(filter(None, [
            str(os.getcwd()), env.get('PYTHONPATH', '')]))
        kw['env'] = env

        popen = subprocess.Popen(cmdargs, stdin=subprocess.PIPE, stdout=stdout, stderr=stderr, **kw)
        popen.stdin.close()

        return popen

    def run(self, *cmdargs):
        """Run a command with arguments.

        Run a process using subprocess.Popen saving the stdout and stderr.

        Returns a :py:class:`RunResult`.

        """
        return self._run(*cmdargs)

    def _run(self, *cmdargs):
        cmdargs = [str(x) for x in cmdargs]
        p1 = self.tmpdir.join("stdout")
        p2 = self.tmpdir.join("stderr")
        print("running:", ' '.join(cmdargs))
        print("     in:", str(py.path.local()))
        f1 = codecs.open(str(p1), "w", encoding="utf8")
        f2 = codecs.open(str(p2), "w", encoding="utf8")
        try:
            now = time.time()
            popen = self.popen(cmdargs, stdout=f1, stderr=f2,
                               close_fds=(sys.platform != "win32"))
            ret = popen.wait()
        finally:
            f1.close()
            f2.close()
        f1 = codecs.open(str(p1), "r", encoding="utf8")
        f2 = codecs.open(str(p2), "r", encoding="utf8")
        try:
            out = f1.read().splitlines()
            err = f2.read().splitlines()
        finally:
            f1.close()
            f2.close()
        self._dump_lines(out, sys.stdout)
        self._dump_lines(err, sys.stderr)
        return RunResult(ret, out, err, time.time() - now)

    def _dump_lines(self, lines, fp):
        try:
            for line in lines:
                print(line, file=fp)
        except UnicodeEncodeError:
            print("couldn't print to %s because of encoding" % (fp,))

    def _getpytestargs(self):
        # we cannot use `(sys.executable, script)` because on Windows the
        # script is e.g. `pytest.exe`
        return (sys.executable, PYTEST_FULLPATH) # noqa

    def runpython(self, script):
        """Run a python script using sys.executable as interpreter.

        Returns a :py:class:`RunResult`.

        """
        return self.run(sys.executable, script)

    def runpython_c(self, command):
        """Run python -c "command", return a :py:class:`RunResult`."""
        return self.run(sys.executable, "-c", command)

    def runpytest_subprocess(self, *args, **kwargs):
        """Run pytest as a subprocess with given arguments.

        Any plugins added to the :py:attr:`plugins` list will added using the
        ``-p`` command line option.  Additionally ``--basetemp`` is used put
        any temporary files and directories in a numbered directory prefixed
        with "runpytest-" so they do not conflict with the normal numbered
        pytest location for temporary files and directories.

        Returns a :py:class:`RunResult`.

        """
        p = py.path.local.make_numbered_dir(prefix="runpytest-",
                                            keep=None, rootdir=self.tmpdir)
        args = ('--basetemp=%s' % p,) + args
        plugins = [x for x in self.plugins if isinstance(x, str)]
        if plugins:
            args = ('-p', plugins[0]) + args
        args = self._getpytestargs() + args
        return self.run(*args)

    def spawn_pytest(self, string, expect_timeout=10.0):
        """Run pytest using pexpect.

        This makes sure to use the right pytest and sets up the temporary
        directory locations.

        The pexpect child is returned.

        """
        basetemp = self.tmpdir.mkdir("temp-pexpect")
        invoke = " ".join(map(str, self._getpytestargs()))
        cmd = "%s --basetemp=%s %s" % (invoke, basetemp, string)
        return self.spawn(cmd, expect_timeout=expect_timeout)

    def spawn(self, cmd, expect_timeout=10.0):
        """Run a command using pexpect.

        The pexpect child is returned.

        """
        pexpect = pytest.importorskip("pexpect", "3.0")
        if hasattr(sys, 'pypy_version_info') and '64' in platform.machine():
            pytest.skip("pypy-64 bit not supported")
        if sys.platform.startswith("freebsd"):
            pytest.xfail("pexpect does not work reliably on freebsd")
        logfile = self.tmpdir.join("spawn.out").open("wb")
        child = pexpect.spawn(cmd, logfile=logfile)
        self.request.addfinalizer(logfile.close)
        child.timeout = expect_timeout
        return child


def getdecoded(out):
    try:
        return out.decode("utf-8")
    except UnicodeDecodeError:
        return "INTERNAL not-utf8-decodeable, truncated string:\n%s" % (
            py.io.saferepr(out),)


class LineComp(object):
    def __init__(self):
        self.stringio = py.io.TextIO()

    def assert_contains_lines(self, lines2):
        """Assert that lines2 are contained (linearly) in lines1.

        Return a list of extralines found.

        """
        __tracebackhide__ = True
        val = self.stringio.getvalue()
        self.stringio.truncate(0)
        self.stringio.seek(0)
        lines1 = val.split("\n")
        return LineMatcher(lines1).fnmatch_lines(lines2)


class LineMatcher(object):
    """Flexible matching of text.

    This is a convenience class to test large texts like the output of
    commands.

    The constructor takes a list of lines without their trailing newlines, i.e.
    ``text.splitlines()``.

    """

    def __init__(self, lines):
        self.lines = lines
        self._log_output = []

    def str(self):
        """Return the entire original text."""
        return "\n".join(self.lines)

    def _getlines(self, lines2):
        if isinstance(lines2, str):
            lines2 = Source(lines2)
        if isinstance(lines2, Source):
            lines2 = lines2.strip().lines
        return lines2

    def fnmatch_lines_random(self, lines2):
        """Check lines exist in the output using in any order.

        Lines are checked using ``fnmatch.fnmatch``. The argument is a list of
        lines which have to occur in the output, in any order.

        """
        self._match_lines_random(lines2, fnmatch)

    def re_match_lines_random(self, lines2):
        """Check lines exist in the output using ``re.match``, in any order.

        The argument is a list of lines which have to occur in the output, in
        any order.

        """
        self._match_lines_random(lines2, lambda name, pat: re.match(pat, name))

    def _match_lines_random(self, lines2, match_func):
        """Check lines exist in the output.

        The argument is a list of lines which have to occur in the output, in
        any order.  Each line can contain glob whildcards.

        """
        lines2 = self._getlines(lines2)
        for line in lines2:
            for x in self.lines:
                if line == x or match_func(x, line):
                    self._log("matched: ", repr(line))
                    break
            else:
                self._log("line %r not found in output" % line)
                raise ValueError(self._log_text)

    def get_lines_after(self, fnline):
        """Return all lines following the given line in the text.

        The given line can contain glob wildcards.

        """
        for i, line in enumerate(self.lines):
            if fnline == line or fnmatch(line, fnline):
                return self.lines[i + 1:]
        raise ValueError("line %r not found in output" % fnline)

    def _log(self, *args):
        self._log_output.append(' '.join((str(x) for x in args)))

    @property
    def _log_text(self):
        return '\n'.join(self._log_output)

    def fnmatch_lines(self, lines2):
        """Search captured text for matching lines using ``fnmatch.fnmatch``.

        The argument is a list of lines which have to match and can use glob
        wildcards.  If they do not match a pytest.fail() is called.  The
        matches and non-matches are also printed on stdout.

        """
        self._match_lines(lines2, fnmatch, 'fnmatch')

    def re_match_lines(self, lines2):
        """Search captured text for matching lines using ``re.match``.

        The argument is a list of lines which have to match using ``re.match``.
        If they do not match a pytest.fail() is called.

        The matches and non-matches are also printed on stdout.

        """
        self._match_lines(lines2, lambda name, pat: re.match(pat, name), 're.match')

    def _match_lines(self, lines2, match_func, match_nickname):
        """Underlying implementation of ``fnmatch_lines`` and ``re_match_lines``.

        :param list[str] lines2: list of string patterns to match. The actual
            format depends on ``match_func``
        :param match_func: a callable ``match_func(line, pattern)`` where line
            is the captured line from stdout/stderr and pattern is the matching
            pattern
        :param str match_nickname: the nickname for the match function that
            will be logged to stdout when a match occurs

        """
        lines2 = self._getlines(lines2)
        lines1 = self.lines[:]
        nextline = None
        extralines = []
        __tracebackhide__ = True
        for line in lines2:
            nomatchprinted = False
            while lines1:
                nextline = lines1.pop(0)
                if line == nextline:
                    self._log("exact match:", repr(line))
                    break
                elif match_func(nextline, line):
                    self._log("%s:" % match_nickname, repr(line))
                    self._log("   with:", repr(nextline))
                    break
                else:
                    if not nomatchprinted:
                        self._log("nomatch:", repr(line))
                        nomatchprinted = True
                    self._log("    and:", repr(nextline))
                extralines.append(nextline)
            else:
                self._log("remains unmatched: %r" % (line,))
                pytest.fail(self._log_text)
