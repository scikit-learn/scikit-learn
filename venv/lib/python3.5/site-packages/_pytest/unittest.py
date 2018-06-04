""" discovery and running of std-library "unittest" style tests. """
from __future__ import absolute_import, division, print_function

import sys
import traceback

# for transferring markers
import _pytest._code
from _pytest.config import hookimpl
from _pytest.outcomes import fail, skip, xfail
from _pytest.python import transfer_markers, Class, Module, Function


def pytest_pycollect_makeitem(collector, name, obj):
    # has unittest been imported and is obj a subclass of its TestCase?
    try:
        if not issubclass(obj, sys.modules["unittest"].TestCase):
            return
    except Exception:
        return
    # yes, so let's collect it
    return UnitTestCase(name, parent=collector)


class UnitTestCase(Class):
    # marker for fixturemanger.getfixtureinfo()
    # to declare that our children do not support funcargs
    nofuncargs = True

    def setup(self):
        cls = self.obj
        if getattr(cls, '__unittest_skip__', False):
            return  # skipped
        setup = getattr(cls, 'setUpClass', None)
        if setup is not None:
            setup()
        teardown = getattr(cls, 'tearDownClass', None)
        if teardown is not None:
            self.addfinalizer(teardown)
        super(UnitTestCase, self).setup()

    def collect(self):
        from unittest import TestLoader
        cls = self.obj
        if not getattr(cls, "__test__", True):
            return
        self.session._fixturemanager.parsefactories(self, unittest=True)
        loader = TestLoader()
        module = self.getparent(Module).obj
        foundsomething = False
        for name in loader.getTestCaseNames(self.obj):
            x = getattr(self.obj, name)
            if not getattr(x, '__test__', True):
                continue
            funcobj = getattr(x, 'im_func', x)
            transfer_markers(funcobj, cls, module)
            yield TestCaseFunction(name, parent=self)
            foundsomething = True

        if not foundsomething:
            runtest = getattr(self.obj, 'runTest', None)
            if runtest is not None:
                ut = sys.modules.get("twisted.trial.unittest", None)
                if ut is None or runtest != ut.TestCase.runTest:
                    yield TestCaseFunction('runTest', parent=self)


class TestCaseFunction(Function):
    _excinfo = None

    def setup(self):
        self._testcase = self.parent.obj(self.name)
        self._fix_unittest_skip_decorator()
        self._obj = getattr(self._testcase, self.name)
        if hasattr(self._testcase, 'setup_method'):
            self._testcase.setup_method(self._obj)
        if hasattr(self, "_request"):
            self._request._fillfixtures()

    def _fix_unittest_skip_decorator(self):
        """
        The @unittest.skip decorator calls functools.wraps(self._testcase)
        The call to functools.wraps() fails unless self._testcase
        has a __name__ attribute. This is usually automatically supplied
        if the test is a function or method, but we need to add manually
        here.

        See issue #1169
        """
        if sys.version_info[0] == 2:
            setattr(self._testcase, "__name__", self.name)

    def teardown(self):
        if hasattr(self._testcase, 'teardown_method'):
            self._testcase.teardown_method(self._obj)
        # Allow garbage collection on TestCase instance attributes.
        self._testcase = None
        self._obj = None

    def startTest(self, testcase):
        pass

    def _addexcinfo(self, rawexcinfo):
        # unwrap potential exception info (see twisted trial support below)
        rawexcinfo = getattr(rawexcinfo, '_rawexcinfo', rawexcinfo)
        try:
            excinfo = _pytest._code.ExceptionInfo(rawexcinfo)
        except TypeError:
            try:
                try:
                    values = traceback.format_exception(*rawexcinfo)
                    values.insert(0, "NOTE: Incompatible Exception Representation, "
                                  "displaying natively:\n\n")
                    fail("".join(values), pytrace=False)
                except (fail.Exception, KeyboardInterrupt):
                    raise
                except:  # noqa
                    fail("ERROR: Unknown Incompatible Exception "
                         "representation:\n%r" % (rawexcinfo,), pytrace=False)
            except KeyboardInterrupt:
                raise
            except fail.Exception:
                excinfo = _pytest._code.ExceptionInfo()
        self.__dict__.setdefault('_excinfo', []).append(excinfo)

    def addError(self, testcase, rawexcinfo):
        self._addexcinfo(rawexcinfo)

    def addFailure(self, testcase, rawexcinfo):
        self._addexcinfo(rawexcinfo)

    def addSkip(self, testcase, reason):
        try:
            skip(reason)
        except skip.Exception:
            self._skipped_by_mark = True
            self._addexcinfo(sys.exc_info())

    def addExpectedFailure(self, testcase, rawexcinfo, reason=""):
        try:
            xfail(str(reason))
        except xfail.Exception:
            self._addexcinfo(sys.exc_info())

    def addUnexpectedSuccess(self, testcase, reason=""):
        self._unexpectedsuccess = reason

    def addSuccess(self, testcase):
        pass

    def stopTest(self, testcase):
        pass

    def _handle_skip(self):
        # implements the skipping machinery (see #2137)
        # analog to pythons Lib/unittest/case.py:run
        testMethod = getattr(self._testcase, self._testcase._testMethodName)
        if (getattr(self._testcase.__class__, "__unittest_skip__", False) or
                getattr(testMethod, "__unittest_skip__", False)):
            # If the class or method was skipped.
            skip_why = (getattr(self._testcase.__class__, '__unittest_skip_why__', '') or
                        getattr(testMethod, '__unittest_skip_why__', ''))
            try:  # PY3, unittest2 on PY2
                self._testcase._addSkip(self, self._testcase, skip_why)
            except TypeError:  # PY2
                if sys.version_info[0] != 2:
                    raise
                self._testcase._addSkip(self, skip_why)
            return True
        return False

    def runtest(self):
        if self.config.pluginmanager.get_plugin("pdbinvoke") is None:
            self._testcase(result=self)
        else:
            # disables tearDown and cleanups for post mortem debugging (see #1890)
            if self._handle_skip():
                return
            self._testcase.debug()

    def _prunetraceback(self, excinfo):
        Function._prunetraceback(self, excinfo)
        traceback = excinfo.traceback.filter(
            lambda x: not x.frame.f_globals.get('__unittest'))
        if traceback:
            excinfo.traceback = traceback


@hookimpl(tryfirst=True)
def pytest_runtest_makereport(item, call):
    if isinstance(item, TestCaseFunction):
        if item._excinfo:
            call.excinfo = item._excinfo.pop(0)
            try:
                del call.result
            except AttributeError:
                pass

# twisted trial support


@hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item):
    if isinstance(item, TestCaseFunction) and \
       'twisted.trial.unittest' in sys.modules:
        ut = sys.modules['twisted.python.failure']
        Failure__init__ = ut.Failure.__init__
        check_testcase_implements_trial_reporter()

        def excstore(self, exc_value=None, exc_type=None, exc_tb=None,
                     captureVars=None):
            if exc_value is None:
                self._rawexcinfo = sys.exc_info()
            else:
                if exc_type is None:
                    exc_type = type(exc_value)
                self._rawexcinfo = (exc_type, exc_value, exc_tb)
            try:
                Failure__init__(self, exc_value, exc_type, exc_tb,
                                captureVars=captureVars)
            except TypeError:
                Failure__init__(self, exc_value, exc_type, exc_tb)

        ut.Failure.__init__ = excstore
        yield
        ut.Failure.__init__ = Failure__init__
    else:
        yield


def check_testcase_implements_trial_reporter(done=[]):
    if done:
        return
    from zope.interface import classImplements
    from twisted.trial.itrial import IReporter
    classImplements(TestCaseFunction, IReporter)
    done.append(1)
