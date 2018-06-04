""" basic collect and runtest protocol implementations """
from __future__ import absolute_import, division, print_function

import bdb
import os
import sys
from time import time

import py
from _pytest._code.code import TerminalRepr, ExceptionInfo
from _pytest.outcomes import skip, Skipped, TEST_OUTCOME

#
# pytest plugin hooks


def pytest_addoption(parser):
    group = parser.getgroup("terminal reporting", "reporting", after="general")
    group.addoption('--durations',
                    action="store", type=int, default=None, metavar="N",
                    help="show N slowest setup/test durations (N=0 for all)."),


def pytest_terminal_summary(terminalreporter):
    durations = terminalreporter.config.option.durations
    if durations is None:
        return
    tr = terminalreporter
    dlist = []
    for replist in tr.stats.values():
        for rep in replist:
            if hasattr(rep, 'duration'):
                dlist.append(rep)
    if not dlist:
        return
    dlist.sort(key=lambda x: x.duration)
    dlist.reverse()
    if not durations:
        tr.write_sep("=", "slowest test durations")
    else:
        tr.write_sep("=", "slowest %s test durations" % durations)
        dlist = dlist[:durations]

    for rep in dlist:
        nodeid = rep.nodeid.replace("::()::", "::")
        tr.write_line("%02.2fs %-8s %s" %
                      (rep.duration, rep.when, nodeid))


def pytest_sessionstart(session):
    session._setupstate = SetupState()


def pytest_sessionfinish(session):
    session._setupstate.teardown_all()


def pytest_runtest_protocol(item, nextitem):
    item.ihook.pytest_runtest_logstart(
        nodeid=item.nodeid, location=item.location,
    )
    runtestprotocol(item, nextitem=nextitem)
    item.ihook.pytest_runtest_logfinish(
        nodeid=item.nodeid, location=item.location,
    )
    return True


def runtestprotocol(item, log=True, nextitem=None):
    hasrequest = hasattr(item, "_request")
    if hasrequest and not item._request:
        item._initrequest()
    rep = call_and_report(item, "setup", log)
    reports = [rep]
    if rep.passed:
        if item.config.option.setupshow:
            show_test_item(item)
        if not item.config.option.setuponly:
            reports.append(call_and_report(item, "call", log))
    reports.append(call_and_report(item, "teardown", log,
                                   nextitem=nextitem))
    # after all teardown hooks have been called
    # want funcargs and request info to go away
    if hasrequest:
        item._request = False
        item.funcargs = None
    return reports


def show_test_item(item):
    """Show test function, parameters and the fixtures of the test item."""
    tw = item.config.get_terminal_writer()
    tw.line()
    tw.write(' ' * 8)
    tw.write(item._nodeid)
    used_fixtures = sorted(item._fixtureinfo.name2fixturedefs.keys())
    if used_fixtures:
        tw.write(' (fixtures used: {0})'.format(', '.join(used_fixtures)))


def pytest_runtest_setup(item):
    _update_current_test_var(item, 'setup')
    item.session._setupstate.prepare(item)


def pytest_runtest_call(item):
    _update_current_test_var(item, 'call')
    sys.last_type, sys.last_value, sys.last_traceback = (None, None, None)
    try:
        item.runtest()
    except Exception:
        # Store trace info to allow postmortem debugging
        type, value, tb = sys.exc_info()
        tb = tb.tb_next  # Skip *this* frame
        sys.last_type = type
        sys.last_value = value
        sys.last_traceback = tb
        del type, value, tb  # Get rid of these in this frame
        raise


def pytest_runtest_teardown(item, nextitem):
    _update_current_test_var(item, 'teardown')
    item.session._setupstate.teardown_exact(item, nextitem)
    _update_current_test_var(item, None)


def _update_current_test_var(item, when):
    """
    Update PYTEST_CURRENT_TEST to reflect the current item and stage.

    If ``when`` is None, delete PYTEST_CURRENT_TEST from the environment.
    """
    var_name = 'PYTEST_CURRENT_TEST'
    if when:
        value = '{0} ({1})'.format(item.nodeid, when)
        # don't allow null bytes on environment variables (see #2644, #2957)
        value = value.replace('\x00', '(null)')
        os.environ[var_name] = value
    else:
        os.environ.pop(var_name)


def pytest_report_teststatus(report):
    if report.when in ("setup", "teardown"):
        if report.failed:
            #      category, shortletter, verbose-word
            return "error", "E", "ERROR"
        elif report.skipped:
            return "skipped", "s", "SKIPPED"
        else:
            return "", "", ""


#
# Implementation

def call_and_report(item, when, log=True, **kwds):
    call = call_runtest_hook(item, when, **kwds)
    hook = item.ihook
    report = hook.pytest_runtest_makereport(item=item, call=call)
    if log:
        hook.pytest_runtest_logreport(report=report)
    if check_interactive_exception(call, report):
        hook.pytest_exception_interact(node=item, call=call, report=report)
    return report


def check_interactive_exception(call, report):
    return call.excinfo and not (
        hasattr(report, "wasxfail") or
        call.excinfo.errisinstance(skip.Exception) or
        call.excinfo.errisinstance(bdb.BdbQuit))


def call_runtest_hook(item, when, **kwds):
    hookname = "pytest_runtest_" + when
    ihook = getattr(item.ihook, hookname)
    return CallInfo(lambda: ihook(item=item, **kwds), when=when)


class CallInfo(object):
    """ Result/Exception info a function invocation. """
    #: None or ExceptionInfo object.
    excinfo = None

    def __init__(self, func, when):
        #: context of invocation: one of "setup", "call",
        #: "teardown", "memocollect"
        self.when = when
        self.start = time()
        try:
            self.result = func()
        except KeyboardInterrupt:
            self.stop = time()
            raise
        except:  # noqa
            self.excinfo = ExceptionInfo()
        self.stop = time()

    def __repr__(self):
        if self.excinfo:
            status = "exception: %s" % str(self.excinfo.value)
        else:
            status = "result: %r" % (self.result,)
        return "<CallInfo when=%r %s>" % (self.when, status)


def getslaveinfoline(node):
    try:
        return node._slaveinfocache
    except AttributeError:
        d = node.slaveinfo
        ver = "%s.%s.%s" % d['version_info'][:3]
        node._slaveinfocache = s = "[%s] %s -- Python %s %s" % (
            d['id'], d['sysplatform'], ver, d['executable'])
        return s


class BaseReport(object):

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def toterminal(self, out):
        if hasattr(self, 'node'):
            out.line(getslaveinfoline(self.node))

        longrepr = self.longrepr
        if longrepr is None:
            return

        if hasattr(longrepr, 'toterminal'):
            longrepr.toterminal(out)
        else:
            try:
                out.line(longrepr)
            except UnicodeEncodeError:
                out.line("<unprintable longrepr>")

    def get_sections(self, prefix):
        for name, content in self.sections:
            if name.startswith(prefix):
                yield prefix, content

    @property
    def longreprtext(self):
        """
        Read-only property that returns the full string representation
        of ``longrepr``.

        .. versionadded:: 3.0
        """
        tw = py.io.TerminalWriter(stringio=True)
        tw.hasmarkup = False
        self.toterminal(tw)
        exc = tw.stringio.getvalue()
        return exc.strip()

    @property
    def caplog(self):
        """Return captured log lines, if log capturing is enabled

        .. versionadded:: 3.5
        """
        return '\n'.join(content for (prefix, content) in self.get_sections('Captured log'))

    @property
    def capstdout(self):
        """Return captured text from stdout, if capturing is enabled

        .. versionadded:: 3.0
        """
        return ''.join(content for (prefix, content) in self.get_sections('Captured stdout'))

    @property
    def capstderr(self):
        """Return captured text from stderr, if capturing is enabled

        .. versionadded:: 3.0
        """
        return ''.join(content for (prefix, content) in self.get_sections('Captured stderr'))

    passed = property(lambda x: x.outcome == "passed")
    failed = property(lambda x: x.outcome == "failed")
    skipped = property(lambda x: x.outcome == "skipped")

    @property
    def fspath(self):
        return self.nodeid.split("::")[0]


def pytest_runtest_makereport(item, call):
    when = call.when
    duration = call.stop - call.start
    keywords = dict([(x, 1) for x in item.keywords])
    excinfo = call.excinfo
    sections = []
    if not call.excinfo:
        outcome = "passed"
        longrepr = None
    else:
        if not isinstance(excinfo, ExceptionInfo):
            outcome = "failed"
            longrepr = excinfo
        elif excinfo.errisinstance(skip.Exception):
            outcome = "skipped"
            r = excinfo._getreprcrash()
            longrepr = (str(r.path), r.lineno, r.message)
        else:
            outcome = "failed"
            if call.when == "call":
                longrepr = item.repr_failure(excinfo)
            else:  # exception in setup or teardown
                longrepr = item._repr_failure_py(excinfo,
                                                 style=item.config.option.tbstyle)
    for rwhen, key, content in item._report_sections:
        sections.append(("Captured %s %s" % (key, rwhen), content))
    return TestReport(item.nodeid, item.location,
                      keywords, outcome, longrepr, when,
                      sections, duration, user_properties=item.user_properties)


class TestReport(BaseReport):
    """ Basic test report object (also used for setup and teardown calls if
    they fail).
    """

    def __init__(self, nodeid, location, keywords, outcome,
                 longrepr, when, sections=(), duration=0, user_properties=(), **extra):
        #: normalized collection node id
        self.nodeid = nodeid

        #: a (filesystempath, lineno, domaininfo) tuple indicating the
        #: actual location of a test item - it might be different from the
        #: collected one e.g. if a method is inherited from a different module.
        self.location = location

        #: a name -> value dictionary containing all keywords and
        #: markers associated with a test invocation.
        self.keywords = keywords

        #: test outcome, always one of "passed", "failed", "skipped".
        self.outcome = outcome

        #: None or a failure representation.
        self.longrepr = longrepr

        #: one of 'setup', 'call', 'teardown' to indicate runtest phase.
        self.when = when

        #: user properties is a list of tuples (name, value) that holds user
        #: defined properties of the test
        self.user_properties = user_properties

        #: list of pairs ``(str, str)`` of extra information which needs to
        #: marshallable. Used by pytest to add captured text
        #: from ``stdout`` and ``stderr``, but may be used by other plugins
        #: to add arbitrary information to reports.
        self.sections = list(sections)

        #: time it took to run just the test
        self.duration = duration

        self.__dict__.update(extra)

    def __repr__(self):
        return "<TestReport %r when=%r outcome=%r>" % (
            self.nodeid, self.when, self.outcome)


class TeardownErrorReport(BaseReport):
    outcome = "failed"
    when = "teardown"

    def __init__(self, longrepr, **extra):
        self.longrepr = longrepr
        self.sections = []
        self.__dict__.update(extra)


def pytest_make_collect_report(collector):
    call = CallInfo(
        lambda: list(collector.collect()),
        'collect')
    longrepr = None
    if not call.excinfo:
        outcome = "passed"
    else:
        from _pytest import nose
        skip_exceptions = (Skipped,) + nose.get_skip_exceptions()
        if call.excinfo.errisinstance(skip_exceptions):
            outcome = "skipped"
            r = collector._repr_failure_py(call.excinfo, "line").reprcrash
            longrepr = (str(r.path), r.lineno, r.message)
        else:
            outcome = "failed"
            errorinfo = collector.repr_failure(call.excinfo)
            if not hasattr(errorinfo, "toterminal"):
                errorinfo = CollectErrorRepr(errorinfo)
            longrepr = errorinfo
    rep = CollectReport(collector.nodeid, outcome, longrepr,
                        getattr(call, 'result', None))
    rep.call = call  # see collect_one_node
    return rep


class CollectReport(BaseReport):
    def __init__(self, nodeid, outcome, longrepr, result,
                 sections=(), **extra):
        self.nodeid = nodeid
        self.outcome = outcome
        self.longrepr = longrepr
        self.result = result or []
        self.sections = list(sections)
        self.__dict__.update(extra)

    @property
    def location(self):
        return (self.fspath, None, self.fspath)

    def __repr__(self):
        return "<CollectReport %r lenresult=%s outcome=%r>" % (
            self.nodeid, len(self.result), self.outcome)


class CollectErrorRepr(TerminalRepr):
    def __init__(self, msg):
        self.longrepr = msg

    def toterminal(self, out):
        out.line(self.longrepr, red=True)


class SetupState(object):
    """ shared state for setting up/tearing down test items or collectors. """

    def __init__(self):
        self.stack = []
        self._finalizers = {}

    def addfinalizer(self, finalizer, colitem):
        """ attach a finalizer to the given colitem.
        if colitem is None, this will add a finalizer that
        is called at the end of teardown_all().
        """
        assert colitem and not isinstance(colitem, tuple)
        assert callable(finalizer)
        # assert colitem in self.stack  # some unit tests don't setup stack :/
        self._finalizers.setdefault(colitem, []).append(finalizer)

    def _pop_and_teardown(self):
        colitem = self.stack.pop()
        self._teardown_with_finalization(colitem)

    def _callfinalizers(self, colitem):
        finalizers = self._finalizers.pop(colitem, None)
        exc = None
        while finalizers:
            fin = finalizers.pop()
            try:
                fin()
            except TEST_OUTCOME:
                # XXX Only first exception will be seen by user,
                #     ideally all should be reported.
                if exc is None:
                    exc = sys.exc_info()
        if exc:
            py.builtin._reraise(*exc)

    def _teardown_with_finalization(self, colitem):
        self._callfinalizers(colitem)
        if hasattr(colitem, "teardown"):
            colitem.teardown()
        for colitem in self._finalizers:
            assert colitem is None or colitem in self.stack \
                or isinstance(colitem, tuple)

    def teardown_all(self):
        while self.stack:
            self._pop_and_teardown()
        for key in list(self._finalizers):
            self._teardown_with_finalization(key)
        assert not self._finalizers

    def teardown_exact(self, item, nextitem):
        needed_collectors = nextitem and nextitem.listchain() or []
        self._teardown_towards(needed_collectors)

    def _teardown_towards(self, needed_collectors):
        while self.stack:
            if self.stack == needed_collectors[:len(self.stack)]:
                break
            self._pop_and_teardown()

    def prepare(self, colitem):
        """ setup objects along the collector chain to the test-method
            and teardown previously setup objects."""
        needed_collectors = colitem.listchain()
        self._teardown_towards(needed_collectors)

        # check if the last collection node has raised an error
        for col in self.stack:
            if hasattr(col, '_prepare_exc'):
                py.builtin._reraise(*col._prepare_exc)
        for col in needed_collectors[len(self.stack):]:
            self.stack.append(col)
            try:
                col.setup()
            except TEST_OUTCOME:
                col._prepare_exc = sys.exc_info()
                raise


def collect_one_node(collector):
    ihook = collector.ihook
    ihook.pytest_collectstart(collector=collector)
    rep = ihook.pytest_make_collect_report(collector=collector)
    call = rep.__dict__.pop("call", None)
    if call and check_interactive_exception(call, rep):
        ihook.pytest_exception_interact(node=collector, call=call, report=rep)
    return rep
