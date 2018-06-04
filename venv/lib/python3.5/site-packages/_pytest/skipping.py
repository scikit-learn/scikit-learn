""" support for skip/xfail functions and markers. """
from __future__ import absolute_import, division, print_function

from _pytest.config import hookimpl
from _pytest.mark import MarkInfo, MarkDecorator
from _pytest.mark.evaluate import MarkEvaluator
from _pytest.outcomes import fail, skip, xfail


def pytest_addoption(parser):
    group = parser.getgroup("general")
    group.addoption('--runxfail',
                    action="store_true", dest="runxfail", default=False,
                    help="run tests even if they are marked xfail")

    parser.addini("xfail_strict",
                  "default for the strict parameter of xfail "
                  "markers when not given explicitly (default: False)",
                  default=False,
                  type="bool")


def pytest_configure(config):
    if config.option.runxfail:
        # yay a hack
        import pytest
        old = pytest.xfail
        config._cleanup.append(lambda: setattr(pytest, "xfail", old))

        def nop(*args, **kwargs):
            pass

        nop.Exception = xfail.Exception
        setattr(pytest, "xfail", nop)

    config.addinivalue_line("markers",
                            "skip(reason=None): skip the given test function with an optional reason. "
                            "Example: skip(reason=\"no way of currently testing this\") skips the "
                            "test."
                            )
    config.addinivalue_line("markers",
                            "skipif(condition): skip the given test function if eval(condition) "
                            "results in a True value.  Evaluation happens within the "
                            "module global context. Example: skipif('sys.platform == \"win32\"') "
                            "skips the test if we are on the win32 platform. see "
                            "http://pytest.org/latest/skipping.html"
                            )
    config.addinivalue_line("markers",
                            "xfail(condition, reason=None, run=True, raises=None, strict=False): "
                            "mark the test function as an expected failure if eval(condition) "
                            "has a True value. Optionally specify a reason for better reporting "
                            "and run=False if you don't even want to execute the test function. "
                            "If only specific exception(s) are expected, you can list them in "
                            "raises, and if the test fails in other ways, it will be reported as "
                            "a true failure. See http://pytest.org/latest/skipping.html"
                            )


@hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    # Check if skip or skipif are specified as pytest marks
    item._skipped_by_mark = False
    skipif_info = item.keywords.get('skipif')
    if isinstance(skipif_info, (MarkInfo, MarkDecorator)):
        eval_skipif = MarkEvaluator(item, 'skipif')
        if eval_skipif.istrue():
            item._skipped_by_mark = True
            skip(eval_skipif.getexplanation())

    skip_info = item.keywords.get('skip')
    if isinstance(skip_info, (MarkInfo, MarkDecorator)):
        item._skipped_by_mark = True
        if 'reason' in skip_info.kwargs:
            skip(skip_info.kwargs['reason'])
        elif skip_info.args:
            skip(skip_info.args[0])
        else:
            skip("unconditional skip")

    item._evalxfail = MarkEvaluator(item, 'xfail')
    check_xfail_no_run(item)


@hookimpl(hookwrapper=True)
def pytest_pyfunc_call(pyfuncitem):
    check_xfail_no_run(pyfuncitem)
    outcome = yield
    passed = outcome.excinfo is None
    if passed:
        check_strict_xfail(pyfuncitem)


def check_xfail_no_run(item):
    """check xfail(run=False)"""
    if not item.config.option.runxfail:
        evalxfail = item._evalxfail
        if evalxfail.istrue():
            if not evalxfail.get('run', True):
                xfail("[NOTRUN] " + evalxfail.getexplanation())


def check_strict_xfail(pyfuncitem):
    """check xfail(strict=True) for the given PASSING test"""
    evalxfail = pyfuncitem._evalxfail
    if evalxfail.istrue():
        strict_default = pyfuncitem.config.getini('xfail_strict')
        is_strict_xfail = evalxfail.get('strict', strict_default)
        if is_strict_xfail:
            del pyfuncitem._evalxfail
            explanation = evalxfail.getexplanation()
            fail('[XPASS(strict)] ' + explanation, pytrace=False)


@hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    evalxfail = getattr(item, '_evalxfail', None)
    # unitttest special case, see setting of _unexpectedsuccess
    if hasattr(item, '_unexpectedsuccess') and rep.when == "call":
        from _pytest.compat import _is_unittest_unexpected_success_a_failure
        if item._unexpectedsuccess:
            rep.longrepr = "Unexpected success: {0}".format(item._unexpectedsuccess)
        else:
            rep.longrepr = "Unexpected success"
        if _is_unittest_unexpected_success_a_failure():
            rep.outcome = "failed"
        else:
            rep.outcome = "passed"
            rep.wasxfail = rep.longrepr
    elif item.config.option.runxfail:
        pass  # don't interefere
    elif call.excinfo and call.excinfo.errisinstance(xfail.Exception):
        rep.wasxfail = "reason: " + call.excinfo.value.msg
        rep.outcome = "skipped"
    elif evalxfail and not rep.skipped and evalxfail.wasvalid() and \
            evalxfail.istrue():
        if call.excinfo:
            if evalxfail.invalidraise(call.excinfo.value):
                rep.outcome = "failed"
            else:
                rep.outcome = "skipped"
                rep.wasxfail = evalxfail.getexplanation()
        elif call.when == "call":
            strict_default = item.config.getini('xfail_strict')
            is_strict_xfail = evalxfail.get('strict', strict_default)
            explanation = evalxfail.getexplanation()
            if is_strict_xfail:
                rep.outcome = "failed"
                rep.longrepr = "[XPASS(strict)] {0}".format(explanation)
            else:
                rep.outcome = "passed"
                rep.wasxfail = explanation
    elif getattr(item, '_skipped_by_mark', False) and rep.skipped and type(rep.longrepr) is tuple:
        # skipped by mark.skipif; change the location of the failure
        # to point to the item definition, otherwise it will display
        # the location of where the skip exception was raised within pytest
        filename, line, reason = rep.longrepr
        filename, line = item.location[:2]
        rep.longrepr = filename, line, reason


# called by terminalreporter progress reporting


def pytest_report_teststatus(report):
    if hasattr(report, "wasxfail"):
        if report.skipped:
            return "xfailed", "x", "xfail"
        elif report.passed:
            return "xpassed", "X", ("XPASS", {'yellow': True})


# called by the terminalreporter instance/plugin


def pytest_terminal_summary(terminalreporter):
    tr = terminalreporter
    if not tr.reportchars:
        # for name in "xfailed skipped failed xpassed":
        #    if not tr.stats.get(name, 0):
        #        tr.write_line("HINT: use '-r' option to see extra "
        #              "summary info about tests")
        #        break
        return

    lines = []
    for char in tr.reportchars:
        action = REPORTCHAR_ACTIONS.get(char, lambda tr, lines: None)
        action(terminalreporter, lines)

    if lines:
        tr._tw.sep("=", "short test summary info")
        for line in lines:
            tr._tw.line(line)


def show_simple(terminalreporter, lines, stat, format):
    failed = terminalreporter.stats.get(stat)
    if failed:
        for rep in failed:
            pos = terminalreporter.config.cwd_relative_nodeid(rep.nodeid)
            lines.append(format % (pos,))


def show_xfailed(terminalreporter, lines):
    xfailed = terminalreporter.stats.get("xfailed")
    if xfailed:
        for rep in xfailed:
            pos = terminalreporter.config.cwd_relative_nodeid(rep.nodeid)
            reason = rep.wasxfail
            lines.append("XFAIL %s" % (pos,))
            if reason:
                lines.append("  " + str(reason))


def show_xpassed(terminalreporter, lines):
    xpassed = terminalreporter.stats.get("xpassed")
    if xpassed:
        for rep in xpassed:
            pos = terminalreporter.config.cwd_relative_nodeid(rep.nodeid)
            reason = rep.wasxfail
            lines.append("XPASS %s %s" % (pos, reason))


def folded_skips(skipped):
    d = {}
    for event in skipped:
        key = event.longrepr
        assert len(key) == 3, (event, key)
        keywords = getattr(event, 'keywords', {})
        # folding reports with global pytestmark variable
        # this is workaround, because for now we cannot identify the scope of a skip marker
        # TODO: revisit after marks scope would be fixed
        when = getattr(event, 'when', None)
        if when == 'setup' and 'skip' in keywords and 'pytestmark' not in keywords:
            key = (key[0], None, key[2])
        d.setdefault(key, []).append(event)
    values = []
    for key, events in d.items():
        values.append((len(events),) + key)
    return values


def show_skipped(terminalreporter, lines):
    tr = terminalreporter
    skipped = tr.stats.get('skipped', [])
    if skipped:
        # if not tr.hasopt('skipped'):
        #    tr.write_line(
        #        "%d skipped tests, specify -rs for more info" %
        #        len(skipped))
        #    return
        fskips = folded_skips(skipped)
        if fskips:
            # tr.write_sep("_", "skipped test summary")
            for num, fspath, lineno, reason in fskips:
                if reason.startswith("Skipped: "):
                    reason = reason[9:]
                if lineno is not None:
                    lines.append(
                        "SKIP [%d] %s:%d: %s" %
                        (num, fspath, lineno + 1, reason))
                else:
                    lines.append(
                        "SKIP [%d] %s: %s" %
                        (num, fspath, reason))


def shower(stat, format):
    def show_(terminalreporter, lines):
        return show_simple(terminalreporter, lines, stat, format)

    return show_


REPORTCHAR_ACTIONS = {
    'x': show_xfailed,
    'X': show_xpassed,
    'f': shower('failed', "FAIL %s"),
    'F': shower('failed', "FAIL %s"),
    's': show_skipped,
    'S': show_skipped,
    'p': shower('passed', "PASSED %s"),
    'E': shower('error', "ERROR %s")

}
