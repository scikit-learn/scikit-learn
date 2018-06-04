""" terminal reporting of the full testing process.

This is a good source for looking at the various reporting hooks.
"""
from __future__ import absolute_import, division, print_function

import itertools
import platform
import sys
import time

import pluggy
import py
import six
from more_itertools import collapse

import pytest
from _pytest import nodes
from _pytest.main import EXIT_OK, EXIT_TESTSFAILED, EXIT_INTERRUPTED, \
    EXIT_USAGEERROR, EXIT_NOTESTSCOLLECTED


import argparse


class MoreQuietAction(argparse.Action):
    """
    a modified copy of the argparse count action which counts down and updates
    the legacy quiet attribute at the same time

    used to unify verbosity handling
    """
    def __init__(self,
                 option_strings,
                 dest,
                 default=None,
                 required=False,
                 help=None):
        super(MoreQuietAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            default=default,
            required=required,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        new_count = getattr(namespace, self.dest, 0) - 1
        setattr(namespace, self.dest, new_count)
        # todo Deprecate config.quiet
        namespace.quiet = getattr(namespace, 'quiet', 0) + 1


def pytest_addoption(parser):
    group = parser.getgroup("terminal reporting", "reporting", after="general")
    group._addoption('-v', '--verbose', action="count", default=0,
                     dest="verbose", help="increase verbosity."),
    group._addoption('-q', '--quiet', action=MoreQuietAction, default=0,
                     dest="verbose", help="decrease verbosity."),
    group._addoption("--verbosity", dest='verbose', type=int, default=0,
                     help="set verbosity")
    group._addoption('-r',
                     action="store", dest="reportchars", default='', metavar="chars",
                     help="show extra test summary info as specified by chars (f)ailed, "
                     "(E)error, (s)skipped, (x)failed, (X)passed, "
                     "(p)passed, (P)passed with output, (a)all except pP. "
                     "Warnings are displayed at all times except when "
                     "--disable-warnings is set")
    group._addoption('--disable-warnings', '--disable-pytest-warnings', default=False,
                     dest='disable_warnings', action='store_true',
                     help='disable warnings summary')
    group._addoption('-l', '--showlocals',
                     action="store_true", dest="showlocals", default=False,
                     help="show locals in tracebacks (disabled by default).")
    group._addoption('--tb', metavar="style",
                     action="store", dest="tbstyle", default='auto',
                     choices=['auto', 'long', 'short', 'no', 'line', 'native'],
                     help="traceback print mode (auto/long/short/line/native/no).")
    group._addoption('--show-capture',
                     action="store", dest="showcapture",
                     choices=['no', 'stdout', 'stderr', 'log', 'all'], default='all',
                     help="Controls how captured stdout/stderr/log is shown on failed tests. "
                          "Default is 'all'.")
    group._addoption('--fulltrace', '--full-trace',
                     action="store_true", default=False,
                     help="don't cut any tracebacks (default is to cut).")
    group._addoption('--color', metavar="color",
                     action="store", dest="color", default='auto',
                     choices=['yes', 'no', 'auto'],
                     help="color terminal output (yes/no/auto).")

    parser.addini("console_output_style",
                  help="console output: classic or with additional progress information (classic|progress).",
                  default='progress')


def pytest_configure(config):
    reporter = TerminalReporter(config, sys.stdout)
    config.pluginmanager.register(reporter, 'terminalreporter')
    if config.option.debug or config.option.traceconfig:
        def mywriter(tags, args):
            msg = " ".join(map(str, args))
            reporter.write_line("[traceconfig] " + msg)
        config.trace.root.setprocessor("pytest:config", mywriter)


def getreportopt(config):
    reportopts = ""
    reportchars = config.option.reportchars
    if not config.option.disable_warnings and 'w' not in reportchars:
        reportchars += 'w'
    elif config.option.disable_warnings and 'w' in reportchars:
        reportchars = reportchars.replace('w', '')
    if reportchars:
        for char in reportchars:
            if char not in reportopts and char != 'a':
                reportopts += char
            elif char == 'a':
                reportopts = 'fEsxXw'
    return reportopts


def pytest_report_teststatus(report):
    if report.passed:
        letter = "."
    elif report.skipped:
        letter = "s"
    elif report.failed:
        letter = "F"
        if report.when != "call":
            letter = "f"
    return report.outcome, letter, report.outcome.upper()


class WarningReport(object):
    """
    Simple structure to hold warnings information captured by ``pytest_logwarning``.
    """

    def __init__(self, code, message, nodeid=None, fslocation=None):
        """
        :param code: unused
        :param str message: user friendly message about the warning
        :param str|None nodeid: node id that generated the warning (see ``get_location``).
        :param tuple|py.path.local fslocation:
            file system location of the source of the warning (see ``get_location``).
        """
        self.code = code
        self.message = message
        self.nodeid = nodeid
        self.fslocation = fslocation

    def get_location(self, config):
        """
        Returns the more user-friendly information about the location
        of a warning, or None.
        """
        if self.nodeid:
            return self.nodeid
        if self.fslocation:
            if isinstance(self.fslocation, tuple) and len(self.fslocation) >= 2:
                filename, linenum = self.fslocation[:2]
                relpath = py.path.local(filename).relto(config.invocation_dir)
                return '%s:%s' % (relpath, linenum)
            else:
                return str(self.fslocation)
        return None


class TerminalReporter(object):
    def __init__(self, config, file=None):
        import _pytest.config
        self.config = config
        self.verbosity = self.config.option.verbose
        self.showheader = self.verbosity >= 0
        self.showfspath = self.verbosity >= 0
        self.showlongtestinfo = self.verbosity > 0
        self._numcollected = 0
        self._session = None

        self.stats = {}
        self.startdir = py.path.local()
        if file is None:
            file = sys.stdout
        self._tw = _pytest.config.create_terminal_writer(config, file)
        # self.writer will be deprecated in pytest-3.4
        self.writer = self._tw
        self._screen_width = self._tw.fullwidth
        self.currentfspath = None
        self.reportchars = getreportopt(config)
        self.hasmarkup = self._tw.hasmarkup
        self.isatty = file.isatty()
        self._progress_nodeids_reported = set()
        self._show_progress_info = self._determine_show_progress_info()

    def _determine_show_progress_info(self):
        """Return True if we should display progress information based on the current config"""
        # do not show progress if we are not capturing output (#3038)
        if self.config.getoption('capture') == 'no':
            return False
        # do not show progress if we are showing fixture setup/teardown
        if self.config.getoption('setupshow'):
            return False
        return self.config.getini('console_output_style') == 'progress'

    def hasopt(self, char):
        char = {'xfailed': 'x', 'skipped': 's'}.get(char, char)
        return char in self.reportchars

    def write_fspath_result(self, nodeid, res):
        fspath = self.config.rootdir.join(nodeid.split("::")[0])
        if fspath != self.currentfspath:
            if self.currentfspath is not None:
                self._write_progress_information_filling_space()
            self.currentfspath = fspath
            fspath = self.startdir.bestrelpath(fspath)
            self._tw.line()
            self._tw.write(fspath + " ")
        self._tw.write(res)

    def write_ensure_prefix(self, prefix, extra="", **kwargs):
        if self.currentfspath != prefix:
            self._tw.line()
            self.currentfspath = prefix
            self._tw.write(prefix)
        if extra:
            self._tw.write(extra, **kwargs)
            self.currentfspath = -2

    def ensure_newline(self):
        if self.currentfspath:
            self._tw.line()
            self.currentfspath = None

    def write(self, content, **markup):
        self._tw.write(content, **markup)

    def write_line(self, line, **markup):
        if not isinstance(line, six.text_type):
            line = six.text_type(line, errors="replace")
        self.ensure_newline()
        self._tw.line(line, **markup)

    def rewrite(self, line, **markup):
        """
        Rewinds the terminal cursor to the beginning and writes the given line.

        :kwarg erase: if True, will also add spaces until the full terminal width to ensure
            previous lines are properly erased.

        The rest of the keyword arguments are markup instructions.
        """
        erase = markup.pop('erase', False)
        if erase:
            fill_count = self._tw.fullwidth - len(line) - 1
            fill = ' ' * fill_count
        else:
            fill = ''
        line = str(line)
        self._tw.write("\r" + line + fill, **markup)

    def write_sep(self, sep, title=None, **markup):
        self.ensure_newline()
        self._tw.sep(sep, title, **markup)

    def section(self, title, sep="=", **kw):
        self._tw.sep(sep, title, **kw)

    def line(self, msg, **kw):
        self._tw.line(msg, **kw)

    def pytest_internalerror(self, excrepr):
        for line in six.text_type(excrepr).split("\n"):
            self.write_line("INTERNALERROR> " + line)
        return 1

    def pytest_logwarning(self, code, fslocation, message, nodeid):
        warnings = self.stats.setdefault("warnings", [])
        warning = WarningReport(code=code, fslocation=fslocation,
                                message=message, nodeid=nodeid)
        warnings.append(warning)

    def pytest_plugin_registered(self, plugin):
        if self.config.option.traceconfig:
            msg = "PLUGIN registered: %s" % (plugin,)
            # XXX this event may happen during setup/teardown time
            #     which unfortunately captures our output here
            #     which garbles our output if we use self.write_line
            self.write_line(msg)

    def pytest_deselected(self, items):
        self.stats.setdefault('deselected', []).extend(items)

    def pytest_runtest_logstart(self, nodeid, location):
        # ensure that the path is printed before the
        # 1st test of a module starts running
        if self.showlongtestinfo:
            line = self._locationline(nodeid, *location)
            self.write_ensure_prefix(line, "")
        elif self.showfspath:
            fsid = nodeid.split("::")[0]
            self.write_fspath_result(fsid, "")

    def pytest_runtest_logreport(self, report):
        rep = report
        res = self.config.hook.pytest_report_teststatus(report=rep)
        cat, letter, word = res
        if isinstance(word, tuple):
            word, markup = word
        else:
            markup = None
        self.stats.setdefault(cat, []).append(rep)
        self._tests_ran = True
        if not letter and not word:
            # probably passed setup/teardown
            return
        running_xdist = hasattr(rep, 'node')
        if self.verbosity <= 0:
            if not running_xdist and self.showfspath:
                self.write_fspath_result(rep.nodeid, letter)
            else:
                self._tw.write(letter)
        else:
            self._progress_nodeids_reported.add(rep.nodeid)
            if markup is None:
                if rep.passed:
                    markup = {'green': True}
                elif rep.failed:
                    markup = {'red': True}
                elif rep.skipped:
                    markup = {'yellow': True}
                else:
                    markup = {}
            line = self._locationline(rep.nodeid, *rep.location)
            if not running_xdist:
                self.write_ensure_prefix(line, word, **markup)
                if self._show_progress_info:
                    self._write_progress_information_filling_space()
            else:
                self.ensure_newline()
                self._tw.write("[%s]" % rep.node.gateway.id)
                if self._show_progress_info:
                    self._tw.write(self._get_progress_information_message() + " ", cyan=True)
                else:
                    self._tw.write(' ')
                self._tw.write(word, **markup)
                self._tw.write(" " + line)
                self.currentfspath = -2

    def pytest_runtest_logfinish(self, nodeid):
        if self.verbosity <= 0 and self._show_progress_info:
            self._progress_nodeids_reported.add(nodeid)
            last_item = len(self._progress_nodeids_reported) == self._session.testscollected
            if last_item:
                self._write_progress_information_filling_space()
            else:
                past_edge = self._tw.chars_on_current_line + self._PROGRESS_LENGTH + 1 >= self._screen_width
                if past_edge:
                    msg = self._get_progress_information_message()
                    self._tw.write(msg + '\n', cyan=True)

    _PROGRESS_LENGTH = len(' [100%]')

    def _get_progress_information_message(self):
        if self.config.getoption('capture') == 'no':
            return ''
        collected = self._session.testscollected
        if collected:
            progress = len(self._progress_nodeids_reported) * 100 // collected
            return ' [{:3d}%]'.format(progress)
        return ' [100%]'

    def _write_progress_information_filling_space(self):
        msg = self._get_progress_information_message()
        fill = ' ' * (self._tw.fullwidth - self._tw.chars_on_current_line - len(msg) - 1)
        self.write(fill + msg, cyan=True)

    def pytest_collection(self):
        if not self.isatty and self.config.option.verbose >= 1:
            self.write("collecting ... ", bold=True)

    def pytest_collectreport(self, report):
        if report.failed:
            self.stats.setdefault("error", []).append(report)
        elif report.skipped:
            self.stats.setdefault("skipped", []).append(report)
        items = [x for x in report.result if isinstance(x, pytest.Item)]
        self._numcollected += len(items)
        if self.isatty:
            # self.write_fspath_result(report.nodeid, 'E')
            self.report_collect()

    def report_collect(self, final=False):
        if self.config.option.verbose < 0:
            return

        errors = len(self.stats.get('error', []))
        skipped = len(self.stats.get('skipped', []))
        deselected = len(self.stats.get('deselected', []))
        if final:
            line = "collected "
        else:
            line = "collecting "
        line += str(self._numcollected) + " item" + ('' if self._numcollected == 1 else 's')
        if errors:
            line += " / %d errors" % errors
        if deselected:
            line += " / %d deselected" % deselected
        if skipped:
            line += " / %d skipped" % skipped
        if self.isatty:
            self.rewrite(line, bold=True, erase=True)
            if final:
                self.write('\n')
        else:
            self.write_line(line)

    @pytest.hookimpl(trylast=True)
    def pytest_collection_modifyitems(self):
        self.report_collect(True)

    @pytest.hookimpl(trylast=True)
    def pytest_sessionstart(self, session):
        self._session = session
        self._sessionstarttime = time.time()
        if not self.showheader:
            return
        self.write_sep("=", "test session starts", bold=True)
        verinfo = platform.python_version()
        msg = "platform %s -- Python %s" % (sys.platform, verinfo)
        if hasattr(sys, 'pypy_version_info'):
            verinfo = ".".join(map(str, sys.pypy_version_info[:3]))
            msg += "[pypy-%s-%s]" % (verinfo, sys.pypy_version_info[3])
        msg += ", pytest-%s, py-%s, pluggy-%s" % (
               pytest.__version__, py.__version__, pluggy.__version__)
        if self.verbosity > 0 or self.config.option.debug or \
           getattr(self.config.option, 'pastebin', None):
            msg += " -- " + str(sys.executable)
        self.write_line(msg)
        lines = self.config.hook.pytest_report_header(
            config=self.config, startdir=self.startdir)
        self._write_report_lines_from_hooks(lines)

    def _write_report_lines_from_hooks(self, lines):
        lines.reverse()
        for line in collapse(lines):
            self.write_line(line)

    def pytest_report_header(self, config):
        inifile = ""
        if config.inifile:
            inifile = " " + config.rootdir.bestrelpath(config.inifile)
        lines = ["rootdir: %s, inifile:%s" % (config.rootdir, inifile)]

        plugininfo = config.pluginmanager.list_plugin_distinfo()
        if plugininfo:

            lines.append(
                "plugins: %s" % ", ".join(_plugin_nameversions(plugininfo)))
        return lines

    def pytest_collection_finish(self, session):
        if self.config.option.collectonly:
            self._printcollecteditems(session.items)
            if self.stats.get('failed'):
                self._tw.sep("!", "collection failures")
                for rep in self.stats.get('failed'):
                    rep.toterminal(self._tw)
                return 1
            return 0
        lines = self.config.hook.pytest_report_collectionfinish(
            config=self.config, startdir=self.startdir, items=session.items)
        self._write_report_lines_from_hooks(lines)

    def _printcollecteditems(self, items):
        # to print out items and their parent collectors
        # we take care to leave out Instances aka ()
        # because later versions are going to get rid of them anyway
        if self.config.option.verbose < 0:
            if self.config.option.verbose < -1:
                counts = {}
                for item in items:
                    name = item.nodeid.split('::', 1)[0]
                    counts[name] = counts.get(name, 0) + 1
                for name, count in sorted(counts.items()):
                    self._tw.line("%s: %d" % (name, count))
            else:
                for item in items:
                    nodeid = item.nodeid
                    nodeid = nodeid.replace("::()::", "::")
                    self._tw.line(nodeid)
            return
        stack = []
        indent = ""
        for item in items:
            needed_collectors = item.listchain()[1:]  # strip root node
            while stack:
                if stack == needed_collectors[:len(stack)]:
                    break
                stack.pop()
            for col in needed_collectors[len(stack):]:
                stack.append(col)
                # if col.name == "()":
                #    continue
                indent = (len(stack) - 1) * "  "
                self._tw.line("%s%s" % (indent, col))

    @pytest.hookimpl(hookwrapper=True)
    def pytest_sessionfinish(self, exitstatus):
        outcome = yield
        outcome.get_result()
        self._tw.line("")
        summary_exit_codes = (
            EXIT_OK, EXIT_TESTSFAILED, EXIT_INTERRUPTED, EXIT_USAGEERROR,
            EXIT_NOTESTSCOLLECTED)
        if exitstatus in summary_exit_codes:
            self.config.hook.pytest_terminal_summary(terminalreporter=self,
                                                     exitstatus=exitstatus)
        if exitstatus == EXIT_INTERRUPTED:
            self._report_keyboardinterrupt()
            del self._keyboardinterrupt_memo
        self.summary_stats()

    @pytest.hookimpl(hookwrapper=True)
    def pytest_terminal_summary(self):
        self.summary_errors()
        self.summary_failures()
        yield
        self.summary_warnings()
        self.summary_passes()

    def pytest_keyboard_interrupt(self, excinfo):
        self._keyboardinterrupt_memo = excinfo.getrepr(funcargs=True)

    def pytest_unconfigure(self):
        if hasattr(self, '_keyboardinterrupt_memo'):
            self._report_keyboardinterrupt()

    def _report_keyboardinterrupt(self):
        excrepr = self._keyboardinterrupt_memo
        msg = excrepr.reprcrash.message
        self.write_sep("!", msg)
        if "KeyboardInterrupt" in msg:
            if self.config.option.fulltrace:
                excrepr.toterminal(self._tw)
            else:
                self._tw.line("to show a full traceback on KeyboardInterrupt use --fulltrace", yellow=True)
                excrepr.reprcrash.toterminal(self._tw)

    def _locationline(self, nodeid, fspath, lineno, domain):
        def mkrel(nodeid):
            line = self.config.cwd_relative_nodeid(nodeid)
            if domain and line.endswith(domain):
                line = line[:-len(domain)]
                values = domain.split("[")
                values[0] = values[0].replace('.', '::')  # don't replace '.' in params
                line += "[".join(values)
            return line
        # collect_fspath comes from testid which has a "/"-normalized path

        if fspath:
            res = mkrel(nodeid).replace("::()", "")  # parens-normalization
            if nodeid.split("::")[0] != fspath.replace("\\", nodes.SEP):
                res += " <- " + self.startdir.bestrelpath(fspath)
        else:
            res = "[location]"
        return res + " "

    def _getfailureheadline(self, rep):
        if hasattr(rep, 'location'):
            fspath, lineno, domain = rep.location
            return domain
        else:
            return "test session"  # XXX?

    def _getcrashline(self, rep):
        try:
            return str(rep.longrepr.reprcrash)
        except AttributeError:
            try:
                return str(rep.longrepr)[:50]
            except AttributeError:
                return ""

    #
    # summaries for sessionfinish
    #
    def getreports(self, name):
        values = []
        for x in self.stats.get(name, []):
            if not hasattr(x, '_pdbshown'):
                values.append(x)
        return values

    def summary_warnings(self):
        if self.hasopt("w"):
            all_warnings = self.stats.get("warnings")
            if not all_warnings:
                return

            grouped = itertools.groupby(all_warnings, key=lambda wr: wr.get_location(self.config))

            self.write_sep("=", "warnings summary", yellow=True, bold=False)
            for location, warning_records in grouped:
                self._tw.line(str(location) or '<undetermined location>')
                for w in warning_records:
                    lines = w.message.splitlines()
                    indented = '\n'.join('  ' + x for x in lines)
                    self._tw.line(indented)
                self._tw.line()
            self._tw.line('-- Docs: http://doc.pytest.org/en/latest/warnings.html')

    def summary_passes(self):
        if self.config.option.tbstyle != "no":
            if self.hasopt("P"):
                reports = self.getreports('passed')
                if not reports:
                    return
                self.write_sep("=", "PASSES")
                for rep in reports:
                    msg = self._getfailureheadline(rep)
                    self.write_sep("_", msg)
                    self._outrep_summary(rep)

    def print_teardown_sections(self, rep):
        for secname, content in rep.sections:
            if 'teardown' in secname:
                self._tw.sep('-', secname)
                if content[-1:] == "\n":
                    content = content[:-1]
                self._tw.line(content)

    def summary_failures(self):
        if self.config.option.tbstyle != "no":
            reports = self.getreports('failed')
            if not reports:
                return
            self.write_sep("=", "FAILURES")
            for rep in reports:
                if self.config.option.tbstyle == "line":
                    line = self._getcrashline(rep)
                    self.write_line(line)
                else:
                    msg = self._getfailureheadline(rep)
                    markup = {'red': True, 'bold': True}
                    self.write_sep("_", msg, **markup)
                    self._outrep_summary(rep)
                    for report in self.getreports(''):
                        if report.nodeid == rep.nodeid and report.when == 'teardown':
                            self.print_teardown_sections(report)

    def summary_errors(self):
        if self.config.option.tbstyle != "no":
            reports = self.getreports('error')
            if not reports:
                return
            self.write_sep("=", "ERRORS")
            for rep in self.stats['error']:
                msg = self._getfailureheadline(rep)
                if not hasattr(rep, 'when'):
                    # collect
                    msg = "ERROR collecting " + msg
                elif rep.when == "setup":
                    msg = "ERROR at setup of " + msg
                elif rep.when == "teardown":
                    msg = "ERROR at teardown of " + msg
                self.write_sep("_", msg)
                self._outrep_summary(rep)

    def _outrep_summary(self, rep):
        rep.toterminal(self._tw)
        showcapture = self.config.option.showcapture
        if showcapture == 'no':
            return
        for secname, content in rep.sections:
            if showcapture != 'all' and showcapture not in secname:
                continue
            self._tw.sep("-", secname)
            if content[-1:] == "\n":
                content = content[:-1]
            self._tw.line(content)

    def summary_stats(self):
        session_duration = time.time() - self._sessionstarttime
        (line, color) = build_summary_stats_line(self.stats)
        msg = "%s in %.2f seconds" % (line, session_duration)
        markup = {color: True, 'bold': True}

        if self.verbosity >= 0:
            self.write_sep("=", msg, **markup)
        if self.verbosity == -1:
            self.write_line(msg, **markup)


def repr_pythonversion(v=None):
    if v is None:
        v = sys.version_info
    try:
        return "%s.%s.%s-%s-%s" % v
    except (TypeError, ValueError):
        return str(v)


def build_summary_stats_line(stats):
    keys = ("failed passed skipped deselected "
            "xfailed xpassed warnings error").split()
    unknown_key_seen = False
    for key in stats.keys():
        if key not in keys:
            if key:  # setup/teardown reports have an empty key, ignore them
                keys.append(key)
                unknown_key_seen = True
    parts = []
    for key in keys:
        val = stats.get(key, None)
        if val:
            parts.append("%d %s" % (len(val), key))

    if parts:
        line = ", ".join(parts)
    else:
        line = "no tests ran"

    if 'failed' in stats or 'error' in stats:
        color = 'red'
    elif 'warnings' in stats or unknown_key_seen:
        color = 'yellow'
    elif 'passed' in stats:
        color = 'green'
    else:
        color = 'yellow'

    return (line, color)


def _plugin_nameversions(plugininfo):
    values = []
    for plugin, dist in plugininfo:
        # gets us name and version!
        name = '{dist.project_name}-{dist.version}'.format(dist=dist)
        # questionable convenience, but it keeps things short
        if name.startswith("pytest-"):
            name = name[7:]
        # we decided to print python package names
        # they can have more than one plugin
        if name not in values:
            values.append(name)
    return values
