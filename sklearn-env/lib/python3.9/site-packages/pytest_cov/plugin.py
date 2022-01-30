"""Coverage plugin for pytest."""
import argparse
import os
import warnings

import coverage
import pytest

from . import compat
from . import embed


class CoverageError(Exception):
    """Indicates that our coverage is too low"""


class PytestCovWarning(pytest.PytestWarning):
    """
    The base for all pytest-cov warnings, never raised directly
    """


class CovDisabledWarning(PytestCovWarning):
    """Indicates that Coverage was manually disabled"""


class CovReportWarning(PytestCovWarning):
    """Indicates that we failed to generate a report"""


def validate_report(arg):
    file_choices = ['annotate', 'html', 'xml']
    term_choices = ['term', 'term-missing']
    term_modifier_choices = ['skip-covered']
    all_choices = term_choices + file_choices
    values = arg.split(":", 1)
    report_type = values[0]
    if report_type not in all_choices + ['']:
        msg = f'invalid choice: "{arg}" (choose from "{all_choices}")'
        raise argparse.ArgumentTypeError(msg)

    if len(values) == 1:
        return report_type, None

    report_modifier = values[1]
    if report_type in term_choices and report_modifier in term_modifier_choices:
        return report_type, report_modifier

    if report_type not in file_choices:
        msg = 'output specifier not supported for: "{}" (choose from "{}")'.format(arg,
                                                                                   file_choices)
        raise argparse.ArgumentTypeError(msg)

    return values


def validate_fail_under(num_str):
    try:
        value = int(num_str)
    except ValueError:
        try:
            value = float(num_str)
        except ValueError:
            raise argparse.ArgumentTypeError('An integer or float value is required.')
    if value > 100:
        raise argparse.ArgumentTypeError('Your desire for over-achievement is admirable but misplaced. '
                                         'The maximum value is 100. Perhaps write more integration tests?')
    return value


def validate_context(arg):
    if coverage.version_info <= (5, 0):
        raise argparse.ArgumentTypeError('Contexts are only supported with coverage.py >= 5.x')
    if arg != "test":
        raise argparse.ArgumentTypeError('The only supported value is "test".')
    return arg


class StoreReport(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        report_type, file = values
        namespace.cov_report[report_type] = file


def pytest_addoption(parser):
    """Add options to control coverage."""

    group = parser.getgroup(
        'cov', 'coverage reporting with distributed testing support')
    group.addoption('--cov', action='append', default=[], metavar='SOURCE',
                    nargs='?', const=True, dest='cov_source',
                    help='Path or package name to measure during execution (multi-allowed). '
                         'Use --cov= to not do any source filtering and record everything.')
    group.addoption('--cov-reset', action='store_const', const=[], dest='cov_source',
                    help='Reset cov sources accumulated in options so far. ')
    group.addoption('--cov-report', action=StoreReport, default={},
                    metavar='TYPE', type=validate_report,
                    help='Type of report to generate: term, term-missing, '
                         'annotate, html, xml (multi-allowed). '
                         'term, term-missing may be followed by ":skip-covered". '
                         'annotate, html and xml may be followed by ":DEST" '
                         'where DEST specifies the output location. '
                         'Use --cov-report= to not generate any output.')
    group.addoption('--cov-config', action='store', default='.coveragerc',
                    metavar='PATH',
                    help='Config file for coverage. Default: .coveragerc')
    group.addoption('--no-cov-on-fail', action='store_true', default=False,
                    help='Do not report coverage if test run fails. '
                         'Default: False')
    group.addoption('--no-cov', action='store_true', default=False,
                    help='Disable coverage report completely (useful for debuggers). '
                         'Default: False')
    group.addoption('--cov-fail-under', action='store', metavar='MIN',
                    type=validate_fail_under,
                    help='Fail if the total coverage is less than MIN.')
    group.addoption('--cov-append', action='store_true', default=False,
                    help='Do not delete coverage but append to current. '
                         'Default: False')
    group.addoption('--cov-branch', action='store_true', default=None,
                    help='Enable branch coverage.')
    group.addoption('--cov-context', action='store', metavar='CONTEXT',
                    type=validate_context,
                    help='Dynamic contexts to use. "test" for now.')


def _prepare_cov_source(cov_source):
    """
    Prepare cov_source so that:

     --cov --cov=foobar is equivalent to --cov (cov_source=None)
     --cov=foo --cov=bar is equivalent to cov_source=['foo', 'bar']
    """
    return None if True in cov_source else [path for path in cov_source if path is not True]


@pytest.mark.tryfirst
def pytest_load_initial_conftests(early_config, parser, args):
    options = early_config.known_args_namespace
    no_cov = options.no_cov_should_warn = False
    for arg in args:
        arg = str(arg)
        if arg == '--no-cov':
            no_cov = True
        elif arg.startswith('--cov') and no_cov:
            options.no_cov_should_warn = True
            break

    if early_config.known_args_namespace.cov_source:
        plugin = CovPlugin(options, early_config.pluginmanager)
        early_config.pluginmanager.register(plugin, '_cov')


class CovPlugin:
    """Use coverage package to produce code coverage reports.

    Delegates all work to a particular implementation based on whether
    this test process is centralised, a distributed master or a
    distributed worker.
    """

    def __init__(self, options, pluginmanager, start=True, no_cov_should_warn=False):
        """Creates a coverage pytest plugin.

        We read the rc file that coverage uses to get the data file
        name.  This is needed since we give coverage through it's API
        the data file name.
        """

        # Our implementation is unknown at this time.
        self.pid = None
        self.cov_controller = None
        self.cov_report = compat.StringIO()
        self.cov_total = None
        self.failed = False
        self._started = False
        self._start_path = None
        self._disabled = False
        self.options = options

        is_dist = (getattr(options, 'numprocesses', False) or
                   getattr(options, 'distload', False) or
                   getattr(options, 'dist', 'no') != 'no')
        if getattr(options, 'no_cov', False):
            self._disabled = True
            return

        if not self.options.cov_report:
            self.options.cov_report = ['term']
        elif len(self.options.cov_report) == 1 and '' in self.options.cov_report:
            self.options.cov_report = {}
        self.options.cov_source = _prepare_cov_source(self.options.cov_source)

        # import engine lazily here to avoid importing
        # it for unit tests that don't need it
        from . import engine

        if is_dist and start:
            self.start(engine.DistMaster)
        elif start:
            self.start(engine.Central)

        # worker is started in pytest hook

    def start(self, controller_cls, config=None, nodeid=None):

        if config is None:
            # fake config option for engine
            class Config:
                option = self.options

            config = Config()

        self.cov_controller = controller_cls(
            self.options.cov_source,
            self.options.cov_report,
            self.options.cov_config,
            self.options.cov_append,
            self.options.cov_branch,
            config,
            nodeid
        )
        self.cov_controller.start()
        self._started = True
        self._start_path = os.getcwd()
        cov_config = self.cov_controller.cov.config
        if self.options.cov_fail_under is None and hasattr(cov_config, 'fail_under'):
            self.options.cov_fail_under = cov_config.fail_under

    def _is_worker(self, session):
        return getattr(session.config, 'workerinput', None) is not None

    def pytest_sessionstart(self, session):
        """At session start determine our implementation and delegate to it."""

        if self.options.no_cov:
            # Coverage can be disabled because it does not cooperate with debuggers well.
            self._disabled = True
            return

        # import engine lazily here to avoid importing
        # it for unit tests that don't need it
        from . import engine

        self.pid = os.getpid()
        if self._is_worker(session):
            nodeid = (
                session.config.workerinput.get('workerid', getattr(session, 'nodeid'))
            )
            self.start(engine.DistWorker, session.config, nodeid)
        elif not self._started:
            self.start(engine.Central)

        if self.options.cov_context == 'test':
            session.config.pluginmanager.register(TestContextPlugin(self.cov_controller.cov), '_cov_contexts')

    def pytest_configure_node(self, node):
        """Delegate to our implementation.

        Mark this hook as optional in case xdist is not installed.
        """
        if not self._disabled:
            self.cov_controller.configure_node(node)
    pytest_configure_node.optionalhook = True

    def pytest_testnodedown(self, node, error):
        """Delegate to our implementation.

        Mark this hook as optional in case xdist is not installed.
        """
        if not self._disabled:
            self.cov_controller.testnodedown(node, error)
    pytest_testnodedown.optionalhook = True

    def _should_report(self):
        return not (self.failed and self.options.no_cov_on_fail)

    def _failed_cov_total(self):
        cov_fail_under = self.options.cov_fail_under
        return cov_fail_under is not None and self.cov_total < cov_fail_under

    # we need to wrap pytest_runtestloop. by the time pytest_sessionfinish
    # runs, it's too late to set testsfailed
    @compat.hookwrapper
    def pytest_runtestloop(self, session):
        yield

        if self._disabled:
            return

        compat_session = compat.SessionWrapper(session)

        self.failed = bool(compat_session.testsfailed)
        if self.cov_controller is not None:
            self.cov_controller.finish()

        if not self._is_worker(session) and self._should_report():

            # import coverage lazily here to avoid importing
            # it for unit tests that don't need it
            from coverage.misc import CoverageException

            try:
                self.cov_total = self.cov_controller.summary(self.cov_report)
            except CoverageException as exc:
                message = 'Failed to generate report: %s\n' % exc
                session.config.pluginmanager.getplugin("terminalreporter").write(
                    'WARNING: %s\n' % message, red=True, bold=True)
                warnings.warn(CovReportWarning(message))
                self.cov_total = 0
            assert self.cov_total is not None, 'Test coverage should never be `None`'
            if self._failed_cov_total():
                # make sure we get the EXIT_TESTSFAILED exit code
                compat_session.testsfailed += 1

    def pytest_terminal_summary(self, terminalreporter):
        if self._disabled:
            if self.options.no_cov_should_warn:
                message = 'Coverage disabled via --no-cov switch!'
                terminalreporter.write('WARNING: %s\n' % message, red=True, bold=True)
                warnings.warn(CovDisabledWarning(message))
            return
        if self.cov_controller is None:
            return

        if self.cov_total is None:
            # we shouldn't report, or report generation failed (error raised above)
            return

        terminalreporter.write('\n' + self.cov_report.getvalue() + '\n')

        if self.options.cov_fail_under is not None and self.options.cov_fail_under > 0:
            failed = self.cov_total < self.options.cov_fail_under
            markup = {'red': True, 'bold': True} if failed else {'green': True}
            message = (
                '{fail}Required test coverage of {required}% {reached}. '
                'Total coverage: {actual:.2f}%\n'
                .format(
                    required=self.options.cov_fail_under,
                    actual=self.cov_total,
                    fail="FAIL " if failed else "",
                    reached="not reached" if failed else "reached"
                )
            )
            terminalreporter.write(message, **markup)

    def pytest_runtest_setup(self, item):
        if os.getpid() != self.pid:
            # test is run in another process than session, run
            # coverage manually
            embed.init()

    def pytest_runtest_teardown(self, item):
        embed.cleanup()

    @compat.hookwrapper
    def pytest_runtest_call(self, item):
        if (item.get_closest_marker('no_cover')
                or 'no_cover' in getattr(item, 'fixturenames', ())):
            self.cov_controller.pause()
            yield
            self.cov_controller.resume()
        else:
            yield


class TestContextPlugin:
    def __init__(self, cov):
        self.cov = cov

    def pytest_runtest_setup(self, item):
        self.switch_context(item, 'setup')

    def pytest_runtest_teardown(self, item):
        self.switch_context(item, 'teardown')

    def pytest_runtest_call(self, item):
        self.switch_context(item, 'run')

    def switch_context(self, item, when):
        context = f"{item.nodeid}|{when}"
        self.cov.switch_context(context)
        os.environ['COV_CORE_CONTEXT'] = context


@pytest.fixture
def no_cover():
    """A pytest fixture to disable coverage."""
    pass


@pytest.fixture
def cov(request):
    """A pytest fixture to provide access to the underlying coverage object."""

    # Check with hasplugin to avoid getplugin exception in older pytest.
    if request.config.pluginmanager.hasplugin('_cov'):
        plugin = request.config.pluginmanager.getplugin('_cov')
        if plugin.cov_controller:
            return plugin.cov_controller.cov
    return None


def pytest_configure(config):
    config.addinivalue_line("markers", "no_cover: disable coverage for this test.")
