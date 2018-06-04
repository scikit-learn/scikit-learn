""" Access and control log capturing. """
from __future__ import absolute_import, division, print_function

import logging
from contextlib import closing, contextmanager
import re
import six

from _pytest.config import create_terminal_writer
import pytest
import py


DEFAULT_LOG_FORMAT = '%(filename)-25s %(lineno)4d %(levelname)-8s %(message)s'
DEFAULT_LOG_DATE_FORMAT = '%H:%M:%S'


class ColoredLevelFormatter(logging.Formatter):
    """
    Colorize the %(levelname)..s part of the log format passed to __init__.
    """

    LOGLEVEL_COLOROPTS = {
        logging.CRITICAL: {'red'},
        logging.ERROR: {'red', 'bold'},
        logging.WARNING: {'yellow'},
        logging.WARN: {'yellow'},
        logging.INFO: {'green'},
        logging.DEBUG: {'purple'},
        logging.NOTSET: set(),
    }
    LEVELNAME_FMT_REGEX = re.compile(r'%\(levelname\)([+-]?\d*s)')

    def __init__(self, terminalwriter, *args, **kwargs):
        super(ColoredLevelFormatter, self).__init__(
            *args, **kwargs)
        if six.PY2:
            self._original_fmt = self._fmt
        else:
            self._original_fmt = self._style._fmt
        self._level_to_fmt_mapping = {}

        levelname_fmt_match = self.LEVELNAME_FMT_REGEX.search(self._fmt)
        if not levelname_fmt_match:
            return
        levelname_fmt = levelname_fmt_match.group()

        for level, color_opts in self.LOGLEVEL_COLOROPTS.items():
            formatted_levelname = levelname_fmt % {
                'levelname': logging.getLevelName(level)}

            # add ANSI escape sequences around the formatted levelname
            color_kwargs = {name: True for name in color_opts}
            colorized_formatted_levelname = terminalwriter.markup(
                formatted_levelname, **color_kwargs)
            self._level_to_fmt_mapping[level] = self.LEVELNAME_FMT_REGEX.sub(
                colorized_formatted_levelname,
                self._fmt)

    def format(self, record):
        fmt = self._level_to_fmt_mapping.get(
            record.levelno, self._original_fmt)
        if six.PY2:
            self._fmt = fmt
        else:
            self._style._fmt = fmt
        return super(ColoredLevelFormatter, self).format(record)


def get_option_ini(config, *names):
    for name in names:
        ret = config.getoption(name)  # 'default' arg won't work as expected
        if ret is None:
            ret = config.getini(name)
        if ret:
            return ret


def pytest_addoption(parser):
    """Add options to control log capturing."""
    group = parser.getgroup('logging')

    def add_option_ini(option, dest, default=None, type=None, **kwargs):
        parser.addini(dest, default=default, type=type,
                      help='default value for ' + option)
        group.addoption(option, dest=dest, **kwargs)

    add_option_ini(
        '--no-print-logs',
        dest='log_print', action='store_const', const=False, default=True,
        type='bool',
        help='disable printing caught logs on failed tests.')
    add_option_ini(
        '--log-level',
        dest='log_level', default=None,
        help='logging level used by the logging module')
    add_option_ini(
        '--log-format',
        dest='log_format', default=DEFAULT_LOG_FORMAT,
        help='log format as used by the logging module.')
    add_option_ini(
        '--log-date-format',
        dest='log_date_format', default=DEFAULT_LOG_DATE_FORMAT,
        help='log date format as used by the logging module.')
    parser.addini(
        'log_cli', default=False, type='bool',
        help='enable log display during test run (also known as "live logging").')
    add_option_ini(
        '--log-cli-level',
        dest='log_cli_level', default=None,
        help='cli logging level.')
    add_option_ini(
        '--log-cli-format',
        dest='log_cli_format', default=None,
        help='log format as used by the logging module.')
    add_option_ini(
        '--log-cli-date-format',
        dest='log_cli_date_format', default=None,
        help='log date format as used by the logging module.')
    add_option_ini(
        '--log-file',
        dest='log_file', default=None,
        help='path to a file when logging will be written to.')
    add_option_ini(
        '--log-file-level',
        dest='log_file_level', default=None,
        help='log file logging level.')
    add_option_ini(
        '--log-file-format',
        dest='log_file_format', default=DEFAULT_LOG_FORMAT,
        help='log format as used by the logging module.')
    add_option_ini(
        '--log-file-date-format',
        dest='log_file_date_format', default=DEFAULT_LOG_DATE_FORMAT,
        help='log date format as used by the logging module.')


@contextmanager
def catching_logs(handler, formatter=None, level=None):
    """Context manager that prepares the whole logging machinery properly."""
    root_logger = logging.getLogger()

    if formatter is not None:
        handler.setFormatter(formatter)
    if level is not None:
        handler.setLevel(level)

    # Adding the same handler twice would confuse logging system.
    # Just don't do that.
    add_new_handler = handler not in root_logger.handlers

    if add_new_handler:
        root_logger.addHandler(handler)
    if level is not None:
        orig_level = root_logger.level
        root_logger.setLevel(level)
    try:
        yield handler
    finally:
        if level is not None:
            root_logger.setLevel(orig_level)
        if add_new_handler:
            root_logger.removeHandler(handler)


class LogCaptureHandler(logging.StreamHandler):
    """A logging handler that stores log records and the log text."""

    def __init__(self):
        """Creates a new log handler."""
        logging.StreamHandler.__init__(self, py.io.TextIO())
        self.records = []

    def emit(self, record):
        """Keep the log records in a list in addition to the log text."""
        self.records.append(record)
        logging.StreamHandler.emit(self, record)

    def reset(self):
        self.records = []
        self.stream = py.io.TextIO()


class LogCaptureFixture(object):
    """Provides access and control of log capturing."""

    def __init__(self, item):
        """Creates a new funcarg."""
        self._item = item
        self._initial_log_levels = {}  # type: Dict[str, int] # dict of log name -> log level

    def _finalize(self):
        """Finalizes the fixture.

        This restores the log levels changed by :meth:`set_level`.
        """
        # restore log levels
        for logger_name, level in self._initial_log_levels.items():
            logger = logging.getLogger(logger_name)
            logger.setLevel(level)

    @property
    def handler(self):
        """
        :rtype: LogCaptureHandler
        """
        return self._item.catch_log_handler

    def get_records(self, when):
        """
        Get the logging records for one of the possible test phases.

        :param str when:
            Which test phase to obtain the records from. Valid values are: "setup", "call" and "teardown".

        :rtype: List[logging.LogRecord]
        :return: the list of captured records at the given stage

        .. versionadded:: 3.4
        """
        handler = self._item.catch_log_handlers.get(when)
        if handler:
            return handler.records
        else:
            return []

    @property
    def text(self):
        """Returns the log text."""
        return self.handler.stream.getvalue()

    @property
    def records(self):
        """Returns the list of log records."""
        return self.handler.records

    @property
    def record_tuples(self):
        """Returns a list of a striped down version of log records intended
        for use in assertion comparison.

        The format of the tuple is:

            (logger_name, log_level, message)
        """
        return [(r.name, r.levelno, r.getMessage()) for r in self.records]

    def clear(self):
        """Reset the list of log records and the captured log text."""
        self.handler.reset()

    def set_level(self, level, logger=None):
        """Sets the level for capturing of logs. The level will be restored to its previous value at the end of
        the test.

        :param int level: the logger to level.
        :param str logger: the logger to update the level. If not given, the root logger level is updated.

        .. versionchanged:: 3.4
            The levels of the loggers changed by this function will be restored to their initial values at the
            end of the test.
        """
        logger_name = logger
        logger = logging.getLogger(logger_name)
        # save the original log-level to restore it during teardown
        self._initial_log_levels.setdefault(logger_name, logger.level)
        logger.setLevel(level)

    @contextmanager
    def at_level(self, level, logger=None):
        """Context manager that sets the level for capturing of logs. After the end of the 'with' statement the
        level is restored to its original value.

        :param int level: the logger to level.
        :param str logger: the logger to update the level. If not given, the root logger level is updated.
        """
        logger = logging.getLogger(logger)
        orig_level = logger.level
        logger.setLevel(level)
        try:
            yield
        finally:
            logger.setLevel(orig_level)


@pytest.fixture
def caplog(request):
    """Access and control log capturing.

    Captured logs are available through the following methods::

    * caplog.text            -> string containing formatted log output
    * caplog.records         -> list of logging.LogRecord instances
    * caplog.record_tuples   -> list of (logger_name, level, message) tuples
    * caplog.clear()         -> clear captured records and formatted log output string
    """
    result = LogCaptureFixture(request.node)
    yield result
    result._finalize()


def get_actual_log_level(config, *setting_names):
    """Return the actual logging level."""

    for setting_name in setting_names:
        log_level = config.getoption(setting_name)
        if log_level is None:
            log_level = config.getini(setting_name)
        if log_level:
            break
    else:
        return

    if isinstance(log_level, six.string_types):
        log_level = log_level.upper()
    try:
        return int(getattr(logging, log_level, log_level))
    except ValueError:
        # Python logging does not recognise this as a logging level
        raise pytest.UsageError(
            "'{0}' is not recognized as a logging level name for "
            "'{1}'. Please consider passing the "
            "logging level num instead.".format(
                log_level,
                setting_name))


def pytest_configure(config):
    config.pluginmanager.register(LoggingPlugin(config), 'logging-plugin')


@contextmanager
def _dummy_context_manager():
    yield


class LoggingPlugin(object):
    """Attaches to the logging module and captures log messages for each test.
    """

    def __init__(self, config):
        """Creates a new plugin to capture log messages.

        The formatter can be safely shared across all handlers so
        create a single one for the entire test session here.
        """
        self._config = config

        # enable verbose output automatically if live logging is enabled
        if self._log_cli_enabled() and not config.getoption('verbose'):
            # sanity check: terminal reporter should not have been loaded at this point
            assert self._config.pluginmanager.get_plugin('terminalreporter') is None
            config.option.verbose = 1

        self.print_logs = get_option_ini(config, 'log_print')
        self.formatter = logging.Formatter(get_option_ini(config, 'log_format'),
                                           get_option_ini(config, 'log_date_format'))
        self.log_level = get_actual_log_level(config, 'log_level')

        log_file = get_option_ini(config, 'log_file')
        if log_file:
            self.log_file_level = get_actual_log_level(config, 'log_file_level')

            log_file_format = get_option_ini(config, 'log_file_format', 'log_format')
            log_file_date_format = get_option_ini(config, 'log_file_date_format', 'log_date_format')
            # Each pytest runtests session will write to a clean logfile
            self.log_file_handler = logging.FileHandler(log_file, mode='w')
            log_file_formatter = logging.Formatter(log_file_format, datefmt=log_file_date_format)
            self.log_file_handler.setFormatter(log_file_formatter)
        else:
            self.log_file_handler = None

        # initialized during pytest_runtestloop
        self.log_cli_handler = None

    def _log_cli_enabled(self):
        """Return True if log_cli should be considered enabled, either explicitly
        or because --log-cli-level was given in the command-line.
        """
        return self._config.getoption('--log-cli-level') is not None or \
            self._config.getini('log_cli')

    @contextmanager
    def _runtest_for(self, item, when):
        """Implements the internals of pytest_runtest_xxx() hook."""
        with catching_logs(LogCaptureHandler(),
                           formatter=self.formatter, level=self.log_level) as log_handler:
            if self.log_cli_handler:
                self.log_cli_handler.set_when(when)

            if item is None:
                yield  # run the test
                return

            if not hasattr(item, 'catch_log_handlers'):
                item.catch_log_handlers = {}
            item.catch_log_handlers[when] = log_handler
            item.catch_log_handler = log_handler
            try:
                yield  # run test
            finally:
                del item.catch_log_handler
                if when == 'teardown':
                    del item.catch_log_handlers

            if self.print_logs:
                # Add a captured log section to the report.
                log = log_handler.stream.getvalue().strip()
                item.add_report_section(when, 'log', log)

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_setup(self, item):
        with self._runtest_for(item, 'setup'):
            yield

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_call(self, item):
        with self._runtest_for(item, 'call'):
            yield

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_teardown(self, item):
        with self._runtest_for(item, 'teardown'):
            yield

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_logstart(self):
        if self.log_cli_handler:
            self.log_cli_handler.reset()
        with self._runtest_for(None, 'start'):
            yield

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_logfinish(self):
        with self._runtest_for(None, 'finish'):
            yield

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtestloop(self, session):
        """Runs all collected test items."""
        self._setup_cli_logging()
        with self.live_logs_context:
            if self.log_file_handler is not None:
                with closing(self.log_file_handler):
                    with catching_logs(self.log_file_handler,
                                       level=self.log_file_level):
                        yield  # run all the tests
            else:
                yield  # run all the tests

    def _setup_cli_logging(self):
        """Sets up the handler and logger for the Live Logs feature, if enabled.

        This must be done right before starting the loop so we can access the terminal reporter plugin.
        """
        terminal_reporter = self._config.pluginmanager.get_plugin('terminalreporter')
        if self._log_cli_enabled() and terminal_reporter is not None:
            capture_manager = self._config.pluginmanager.get_plugin('capturemanager')
            log_cli_handler = _LiveLoggingStreamHandler(terminal_reporter, capture_manager)
            log_cli_format = get_option_ini(self._config, 'log_cli_format', 'log_format')
            log_cli_date_format = get_option_ini(self._config, 'log_cli_date_format', 'log_date_format')
            if self._config.option.color != 'no' and ColoredLevelFormatter.LEVELNAME_FMT_REGEX.search(log_cli_format):
                log_cli_formatter = ColoredLevelFormatter(create_terminal_writer(self._config),
                                                          log_cli_format, datefmt=log_cli_date_format)
            else:
                log_cli_formatter = logging.Formatter(log_cli_format, datefmt=log_cli_date_format)
            log_cli_level = get_actual_log_level(self._config, 'log_cli_level', 'log_level')
            self.log_cli_handler = log_cli_handler
            self.live_logs_context = catching_logs(log_cli_handler, formatter=log_cli_formatter, level=log_cli_level)
        else:
            self.live_logs_context = _dummy_context_manager()


class _LiveLoggingStreamHandler(logging.StreamHandler):
    """
    Custom StreamHandler used by the live logging feature: it will write a newline before the first log message
    in each test.

    During live logging we must also explicitly disable stdout/stderr capturing otherwise it will get captured
    and won't appear in the terminal.
    """

    def __init__(self, terminal_reporter, capture_manager):
        """
        :param _pytest.terminal.TerminalReporter terminal_reporter:
        :param _pytest.capture.CaptureManager capture_manager:
        """
        logging.StreamHandler.__init__(self, stream=terminal_reporter)
        self.capture_manager = capture_manager
        self.reset()
        self.set_when(None)
        self._test_outcome_written = False

    def reset(self):
        """Reset the handler; should be called before the start of each test"""
        self._first_record_emitted = False

    def set_when(self, when):
        """Prepares for the given test phase (setup/call/teardown)"""
        self._when = when
        self._section_name_shown = False
        if when == 'start':
            self._test_outcome_written = False

    def emit(self, record):
        if self.capture_manager is not None:
            self.capture_manager.suspend_global_capture()
        try:
            if not self._first_record_emitted:
                self.stream.write('\n')
                self._first_record_emitted = True
            elif self._when in ('teardown', 'finish'):
                if not self._test_outcome_written:
                    self._test_outcome_written = True
                    self.stream.write('\n')
            if not self._section_name_shown and self._when:
                self.stream.section('live log ' + self._when, sep='-', bold=True)
                self._section_name_shown = True
            logging.StreamHandler.emit(self, record)
        finally:
            if self.capture_manager is not None:
                self.capture_manager.resume_global_capture()
