# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2014 The Meson development team

"""This is (mostly) a standalone module used to write logging
information about Meson runs. Some output goes to screen,
some to logging dir and some goes to both."""

from __future__ import annotations

import enum
import os
import io
import sys
import time
import platform
import shlex
import subprocess
import shutil
import typing as T
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

if T.TYPE_CHECKING:
    from typing_extensions import Literal

    from ._typing import StringProtocol, SizedStringProtocol
    from .mparser import BaseNode

    TV_Loggable = T.Union[str, 'AnsiDecorator', StringProtocol]
    TV_LoggableList = T.List[TV_Loggable]

def is_windows() -> bool:
    platname = platform.system().lower()
    return platname == 'windows'

def _windows_ansi() -> bool:
    # windll only exists on windows, so mypy will get mad
    from ctypes import windll, byref  # type: ignore
    from ctypes.wintypes import DWORD

    kernel = windll.kernel32
    stdout = kernel.GetStdHandle(-11)
    mode = DWORD()
    if not kernel.GetConsoleMode(stdout, byref(mode)):
        return False
    # ENABLE_VIRTUAL_TERMINAL_PROCESSING == 0x4
    # If the call to enable VT processing fails (returns 0), we fallback to
    # original behavior
    return bool(kernel.SetConsoleMode(stdout, mode.value | 0x4) or os.environ.get('ANSICON'))

_in_ci = 'CI' in os.environ
_ci_is_github = 'GITHUB_ACTIONS' in os.environ


class _Severity(enum.Enum):

    NOTICE = enum.auto()
    WARNING = enum.auto()
    ERROR = enum.auto()
    DEPRECATION = enum.auto()

@dataclass
class _Logger:

    log_dir: T.Optional[str] = None
    log_depth: T.List[str] = field(default_factory=list)
    log_to_stderr: bool = False
    log_file: T.Optional[T.TextIO] = None
    log_timestamp_start: T.Optional[float] = None
    log_fatal_warnings = False
    log_disable_stdout = False
    log_errors_only = False
    logged_once: T.Set[T.Tuple[str, ...]] = field(default_factory=set)
    log_warnings_counter = 0
    log_pager: T.Optional['subprocess.Popen'] = None

    _LOG_FNAME: T.ClassVar[str] = 'meson-log.txt'

    @contextmanager
    def no_logging(self) -> T.Iterator[None]:
        self.log_disable_stdout = True
        try:
            yield
        finally:
            self.log_disable_stdout = False

    @contextmanager
    def force_logging(self) -> T.Iterator[None]:
        restore = self.log_disable_stdout
        self.log_disable_stdout = False
        try:
            yield
        finally:
            self.log_disable_stdout = restore

    def set_quiet(self) -> None:
        self.log_errors_only = True

    def set_verbose(self) -> None:
        self.log_errors_only = False

    def set_timestamp_start(self, start: float) -> None:
        self.log_timestamp_start = start

    def shutdown(self) -> T.Optional[str]:
        if self.log_file is not None:
            path = self.log_file.name
            exception_around_goer = self.log_file
            self.log_file = None
            exception_around_goer.close()
            return path
        self.stop_pager()
        return None

    def start_pager(self) -> None:
        if not self.colorize_console():
            return
        pager_cmd = []
        if 'PAGER' in os.environ:
            pager_cmd = shlex.split(os.environ['PAGER'])
        else:
            less = shutil.which('less')
            if not less and is_windows():
                git = shutil.which('git')
                if git:
                    path = Path(git).parents[1] / 'usr' / 'bin'
                    less = shutil.which('less', path=str(path))
            if less:
                pager_cmd = [less]
        if not pager_cmd:
            return
        try:
            # Set 'LESS' environment variable, rather than arguments in
            # pager_cmd, to also support the case where the user has 'PAGER'
            # set to 'less'. Arguments set are:
            # "R" : support color
            # "X" : do not clear the screen when leaving the pager
            # "F" : skip the pager if content fits into the screen
            env = os.environ.copy()
            if 'LESS' not in env:
                env['LESS'] = 'RXF'
            # Set "-c" for lv to support color
            if 'LV' not in env:
                env['LV'] = '-c'
            self.log_pager = subprocess.Popen(pager_cmd, stdin=subprocess.PIPE,
                                              text=True, encoding='utf-8', env=env)
        except Exception as e:
            # Ignore errors, unless it is a user defined pager.
            if 'PAGER' in os.environ:
                from .mesonlib import MesonException
                raise MesonException(f'Failed to start pager: {str(e)}')

    def stop_pager(self) -> None:
        if self.log_pager:
            try:
                self.log_pager.stdin.flush()
                self.log_pager.stdin.close()
            except OSError:
                pass
            self.log_pager.wait()
            self.log_pager = None

    def initialize(self, logdir: str, fatal_warnings: bool = False) -> None:
        self.log_dir = logdir
        self.log_file = open(os.path.join(logdir, self._LOG_FNAME), 'w', encoding='utf-8')
        self.log_fatal_warnings = fatal_warnings

    def process_markup(self, args: T.Sequence[TV_Loggable], keep: bool, display_timestamp: bool = True) -> T.List[str]:
        arr: T.List[str] = []
        if self.log_timestamp_start is not None and display_timestamp:
            arr = ['[{:.3f}]'.format(time.monotonic() - self.log_timestamp_start)]
        for arg in args:
            if arg is None:
                continue
            if isinstance(arg, str):
                arr.append(arg)
            elif isinstance(arg, AnsiDecorator):
                arr.append(arg.get_text(keep))
            else:
                arr.append(str(arg))
        return arr

    def force_print(self, *args: str, nested: bool, sep: T.Optional[str] = None,
                    end: T.Optional[str] = None) -> None:
        if self.log_disable_stdout:
            return
        iostr = io.StringIO()
        print(*args, sep=sep, end=end, file=iostr)

        raw = iostr.getvalue()
        if self.log_depth:
            prepend = self.log_depth[-1] + '| ' if nested else ''
            lines = []
            for l in raw.split('\n'):
                l = l.strip()
                lines.append(prepend + l if l else '')
            raw = '\n'.join(lines)

        # _Something_ is going to get printed.
        if self.log_pager:
            output = self.log_pager.stdin
        elif self.log_to_stderr:
            output = sys.stderr
        else:
            output = sys.stdout
        try:
            print(raw, end='', file=output)
        except UnicodeEncodeError:
            cleaned = raw.encode('ascii', 'replace').decode('ascii')
            print(cleaned, end='', file=output)

    def debug(self, *args: TV_Loggable, sep: T.Optional[str] = None,
              end: T.Optional[str] = None, display_timestamp: bool = True) -> None:
        arr = process_markup(args, False, display_timestamp)
        if self.log_file is not None:
            print(*arr, file=self.log_file, sep=sep, end=end)
            self.log_file.flush()

    def _log(self, *args: TV_Loggable, is_error: bool = False,
             nested: bool = True, sep: T.Optional[str] = None,
             end: T.Optional[str] = None, display_timestamp: bool = True) -> None:
        arr = process_markup(args, False, display_timestamp)
        if self.log_file is not None:
            print(*arr, file=self.log_file, sep=sep, end=end)
            self.log_file.flush()
        if self.colorize_console():
            arr = process_markup(args, True, display_timestamp)
        if not self.log_errors_only or is_error:
            force_print(*arr, nested=nested, sep=sep, end=end)

    def _debug_log_cmd(self, cmd: str, args: T.List[str]) -> None:
        if not _in_ci:
            return
        args = [f'"{x}"' for x in args]  # Quote all args, just in case
        self.debug('!meson_ci!/{} {}'.format(cmd, ' '.join(args)))

    def cmd_ci_include(self, file: str) -> None:
        self._debug_log_cmd('ci_include', [file])

    def log(self, *args: TV_Loggable, is_error: bool = False,
            once: bool = False, nested: bool = True,
            sep: T.Optional[str] = None,
            end: T.Optional[str] = None,
            display_timestamp: bool = True) -> None:
        if self._should_log(*args, once=once):
            self._log(*args, is_error=is_error, nested=nested, sep=sep, end=end, display_timestamp=display_timestamp)

    def log_timestamp(self, *args: TV_Loggable) -> None:
        if self.log_timestamp_start:
            self.log(*args)

    def _should_log(self, *args: TV_Loggable, once: bool) -> bool:
        def to_str(x: TV_Loggable) -> str:
            if isinstance(x, str):
                return x
            if isinstance(x, AnsiDecorator):
                return x.text
            return str(x)
        if not once:
            return True
        t = tuple(to_str(a) for a in args)
        if t in self.logged_once:
            return False
        self.logged_once.add(t)
        return True

    def _log_error(self, severity: _Severity, *rargs: TV_Loggable,
                   once: bool = False, fatal: bool = True,
                   location: T.Optional[BaseNode] = None,
                   nested: bool = True, sep: T.Optional[str] = None,
                   end: T.Optional[str] = None,
                   is_error: bool = True) -> None:
        from .mesonlib import MesonException, relpath

        # The typing requirements here are non-obvious. Lists are invariant,
        # therefore T.List[A] and T.List[T.Union[A, B]] are not able to be joined
        if severity is _Severity.NOTICE:
            label: TV_LoggableList = [bold('NOTICE:')]
        elif severity is _Severity.WARNING:
            label = [yellow('WARNING:')]
        elif severity is _Severity.ERROR:
            label = [red('ERROR:')]
        elif severity is _Severity.DEPRECATION:
            label = [red('DEPRECATION:')]
        # rargs is a tuple, not a list
        args = label + list(rargs)

        if not self._should_log(*args, once=once):
            return

        if location is not None:
            location_file = relpath(location.filename, os.getcwd())
            location_str = get_error_location_string(location_file, location.lineno)
            # Unions are frankly awful, and we have to T.cast here to get mypy
            # to understand that the list concatenation is safe
            location_list = T.cast('TV_LoggableList', [location_str])
            args = location_list + args

        self._log(*args, nested=nested, sep=sep, end=end, is_error=is_error)

        self.log_warnings_counter += 1

        if self.log_fatal_warnings and fatal:
            raise MesonException("Fatal warnings enabled, aborting")

    def error(self, *args: TV_Loggable,
              once: bool = False, fatal: bool = True,
              location: T.Optional[BaseNode] = None,
              nested: bool = True, sep: T.Optional[str] = None,
              end: T.Optional[str] = None) -> None:
        return self._log_error(_Severity.ERROR, *args, once=once, fatal=fatal, location=location,
                               nested=nested, sep=sep, end=end, is_error=True)

    def warning(self, *args: TV_Loggable,
                once: bool = False, fatal: bool = True,
                location: T.Optional[BaseNode] = None,
                nested: bool = True, sep: T.Optional[str] = None,
                end: T.Optional[str] = None) -> None:
        return self._log_error(_Severity.WARNING, *args, once=once, fatal=fatal, location=location,
                               nested=nested, sep=sep, end=end, is_error=True)

    def deprecation(self, *args: TV_Loggable,
                    once: bool = False, fatal: bool = True,
                    location: T.Optional[BaseNode] = None,
                    nested: bool = True, sep: T.Optional[str] = None,
                    end: T.Optional[str] = None) -> None:
        return self._log_error(_Severity.DEPRECATION, *args, once=once, fatal=fatal, location=location,
                               nested=nested, sep=sep, end=end, is_error=True)

    def notice(self, *args: TV_Loggable,
               once: bool = False, fatal: bool = True,
               location: T.Optional[BaseNode] = None,
               nested: bool = True, sep: T.Optional[str] = None,
               end: T.Optional[str] = None) -> None:
        return self._log_error(_Severity.NOTICE, *args, once=once, fatal=fatal, location=location,
                               nested=nested, sep=sep, end=end, is_error=False)

    def exception(self, e: Exception, prefix: T.Optional[AnsiDecorator] = None) -> None:
        if prefix is None:
            prefix = red('ERROR:')
        self.log()
        args: T.List[T.Union[AnsiDecorator, str]] = []
        if all(getattr(e, a, None) is not None for a in ['file', 'lineno', 'colno']):
            # Mypy doesn't follow hasattr, and it's pretty easy to visually inspect
            # that this is correct, so we'll just ignore it.
            path = get_relative_path(Path(e.file), Path(os.getcwd()))  # type: ignore
            args.append(f'{path}:{e.lineno}:{e.colno}:')  # type: ignore
        if prefix:
            args.append(prefix)
        args.append(str(e))

        with self.force_logging():
            self.log(*args, is_error=True)

    @contextmanager
    def nested(self, name: str = '') -> T.Generator[None, None, None]:
        self.log_depth.append(name)
        try:
            yield
        finally:
            self.log_depth.pop()

    def get_log_dir(self) -> str:
        return self.log_dir

    def get_log_depth(self) -> int:
        return len(self.log_depth)

    @contextmanager
    def nested_warnings(self) -> T.Iterator[None]:
        old = self.log_warnings_counter
        self.log_warnings_counter = 0
        try:
            yield
        finally:
            self.log_warnings_counter = old

    def get_warning_count(self) -> int:
        return self.log_warnings_counter

    def redirect(self, to_stderr: bool) -> None:
        self.log_to_stderr = to_stderr

    def colorize_console(self) -> bool:
        output = sys.stderr if self.log_to_stderr else sys.stdout
        _colorize_console: bool = getattr(output, 'colorize_console', None)
        if _colorize_console is not None:
            return _colorize_console
        try:
            if is_windows():
                _colorize_console = os.isatty(output.fileno()) and _windows_ansi()
            else:
                _colorize_console = os.isatty(output.fileno()) and os.environ.get('TERM', 'dumb') != 'dumb'
        except Exception:
            _colorize_console = False
        output.colorize_console = _colorize_console  # type: ignore[attr-defined]
        return _colorize_console

    def setup_console(self) -> None:
        # on Windows, a subprocess might call SetConsoleMode() on the console
        # connected to stdout and turn off ANSI escape processing. Call this after
        # running a subprocess to ensure we turn it on again.
        output = sys.stderr if self.log_to_stderr else sys.stdout
        if is_windows():
            try:
                delattr(output, 'colorize_console')
            except AttributeError:
                pass

_logger = _Logger()
cmd_ci_include = _logger.cmd_ci_include
colorize_console = _logger.colorize_console
debug = _logger.debug
deprecation = _logger.deprecation
error = _logger.error
exception = _logger.exception
force_print = _logger.force_print
get_log_depth = _logger.get_log_depth
get_log_dir = _logger.get_log_dir
get_warning_count = _logger.get_warning_count
initialize = _logger.initialize
log = _logger.log
log_timestamp = _logger.log_timestamp
nested = _logger.nested
nested_warnings = _logger.nested_warnings
no_logging = _logger.no_logging
notice = _logger.notice
process_markup = _logger.process_markup
redirect = _logger.redirect
set_quiet = _logger.set_quiet
set_timestamp_start = _logger.set_timestamp_start
set_verbose = _logger.set_verbose
setup_console = _logger.setup_console
shutdown = _logger.shutdown
start_pager = _logger.start_pager
stop_pager = _logger.stop_pager
warning = _logger.warning

class AnsiDecorator:
    plain_code = "\033[0m"

    def __init__(self, text: str, code: str, quoted: bool = False):
        self.text = text
        self.code = code
        self.quoted = quoted

    def get_text(self, with_codes: bool) -> str:
        text = self.text
        if with_codes and self.code:
            text = self.code + self.text + AnsiDecorator.plain_code
        if self.quoted:
            text = f'"{text}"'
        return text

    def __len__(self) -> int:
        return len(self.text)

    def __str__(self) -> str:
        return self.get_text(colorize_console())

class AnsiText:
    def __init__(self, *args: 'SizedStringProtocol'):
        self.args = args

    def __len__(self) -> int:
        return sum(len(x) for x in self.args)

    def __str__(self) -> str:
        return ''.join(str(x) for x in self.args)


def bold(text: str, quoted: bool = False) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1m", quoted=quoted)

def italic(text: str, quoted: bool = False) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[3m", quoted=quoted)

def plain(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "")

def red(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1;31m")

def green(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1;32m")

def yellow(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1;33m")

def blue(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1;34m")

def cyan(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[1;36m")

def normal_red(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[31m")

def normal_green(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[32m")

def normal_yellow(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[33m")

def normal_blue(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[34m")

def normal_cyan(text: str) -> AnsiDecorator:
    return AnsiDecorator(text, "\033[36m")

def get_error_location_string(fname: StringProtocol, lineno: int) -> str:
    return f'{fname}:{lineno}:'

def get_relative_path(target: Path, current: Path) -> Path:
    """Get the path to target from current"""
    # Go up "current" until we find a common ancestor to target
    acc = ['.']
    for part in [current, *current.parents]:
        try:
            path = target.relative_to(part)
            return Path(*acc, path)
        except ValueError:
            pass
        acc += ['..']

    # we failed, should not get here
    return target

# Format a list for logging purposes as a string. It separates
# all but the last item with commas, and the last with 'and'.
def format_list(input_list: T.List[str]) -> str:
    l = len(input_list)
    if l > 2:
        return ' and '.join([', '.join(input_list[:-1]), input_list[-1]])
    elif l == 2:
        return ' and '.join(input_list)
    elif l == 1:
        return input_list[0]
    else:
        return ''


def code_line(text: str, line: str, colno: int) -> str:
    """Print a line with a caret pointing to the colno

    :param text: A message to display before the line
    :param line: The line of code to be pointed to
    :param colno: The column number to point at
    :return: A formatted string of the text, line, and a caret
    """
    return f'{text}\n{line}\n{" " * colno}^'

@T.overload
def ci_fold_file(fname: T.Union[str, os.PathLike], banner: str, force: Literal[True] = True) -> str: ...

@T.overload
def ci_fold_file(fname: T.Union[str, os.PathLike], banner: str, force: Literal[False] = False) -> T.Optional[str]: ...

def ci_fold_file(fname: T.Union[str, os.PathLike], banner: str, force: bool = False) -> T.Optional[str]:
    if not _in_ci and not force:
        return None

    if _ci_is_github:
        header = f'::group::==== {banner} ===='
        footer = '::endgroup::'
    elif force:
        header = banner
        footer = ''
    elif 'MESON_FORCE_SHOW_LOGS' in os.environ:
        header = f'==== Forcing display of logs for {os.path.basename(fname)} ===='
        footer = ''
    else:
        # only github is implemented
        return None

    with open(fname, 'r', encoding='utf-8') as f:
        data = f.read()
    return f'{header}\n{data}\n{footer}\n'
