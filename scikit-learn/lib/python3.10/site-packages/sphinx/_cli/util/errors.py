from __future__ import annotations

import re
import sys
import tempfile
from typing import TYPE_CHECKING, TextIO

from sphinx.errors import SphinxParallelError

if TYPE_CHECKING:
    from sphinx.application import Sphinx

_ANSI_COLOUR_CODES: re.Pattern[str] = re.compile('\x1b.*?m')


def terminal_safe(s: str, /) -> str:
    """Safely encode a string for printing to the terminal."""
    return s.encode('ascii', 'backslashreplace').decode('ascii')


def strip_colors(s: str, /) -> str:
    return _ANSI_COLOUR_CODES.sub('', s).strip()


def error_info(messages: str, extensions: str, traceback: str) -> str:
    import platform

    import docutils
    import jinja2
    import pygments

    import sphinx

    return f"""\
Versions
========

* Platform:         {sys.platform}; ({platform.platform()})
* Python version:   {platform.python_version()} ({platform.python_implementation()})
* Sphinx version:   {sphinx.__display_version__}
* Docutils version: {docutils.__version__}
* Jinja2 version:   {jinja2.__version__}
* Pygments version: {pygments.__version__}

Last Messages
=============

{messages}

Loaded Extensions
=================

{extensions}

Traceback
=========

{traceback}
"""


def save_traceback(app: Sphinx | None, exc: BaseException) -> str:
    """Save the given exception's traceback in a temporary file."""
    if isinstance(exc, SphinxParallelError):
        exc_format = '(Error in parallel process)\n' + exc.traceback
    else:
        import traceback

        exc_format = traceback.format_exc()

    last_msgs = exts_list = ''
    if app is not None:
        extensions = app.extensions.values()
        last_msgs = '\n'.join(f'* {strip_colors(s)}' for s in app.messagelog)
        exts_list = '\n'.join(
            f'* {ext.name} ({ext.version})'
            for ext in extensions
            if ext.version != 'builtin'
        )

    with tempfile.NamedTemporaryFile(
        suffix='.log', prefix='sphinx-err-', delete=False
    ) as f:
        f.write(error_info(last_msgs, exts_list, exc_format).encode('utf-8'))

    return f.name


def handle_exception(
    exception: BaseException,
    /,
    *,
    stderr: TextIO = sys.stderr,
    use_pdb: bool = False,
    print_traceback: bool = False,
    app: Sphinx | None = None,
) -> None:
    from bdb import BdbQuit
    from traceback import TracebackException, print_exc

    from docutils.utils import SystemMessage

    from sphinx._cli.util.colour import red
    from sphinx.errors import SphinxError
    from sphinx.locale import __

    if isinstance(exception, BdbQuit):
        return

    def print_err(*values: str) -> None:
        print(*values, file=stderr)

    def print_red(*values: str) -> None:
        print_err(*map(red, values))

    print_err()
    if print_traceback or use_pdb:
        print_exc(file=stderr)
        print_err()

    if use_pdb:
        from pdb import post_mortem

        print_red(__('Exception occurred, starting debugger:'))
        post_mortem()
        return

    if isinstance(exception, KeyboardInterrupt):
        print_err(__('Interrupted!'))
        return

    if isinstance(exception, SystemMessage):
        print_red(__('reStructuredText markup error:'))
        print_err(str(exception))
        return

    if isinstance(exception, SphinxError):
        print_red(f'{exception.category}:')
        print_err(str(exception))
        return

    if isinstance(exception, UnicodeError):
        print_red(__('Encoding error:'))
        print_err(str(exception))
        return

    if isinstance(exception, RecursionError):
        print_red(__('Recursion error:'))
        print_err(str(exception))
        print_err()
        print_err(
            __(
                'This can happen with very large or deeply nested source '
                'files. You can carefully increase the default Python '
                'recursion limit of 1000 in conf.py with e.g.:'
            )
        )
        print_err('\n    import sys\n    sys.setrecursionlimit(1_500)\n')
        return

    # format an exception with traceback, but only the last frame.
    te = TracebackException.from_exception(exception, limit=-1)
    formatted_tb = te.stack.format()[-1] + ''.join(te.format_exception_only()).rstrip()

    print_red(__('Exception occurred:'))
    print_err(formatted_tb)
    traceback_info_path = save_traceback(app, exception)
    print_err(__('The full traceback has been saved in:'))
    print_err(traceback_info_path)
    print_err()
    print_err(
        __(
            'To report this error to the developers, please open an issue '
            'at <https://github.com/sphinx-doc/sphinx/issues/>. Thanks!'
        )
    )
    print_err(
        __(
            'Please also report this if it was a user error, so '
            'that a better error message can be provided next time.'
        )
    )
