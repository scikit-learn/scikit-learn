# encoding: utf-8
"""
Utilities for working with terminals.

Authors:

* Brian E. Granger
* Fernando Perez
* Alexander Belchenko (e-mail: bialix AT ukr.net)
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import sys
import warnings
from shutil import get_terminal_size as _get_terminal_size

# This variable is part of the expected API of the module:
ignore_termtitle = True



if os.name == 'posix':
    def _term_clear():
        os.system('clear')
elif sys.platform == 'win32':
    def _term_clear():
        os.system('cls')
else:
    def _term_clear():
        pass



def toggle_set_term_title(val):
    """Control whether set_term_title is active or not.

    set_term_title() allows writing to the console titlebar.  In embedded
    widgets this can cause problems, so this call can be used to toggle it on
    or off as needed.

    The default state of the module is for the function to be disabled.

    Parameters
    ----------
      val : bool
        If True, set_term_title() actually writes to the terminal (using the
        appropriate platform-specific module).  If False, it is a no-op.
    """
    global ignore_termtitle
    ignore_termtitle = not(val)


def _set_term_title(*args,**kw):
    """Dummy no-op."""
    pass


def _set_term_title_xterm(title):
    """ Change virtual terminal title in xterm-workalikes """
    sys.stdout.write('\033]0;%s\007' % title)

if os.name == 'posix':
    TERM = os.environ.get('TERM','')
    if TERM.startswith('xterm'):
        _set_term_title = _set_term_title_xterm
elif sys.platform == 'win32':
    try:
        import ctypes

        SetConsoleTitleW = ctypes.windll.kernel32.SetConsoleTitleW
        SetConsoleTitleW.argtypes = [ctypes.c_wchar_p]
    
        def _set_term_title(title):
            """Set terminal title using ctypes to access the Win32 APIs."""
            SetConsoleTitleW(title)
    except ImportError:
        def _set_term_title(title):
            """Set terminal title using the 'title' command."""
            global ignore_termtitle

            try:
                # Cannot be on network share when issuing system commands
                curr = os.getcwd()
                os.chdir("C:")
                ret = os.system("title " + title)
            finally:
                os.chdir(curr)
            if ret:
                # non-zero return code signals error, don't try again
                ignore_termtitle = True


def set_term_title(title):
    """Set terminal title using the necessary platform-dependent calls."""
    if ignore_termtitle:
        return
    _set_term_title(title)


def freeze_term_title():
    warnings.warn("This function is deprecated, use toggle_set_term_title()")
    global ignore_termtitle
    ignore_termtitle = True


def get_terminal_size(defaultx=80, defaulty=25):
    return _get_terminal_size((defaultx, defaulty))
