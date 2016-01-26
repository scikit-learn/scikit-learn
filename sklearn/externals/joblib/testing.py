"""
Helper for testing.
"""

import sys
import warnings
import os.path
import re
import subprocess
import threading

from sklearn.externals.joblib._compat import PY3_OR_LATER


def warnings_to_stdout():
    """ Redirect all warnings to stdout.
    """
    showwarning_orig = warnings.showwarning

    def showwarning(msg, cat, fname, lno, file=None, line=0):
        showwarning_orig(msg, cat, os.path.basename(fname), line, sys.stdout)

    warnings.showwarning = showwarning
    #warnings.simplefilter('always')


try:
    from nose.tools import assert_raises_regex
except ImportError:
    # For Python 2.7
    try:
        from nose.tools import assert_raises_regexp as assert_raises_regex
    except ImportError:
        # for Python 2.6
        def assert_raises_regex(expected_exception, expected_regexp,
                                callable_obj=None, *args, **kwargs):
            """Helper function to check for message patterns in exceptions"""

            not_raised = False
            try:
                callable_obj(*args, **kwargs)
                not_raised = True
            except Exception as e:
                error_message = str(e)
                if not re.compile(expected_regexp).search(error_message):
                    raise AssertionError("Error message should match pattern "
                                         "%r. %r does not." %
                                         (expected_regexp, error_message))
            if not_raised:
                raise AssertionError("Should have raised %r" %
                                     expected_exception(expected_regexp))


def check_subprocess_call(cmd, timeout=1, stdout_regex=None):
    """Runs a command in a subprocess with timeout in seconds.

    Also checks returncode is zero and stdout if stdout_regex is set.
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    def kill_process():
        proc.kill()

    timer = threading.Timer(timeout, kill_process)
    try:
        timer.start()
        stdout, stderr = proc.communicate()

        if PY3_OR_LATER:
            stdout, stderr = stdout.decode(), stderr.decode()
        if proc.returncode != 0:
            message = (
                'Non-zero return code: {0}.\nStdout:\n{1}\n'
                'Stderr:\n{2}').format(
                    proc.returncode, stdout, stderr)
            raise ValueError(message)

        if (stdout_regex is not None and
                not re.search(stdout_regex, stdout)):
            raise ValueError(
                "Unexpected output: '{0!r}' does not match:\n{1!r}".format(
                    stdout_regex, stdout))
    finally:
        timer.cancel()
