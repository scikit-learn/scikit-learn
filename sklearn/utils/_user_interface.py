# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import timeit
from contextlib import contextmanager


def _message_with_time(source, message, time):
    """Create one line message for logging purposes.

    Parameters
    ----------
    source : str
        String indicating the source or the reference of the message.

    message : str
        Short message.

    time : int
        Time in seconds.
    """
    start_message = "[%s] " % source

    # adapted from joblib.logger.short_format_time without the Windows -.1s
    # adjustment
    if time > 60:
        time_str = "%4.1fmin" % (time / 60)
    else:
        time_str = " %5.1fs" % time
    end_message = " %s, total=%s" % (message, time_str)
    dots_len = 70 - len(start_message) - len(end_message)
    return "%s%s%s" % (start_message, dots_len * ".", end_message)


@contextmanager
def _print_elapsed_time(source, message=None):
    """Log elapsed time to stdout when the context is exited.

    Parameters
    ----------
    source : str
        String indicating the source or the reference of the message.

    message : str, default=None
        Short message. If None, nothing will be printed.

    Returns
    -------
    context_manager
        Prints elapsed time upon exit if verbose.
    """
    if message is None:
        yield
    else:
        start = timeit.default_timer()
        yield
        print(_message_with_time(source, message, timeit.default_timer() - start))
