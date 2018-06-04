# Copyright 2017 Palantir Technologies, Inc.
import functools
import inspect
import logging
import os
import threading

log = logging.getLogger(__name__)


def debounce(interval_s, keyed_by=None):
    """Debounce calls to this function until interval_s seconds have passed."""
    def wrapper(func):
        timers = {}
        lock = threading.Lock()

        @functools.wraps(func)
        def debounced(*args, **kwargs):
            call_args = inspect.getcallargs(func, *args, **kwargs)
            key = call_args[keyed_by] if keyed_by else None

            def run():
                with lock:
                    del timers[key]
                return func(*args, **kwargs)

            with lock:
                old_timer = timers.get(key)
                if old_timer:
                    old_timer.cancel()

                timer = threading.Timer(interval_s, run)
                timers[key] = timer
                timer.start()
        return debounced
    return wrapper


def find_parents(root, path, names):
    """Find files matching the given names relative to the given path.

    Args:
        path (str): The file path to start searching up from.
        names (List[str]): The file/directory names to look for.
        root (str): The directory at which to stop recursing upwards.

    Note:
        The path MUST be within the root.
    """
    if not root:
        return []

    if not os.path.commonprefix((root, path)):
        log.warning("Path %s not in %s", path, root)
        return []

    # Split the relative by directory, generate all the parent directories, then check each of them.
    # This avoids running a loop that has different base-cases for unix/windows
    # e.g. /a/b and /a/b/c/d/e.py -> ['/a/b', 'c', 'd']
    dirs = [root] + os.path.relpath(os.path.dirname(path), root).split(os.path.sep)

    # Search each of /a/b/c, /a/b, /a
    while dirs:
        search_dir = os.path.join(*dirs)
        existing = list(filter(os.path.exists, [os.path.join(search_dir, n) for n in names]))
        if existing:
            return existing
        dirs.pop()

    # Otherwise nothing
    return []


def list_to_string(value):
    return ",".join(value) if isinstance(value, list) else value


def merge_dicts(dict_a, dict_b):
    """Recursively merge dictionary b into dictionary a.

    If override_nones is True, then
    """
    def _merge_dicts_(a, b):
        for key in set(a.keys()).union(b.keys()):
            if key in a and key in b:
                if isinstance(a[key], dict) and isinstance(b[key], dict):
                    yield (key, dict(_merge_dicts_(a[key], b[key])))
                elif b[key] is not None:
                    yield (key, b[key])
                else:
                    yield (key, a[key])
            elif key in a:
                yield (key, a[key])
            elif b[key] is not None:
                yield (key, b[key])
    return dict(_merge_dicts_(dict_a, dict_b))


def format_docstring(contents):
    """Python doc strings come in a number of formats, but LSP wants markdown.

    Until we can find a fast enough way of discovering and parsing each format,
    we can do a little better by at least preserving indentation.
    """
    contents = contents.replace('\t', u'\u00A0' * 4)
    contents = contents.replace('  ', u'\u00A0' * 2)
    contents = contents.replace('*', '\\*')
    return contents


def clip_column(column, lines, line_number):
    # Normalise the position as per the LSP that accepts character positions > line length
    # https://github.com/Microsoft/language-server-protocol/blob/master/protocol.md#position
    max_column = len(lines[line_number].rstrip('\r\n')) if len(lines) > line_number else 0
    return min(column, max_column)
