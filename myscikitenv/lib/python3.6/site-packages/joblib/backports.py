"""
Backports of fixes for joblib dependencies
"""
import os
import time
import ctypes
import sys

from distutils.version import LooseVersion

try:
    import numpy as np

    def make_memmap(filename, dtype='uint8', mode='r+', offset=0,
                    shape=None, order='C'):
        """Backport of numpy memmap offset fix.

        See https://github.com/numpy/numpy/pull/8443 for more details.

        The numpy fix will be available in numpy 1.13.
        """
        mm = np.memmap(filename, dtype=dtype, mode=mode, offset=offset,
                       shape=shape, order=order)
        if LooseVersion(np.__version__) < '1.13':
            mm.offset = offset
        return mm
except ImportError:
    def make_memmap(filename, dtype='uint8', mode='r+', offset=0,
                    shape=None, order='C'):
        raise NotImplementedError(
            "'joblib.backports.make_memmap' should not be used "
            'if numpy is not installed.')


if os.name == 'nt':
    # https://github.com/joblib/joblib/issues/540
    access_denied_errors = (5, 13)
    try:
        from os import replace
    except ImportError:
        # Python 2.7
        def replace(src, dst):
            if not isinstance(src, unicode):  # noqa
                src = unicode(src, sys.getfilesystemencoding())  # noqa
            if not isinstance(dst, unicode):  # noqa
                dst = unicode(dst, sys.getfilesystemencoding())  # noqa

            movefile_replace_existing = 0x1
            return_value = ctypes.windll.kernel32.MoveFileExW(
                src, dst, movefile_replace_existing)
            if return_value == 0:
                raise ctypes.WinError()

    def concurrency_safe_rename(src, dst):
        """Renames ``src`` into ``dst`` overwriting ``dst`` if it exists.

        On Windows os.replace (or for Python 2.7 its implementation
        through MoveFileExW) can yield permission errors if executed by
        two different processes.
        """
        max_sleep_time = 1
        total_sleep_time = 0
        sleep_time = 0.001
        while total_sleep_time < max_sleep_time:
            try:
                replace(src, dst)
                break
            except Exception as exc:
                if getattr(exc, 'winerror', None) in access_denied_errors:
                    time.sleep(sleep_time)
                    total_sleep_time += sleep_time
                    sleep_time *= 2
                else:
                    raise
        else:
            raise
else:
    try:
        from os import replace as concurrency_safe_rename
    except ImportError:
        from os import rename as concurrency_safe_rename  # noqa
