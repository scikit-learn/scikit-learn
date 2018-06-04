from __future__ import unicode_literals
import time

__all__ = (
    'TimeIt',
)


class TimeIt(object):
    """
    Context manager that times the duration of the code body.
    The `duration` attribute will contain the execution time in seconds.
    """
    def __init__(self):
        self.duration = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.duration = self.end - self.start
