"""Thread-safe implementation for warnings.catch_warnings

The default implementation of warnings.catch_warnings is not thread-safe [1]
which can cause issues in scikit-learn as we sometime need to customize the
warnings behavior locally in the library rendering it thread-unsafe.

[1] http://docs.python.org/2/library/warnings.html#available-context-managers

"""
# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import warnings
from threading import RLock


_lock = RLock()


class catch_warnings(warnings.catch_warnings):
    """Thread-safe wrapper for warnings.catch_warnings

    If record_category is a warning class, all instances of that category
    raised inside the context manager block are guaranteed to be recorded. The
    filter configuration is restaured when exiting the block.

    This implementation uses a threading.RLock which guarantees thread-safety
    but can harm the concurrency performance of multi-threaded applications.

    sklearn code should therefore avoid to use catch_warnings in multi-threaded
    performance critical code sections.

    """

    def __init__(self, record=False, module=None, record_category=None):
        if record_category is not None:
            record = True
        super(catch_warnings, self).__init__(record=record, module=module)
        self.record_category = record_category

    def __enter__(self):
        _lock.acquire()
        if self.record_category is not None:
            self._module.simplefilter('always', self.record_category)
        return super(catch_warnings, self).__enter__()

    def __exit__(self, *exc_info):
        out = super(catch_warnings, self).__exit__(*exc_info)
        if self.record_category is not None:
            self._module.filters.pop(0)
        _lock.release()
        return out
