"""General shared utilities."""
import logging
from itertools import tee
try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest


# Do not update the version manually - it is managed by `bumpversion`.
__version__ = '2.1.1'
log = logging.getLogger(__name__)


def is_blank(string):
    """Return True iff the string contains only whitespaces."""
    return not string.strip()


def pairwise(iterable, default_value):
    """Return pairs of items from `iterable`.

    pairwise([1, 2, 3], default_value=None) -> (1, 2) (2, 3), (3, None)
    """
    a, b = tee(iterable)
    _ = next(b, default_value)
    return zip_longest(a, b, fillvalue=default_value)
