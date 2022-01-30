import sys
from contextlib import contextmanager
from typing import Iterator

if sys.version_info < (3, 6):
    from collections import OrderedDict as OrderedDict  # noqa: F401
else:
    # OrderedDict is kind of slow, so for most of our uses in Python 3.6
    # and later we'd rather just use dict
    OrderedDict = dict


if sys.version_info < (3, 7):
    @contextmanager
    def nullcontext() -> Iterator[None]:
        yield
else:
    from contextlib import nullcontext as nullcontext  # noqa: F401
