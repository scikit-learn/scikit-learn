"""Example images and datasets.

A curated set of general purpose and scientific images used in tests, examples,
and documentation.

Newer datasets are no longer included as part of the package, but are
downloaded on demand. To make data available offline, use :func:`download_all`.

"""

import lazy_loader as _lazy

__getattr__, __dir__, __all__ = _lazy.attach_stub(__name__, __file__)
