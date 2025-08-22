"""Functionality with an experimental API.

.. warning::
    Although you can count on the functions in this package being
    around in the future, the API may change with any version update
    **and will not follow the skimage two-version deprecation path**.
    Therefore, use the functions herein with care, and do not use them
    in production code that will depend on updated skimage versions.
"""

import lazy_loader as _lazy

__getattr__, __dir__, __all__ = _lazy.attach_stub(__name__, __file__)
