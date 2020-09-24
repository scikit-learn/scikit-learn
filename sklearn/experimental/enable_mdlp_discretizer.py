"""Enables MDLP discretization.

The API and results of these estimators might
change without any deprecation cycle.

Importing this file dynamically sets the
:class:`~sklearn.preprocessing.MLDPDiscretizer`
as attributes of the `preprocessing` module::

    >>> # Explicitly require this experimental feature
    >>> from sklearn.experimental import enable_mdlp_discretizer  # noqa
    >>> # Now you can import normally from preprocessing
    >>> from sklearn.preprocessing import MDLPDiscretizer

The `# noqa` comment comment can be removed: it just tells linters
like `flake8` to ignore the import, which appears as unused.
"""

from ..preprocessing._discretization import MDLPDiscretizer
from .. import preprocessing

# Use settattr to avoid mypy errors when monkeypatching
setattr(preprocessing, "MDLPDiscretizer", MDLPDiscretizer)

preprocessing.__all__ += ["MDLPDiscretizer"]
