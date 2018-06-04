# -*- coding: utf-8 -*-
"""
    sphinx.quickstart
    ~~~~~~~~~~~~~~~~~

    This file has moved to :py:mod:`sphinx.cmd.quickstart`.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import warnings

from sphinx.cmd.quickstart import main as _main
from sphinx.deprecation import RemovedInSphinx20Warning


def main(*args, **kwargs):
    warnings.warn(
        '`sphinx.quickstart.main()` has moved to `sphinx.cmd.quickstart.'
        'main()`.',
        RemovedInSphinx20Warning,
        stacklevel=2,
    )
    args = args[1:]  # skip first argument to adjust arguments (refs: #4615)
    _main(*args, **kwargs)


# So program can be started with "python -m sphinx.quickstart ..."
if __name__ == "__main__":
    warnings.warn(
        '`sphinx.quickstart` has moved to `sphinx.cmd.quickstart`.',
        RemovedInSphinx20Warning,
        stacklevel=2,
    )
    main()
