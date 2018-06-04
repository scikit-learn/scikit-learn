# -*- coding: utf-8 -*-
"""
    sphinx.cmd.build
    ~~~~~~~~~~~~~~~~

    Build documentation from a provided source.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import sys

if False:
    # For type annotation
    from typing import List  # NOQA


def build_main(argv=sys.argv[1:]):
    # type: (List[str]) -> int
    """Sphinx build "main" command-line entry."""
    from sphinx import cmdline
    return cmdline.main(argv)  # type: ignore


def make_main(argv=sys.argv[1:]):
    # type: (List[str]) -> int
    """Sphinx build "make mode" entry."""
    from sphinx import make_mode
    return make_mode.run_make_mode(argv[1:])  # type: ignore


def main(argv=sys.argv[1:]):
    # type: (List[str]) -> int
    if sys.argv[1:2] == ['-M']:
        return make_main(argv)
    else:
        return build_main(argv)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
