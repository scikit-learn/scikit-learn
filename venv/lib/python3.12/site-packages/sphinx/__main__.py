"""The Sphinx documentation toolchain."""

from __future__ import annotations

import sys

from sphinx.cmd.build import main

raise SystemExit(main(sys.argv[1:]))
