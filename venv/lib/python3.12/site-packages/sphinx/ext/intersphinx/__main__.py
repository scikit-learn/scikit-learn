"""Command line interface for the intersphinx extension."""

from __future__ import annotations

import logging as _logging
import sys

from sphinx.ext.intersphinx._cli import inspect_main

_logging.basicConfig()

raise SystemExit(inspect_main(sys.argv[1:]))
