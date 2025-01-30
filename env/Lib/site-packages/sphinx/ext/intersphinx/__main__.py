"""Command line interface for the intersphinx extension."""

import logging as _logging
import sys

from sphinx.ext.intersphinx import inspect_main

_logging.basicConfig()

raise SystemExit(inspect_main(sys.argv[1:]))
