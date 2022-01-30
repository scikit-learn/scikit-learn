# encoding: utf-8
"""
The *pathspec* package provides pattern matching for file paths. So far
this only includes Git's wildmatch pattern matching (the style used for
".gitignore" files).

The following classes are imported and made available from the root of
the `pathspec` package:

- :class:`pathspec.pathspec.PathSpec`

- :class:`pathspec.pattern.Pattern`

- :class:`pathspec.pattern.RegexPattern`

- :class:`pathspec.util.RecursionError`

The following functions are also imported:

- :func:`pathspec.util.iter_tree`
- :func:`pathspec.util.lookup_pattern`
- :func:`pathspec.util.match_files`
"""
from __future__ import unicode_literals

from .pathspec import PathSpec
from .pattern import Pattern, RegexPattern
from .util import iter_tree, lookup_pattern, match_files, RecursionError

from ._meta import (
	__author__,
	__copyright__,
	__credits__,
	__license__,
	__version__,
)

# Load pattern implementations.
from . import patterns

# Expose `GitIgnorePattern` class in the root module for backward
# compatibility with v0.4.
from .patterns.gitwildmatch import GitIgnorePattern
