# encoding: utf-8
"""
This module provides an object oriented interface for pattern matching
of files.
"""

try:
	from typing import (
		Any,
		AnyStr,
		Callable,
		Iterable,
		Iterator,
		Optional,
		Text,
		Union)
except ImportError:
	pass

try:
	# Python 3.6+ type hints.
	from os import PathLike
	from typing import Collection
except ImportError:
	pass

from . import util
from .compat import (
	CollectionType,
	iterkeys,
	izip_longest,
	string_types)
from .pattern import Pattern
from .util import TreeEntry


class PathSpec(object):
	"""
	The :class:`PathSpec` class is a wrapper around a list of compiled
	:class:`.Pattern` instances.
	"""

	def __init__(self, patterns):
		# type: (Iterable[Pattern]) -> None
		"""
		Initializes the :class:`PathSpec` instance.

		*patterns* (:class:`~collections.abc.Collection` or :class:`~collections.abc.Iterable`)
		yields each compiled pattern (:class:`.Pattern`).
		"""

		self.patterns = patterns if isinstance(patterns, CollectionType) else list(patterns)
		"""
		*patterns* (:class:`~collections.abc.Collection` of :class:`.Pattern`)
		contains the compiled patterns.
		"""

	def __eq__(self, other):
		# type: (PathSpec) -> bool
		"""
		Tests the equality of this path-spec with *other* (:class:`PathSpec`)
		by comparing their :attr:`~PathSpec.patterns` attributes.
		"""
		if isinstance(other, PathSpec):
			paired_patterns = izip_longest(self.patterns, other.patterns)
			return all(a == b for a, b in paired_patterns)
		else:
			return NotImplemented

	def __len__(self):
		"""
		Returns the number of compiled patterns this path-spec contains
		(:class:`int`).
		"""
		return len(self.patterns)

	def __add__(self, other):
		# type: (PathSpec) -> PathSpec
		"""
		Combines the :attr:`Pathspec.patterns` patterns from two
		:class:`PathSpec` instances.
		"""
		if isinstance(other, PathSpec):
			return PathSpec(self.patterns + other.patterns)
		else:
			return NotImplemented

	def __iadd__(self, other):
		# type: (PathSpec) -> PathSpec
		"""
		Adds the :attr:`Pathspec.patterns` patterns from one :class:`PathSpec`
		instance to this instance.
		"""
		if isinstance(other, PathSpec):
			self.patterns += other.patterns
			return self
		else:
			return NotImplemented

	@classmethod
	def from_lines(cls, pattern_factory, lines):
		# type: (Union[Text, Callable[[AnyStr], Pattern]], Iterable[AnyStr]) -> PathSpec
		"""
		Compiles the pattern lines.

		*pattern_factory* can be either the name of a registered pattern
		factory (:class:`str`), or a :class:`~collections.abc.Callable` used
		to compile patterns. It must accept an uncompiled pattern (:class:`str`)
		and return the compiled pattern (:class:`.Pattern`).

		*lines* (:class:`~collections.abc.Iterable`) yields each uncompiled
		pattern (:class:`str`). This simply has to yield each line so it can
		be a :class:`file` (e.g., from :func:`open` or :class:`io.StringIO`)
		or the result from :meth:`str.splitlines`.

		Returns the :class:`PathSpec` instance.
		"""
		if isinstance(pattern_factory, string_types):
			pattern_factory = util.lookup_pattern(pattern_factory)
		if not callable(pattern_factory):
			raise TypeError("pattern_factory:{!r} is not callable.".format(pattern_factory))

		if not util._is_iterable(lines):
			raise TypeError("lines:{!r} is not an iterable.".format(lines))

		patterns = [pattern_factory(line) for line in lines if line]
		return cls(patterns)

	def match_file(self, file, separators=None):
		# type: (Union[Text, PathLike], Optional[Collection[Text]]) -> bool
		"""
		Matches the file to this path-spec.

		*file* (:class:`str` or :class:`~pathlib.PurePath`) is the file path
		to be matched against :attr:`self.patterns <PathSpec.patterns>`.

		*separators* (:class:`~collections.abc.Collection` of :class:`str`)
		optionally contains the path separators to normalize. See
		:func:`~pathspec.util.normalize_file` for more information.

		Returns :data:`True` if *file* matched; otherwise, :data:`False`.
		"""
		norm_file = util.normalize_file(file, separators=separators)
		return util.match_file(self.patterns, norm_file)

	def match_entries(self, entries, separators=None):
		# type: (Iterable[TreeEntry], Optional[Collection[Text]]) -> Iterator[TreeEntry]
		"""
		Matches the entries to this path-spec.

		*entries* (:class:`~collections.abc.Iterable` of :class:`~util.TreeEntry`)
		contains the entries to be matched against :attr:`self.patterns <PathSpec.patterns>`.

		*separators* (:class:`~collections.abc.Collection` of :class:`str`;
		or :data:`None`) optionally contains the path separators to
		normalize. See :func:`~pathspec.util.normalize_file` for more
		information.

		Returns the matched entries (:class:`~collections.abc.Iterator` of
		:class:`~util.TreeEntry`).
		"""
		if not util._is_iterable(entries):
			raise TypeError("entries:{!r} is not an iterable.".format(entries))

		entry_map = util._normalize_entries(entries, separators=separators)
		match_paths = util.match_files(self.patterns, iterkeys(entry_map))
		for path in match_paths:
			yield entry_map[path]

	def match_files(self, files, separators=None):
		# type: (Iterable[Union[Text, PathLike]], Optional[Collection[Text]]) -> Iterator[Union[Text, PathLike]]
		"""
		Matches the files to this path-spec.

		*files* (:class:`~collections.abc.Iterable` of :class:`str; or
		:class:`pathlib.PurePath`) contains the file paths to be matched
		against :attr:`self.patterns <PathSpec.patterns>`.

		*separators* (:class:`~collections.abc.Collection` of :class:`str`;
		or :data:`None`) optionally contains the path separators to
		normalize. See :func:`~pathspec.util.normalize_file` for more
		information.

		Returns the matched files (:class:`~collections.abc.Iterator` of
		:class:`str` or :class:`pathlib.PurePath`).
		"""
		if not util._is_iterable(files):
			raise TypeError("files:{!r} is not an iterable.".format(files))

		file_map = util.normalize_files(files, separators=separators)
		matched_files = util.match_files(self.patterns, iterkeys(file_map))
		for norm_file in matched_files:
			for orig_file in file_map[norm_file]:
				yield orig_file

	def match_tree_entries(self, root, on_error=None, follow_links=None):
		# type: (Text, Optional[Callable], Optional[bool]) -> Iterator[TreeEntry]
		"""
		Walks the specified root path for all files and matches them to this
		path-spec.

		*root* (:class:`str`; or :class:`pathlib.PurePath`) is the root
		directory to search.

		*on_error* (:class:`~collections.abc.Callable` or :data:`None`)
		optionally is the error handler for file-system exceptions. See
		:func:`~pathspec.util.iter_tree_entries` for more information.

		*follow_links* (:class:`bool` or :data:`None`) optionally is whether
		to walk symbolic links that resolve to directories. See
		:func:`~pathspec.util.iter_tree_files` for more information.

		Returns the matched files (:class:`~collections.abc.Iterator` of
		:class:`.TreeEntry`).
		"""
		entries = util.iter_tree_entries(root, on_error=on_error, follow_links=follow_links)
		return self.match_entries(entries)

	def match_tree_files(self, root, on_error=None, follow_links=None):
		# type: (Text, Optional[Callable], Optional[bool]) -> Iterator[Text]
		"""
		Walks the specified root path for all files and matches them to this
		path-spec.

		*root* (:class:`str`; or :class:`pathlib.PurePath`) is the root
		directory to search for files.

		*on_error* (:class:`~collections.abc.Callable` or :data:`None`)
		optionally is the error handler for file-system exceptions. See
		:func:`~pathspec.util.iter_tree_files` for more information.

		*follow_links* (:class:`bool` or :data:`None`) optionally is whether
		to walk symbolic links that resolve to directories. See
		:func:`~pathspec.util.iter_tree_files` for more information.

		Returns the matched files (:class:`~collections.abc.Iterable` of
		:class:`str`).
		"""
		files = util.iter_tree_files(root, on_error=on_error, follow_links=follow_links)
		return self.match_files(files)

	# Alias `match_tree_files()` as `match_tree()`.
	match_tree = match_tree_files
