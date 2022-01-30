# encoding: utf-8
"""
This module provides utility methods for dealing with path-specs.
"""

import os
import os.path
import posixpath
import stat
try:
	from typing import (
		Any,
		AnyStr,
		Callable,
		Dict,
		Iterable,
		Iterator,
		List,
		Optional,
		Sequence,
		Set,
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

from .compat import (
	CollectionType,
	IterableType,
	string_types,
	unicode)
from .pattern import Pattern

NORMALIZE_PATH_SEPS = [sep for sep in [os.sep, os.altsep] if sep and sep != posixpath.sep]
"""
*NORMALIZE_PATH_SEPS* (:class:`list` of :class:`str`) contains the path
separators that need to be normalized to the POSIX separator for the
current operating system. The separators are determined by examining
:data:`os.sep` and :data:`os.altsep`.
"""

_registered_patterns = {}
"""
*_registered_patterns* (:class:`dict`) maps a name (:class:`str`) to the
registered pattern factory (:class:`~collections.abc.Callable`).
"""


def detailed_match_files(patterns, files, all_matches=None):
	# type: (Iterable[Pattern], Iterable[Text], Optional[bool]) -> Dict[Text, 'MatchDetail']
	"""
	Matches the files to the patterns, and returns which patterns matched
	the files.

	*patterns* (:class:`~collections.abc.Iterable` of :class:`~pathspec.pattern.Pattern`)
	contains the patterns to use.

	*files* (:class:`~collections.abc.Iterable` of :class:`str`) contains
	the normalized file paths to be matched against *patterns*.

	*all_matches* (:class:`boot` or :data:`None`) is whether to return all
	matches patterns (:data:`True`), or only the last matched pattern
	(:data:`False`). Default is :data:`None` for :data:`False`.

	Returns the matched files (:class:`dict`) which maps each matched file
	(:class:`str`) to the patterns that matched in order (:class:`.MatchDetail`).
	"""
	all_files = files if isinstance(files, CollectionType) else list(files)
	return_files = {}
	for pattern in patterns:
		if pattern.include is not None:
			result_files = pattern.match(all_files)
			if pattern.include:
				# Add files and record pattern.
				for result_file in result_files:
					if result_file in return_files:
						if all_matches:
							return_files[result_file].patterns.append(pattern)
						else:
							return_files[result_file].patterns[0] = pattern
					else:
						return_files[result_file] = MatchDetail([pattern])

			else:
				# Remove files.
				for file in result_files:
					del return_files[file]

	return return_files


def _is_iterable(value):
	# type: (Any) -> bool
	"""
	Check whether the value is an iterable (excludes strings).

	*value* is the value to check,

	Returns whether *value* is a iterable (:class:`bool`).
	"""
	return isinstance(value, IterableType) and not isinstance(value, (unicode, bytes))


def iter_tree_entries(root, on_error=None, follow_links=None):
	# type: (Text, Optional[Callable], Optional[bool]) -> Iterator['TreeEntry']
	"""
	Walks the specified directory for all files and directories.

	*root* (:class:`str`) is the root directory to search.

	*on_error* (:class:`~collections.abc.Callable` or :data:`None`)
	optionally is the error handler for file-system exceptions. It will be
	called with the exception (:exc:`OSError`). Reraise the exception to
	abort the walk. Default is :data:`None` to ignore file-system
	exceptions.

	*follow_links* (:class:`bool` or :data:`None`) optionally is whether
	to walk symbolic links that resolve to directories. Default is
	:data:`None` for :data:`True`.

	Raises :exc:`RecursionError` if recursion is detected.

	Returns an :class:`~collections.abc.Iterator` yielding each file or
	directory entry (:class:`.TreeEntry`) relative to *root*.
	"""
	if on_error is not None and not callable(on_error):
		raise TypeError("on_error:{!r} is not callable.".format(on_error))

	if follow_links is None:
		follow_links = True

	for entry in _iter_tree_entries_next(os.path.abspath(root), '', {}, on_error, follow_links):
		yield entry


def iter_tree_files(root, on_error=None, follow_links=None):
	# type: (Text, Optional[Callable], Optional[bool]) -> Iterator[Text]
	"""
	Walks the specified directory for all files.

	*root* (:class:`str`) is the root directory to search for files.

	*on_error* (:class:`~collections.abc.Callable` or :data:`None`)
	optionally is the error handler for file-system exceptions. It will be
	called with the exception (:exc:`OSError`). Reraise the exception to
	abort the walk. Default is :data:`None` to ignore file-system
	exceptions.

	*follow_links* (:class:`bool` or :data:`None`) optionally is whether
	to walk symbolic links that resolve to directories. Default is
	:data:`None` for :data:`True`.

	Raises :exc:`RecursionError` if recursion is detected.

	Returns an :class:`~collections.abc.Iterator` yielding the path to
	each file (:class:`str`) relative to *root*.
	"""
	if on_error is not None and not callable(on_error):
		raise TypeError("on_error:{!r} is not callable.".format(on_error))

	if follow_links is None:
		follow_links = True

	for entry in _iter_tree_entries_next(os.path.abspath(root), '', {}, on_error, follow_links):
		if not entry.is_dir(follow_links):
			yield entry.path


# Alias `iter_tree_files()` as `iter_tree()`.
iter_tree = iter_tree_files


def _iter_tree_entries_next(root_full, dir_rel, memo, on_error, follow_links):
	# type: (Text, Text, Dict[Text, Text], Callable, bool) -> Iterator['TreeEntry']
	"""
	Scan the directory for all descendant files.

	*root_full* (:class:`str`) the absolute path to the root directory.

	*dir_rel* (:class:`str`) the path to the directory to scan relative to
	*root_full*.

	*memo* (:class:`dict`) keeps track of ancestor directories
	encountered. Maps each ancestor real path (:class:`str`) to relative
	path (:class:`str`).

	*on_error* (:class:`~collections.abc.Callable` or :data:`None`)
	optionally is the error handler for file-system exceptions.

	*follow_links* (:class:`bool`) is whether to walk symbolic links that
	resolve to directories.

	Yields each entry (:class:`.TreeEntry`).
	"""
	dir_full = os.path.join(root_full, dir_rel)
	dir_real = os.path.realpath(dir_full)

	# Remember each encountered ancestor directory and its canonical
	# (real) path. If a canonical path is encountered more than once,
	# recursion has occurred.
	if dir_real not in memo:
		memo[dir_real] = dir_rel
	else:
		raise RecursionError(real_path=dir_real, first_path=memo[dir_real], second_path=dir_rel)

	for node_name in os.listdir(dir_full):
		node_rel = os.path.join(dir_rel, node_name)
		node_full = os.path.join(root_full, node_rel)

		# Inspect child node.
		try:
			node_lstat = os.lstat(node_full)
		except OSError as e:
			if on_error is not None:
				on_error(e)
			continue

		if stat.S_ISLNK(node_lstat.st_mode):
			# Child node is a link, inspect the target node.
			is_link = True
			try:
				node_stat = os.stat(node_full)
			except OSError as e:
				if on_error is not None:
					on_error(e)
				continue
		else:
			is_link = False
			node_stat = node_lstat

		if stat.S_ISDIR(node_stat.st_mode) and (follow_links or not is_link):
			# Child node is a directory, recurse into it and yield its
			# descendant files.
			yield TreeEntry(node_name, node_rel, node_lstat, node_stat)

			for entry in _iter_tree_entries_next(root_full, node_rel, memo, on_error, follow_links):
				yield entry

		elif stat.S_ISREG(node_stat.st_mode) or is_link:
			# Child node is either a file or an unfollowed link, yield it.
			yield TreeEntry(node_name, node_rel, node_lstat, node_stat)

	# NOTE: Make sure to remove the canonical (real) path of the directory
	# from the ancestors memo once we are done with it. This allows the
	# same directory to appear multiple times. If this is not done, the
	# second occurrence of the directory will be incorrectly interpreted
	# as a recursion. See <https://github.com/cpburnz/python-path-specification/pull/7>.
	del memo[dir_real]


def lookup_pattern(name):
	# type: (Text) -> Callable[[AnyStr], Pattern]
	"""
	Lookups a registered pattern factory by name.

	*name* (:class:`str`) is the name of the pattern factory.

	Returns the registered pattern factory (:class:`~collections.abc.Callable`).
	If no pattern factory is registered, raises :exc:`KeyError`.
	"""
	return _registered_patterns[name]


def match_file(patterns, file):
	# type: (Iterable[Pattern], Text) -> bool
	"""
	Matches the file to the patterns.

	*patterns* (:class:`~collections.abc.Iterable` of :class:`~pathspec.pattern.Pattern`)
	contains the patterns to use.

	*file* (:class:`str`) is the normalized file path to be matched
	against *patterns*.

	Returns :data:`True` if *file* matched; otherwise, :data:`False`.
	"""
	matched = False
	for pattern in patterns:
		if pattern.include is not None:
			if file in pattern.match((file,)):
				matched = pattern.include
	return matched


def match_files(patterns, files):
	# type: (Iterable[Pattern], Iterable[Text]) -> Set[Text]
	"""
	Matches the files to the patterns.

	*patterns* (:class:`~collections.abc.Iterable` of :class:`~pathspec.pattern.Pattern`)
	contains the patterns to use.

	*files* (:class:`~collections.abc.Iterable` of :class:`str`) contains
	the normalized file paths to be matched against *patterns*.

	Returns the matched files (:class:`set` of :class:`str`).
	"""
	all_files = files if isinstance(files, CollectionType) else list(files)
	return_files = set()
	for pattern in patterns:
		if pattern.include is not None:
			result_files = pattern.match(all_files)
			if pattern.include:
				return_files.update(result_files)
			else:
				return_files.difference_update(result_files)
	return return_files


def _normalize_entries(entries, separators=None):
	# type: (Iterable['TreeEntry'], Optional[Collection[Text]]) -> Dict[Text, 'TreeEntry']
	"""
	Normalizes the entry paths to use the POSIX path separator.

	*entries* (:class:`~collections.abc.Iterable` of :class:`.TreeEntry`)
	contains the entries to be normalized.

	*separators* (:class:`~collections.abc.Collection` of :class:`str`; or
	:data:`None`) optionally contains the path separators to normalize.
	See :func:`normalize_file` for more information.

	Returns a :class:`dict` mapping the each normalized file path (:class:`str`)
	to the entry (:class:`.TreeEntry`)
	"""
	norm_files = {}
	for entry in entries:
		norm_files[normalize_file(entry.path, separators=separators)] = entry
	return norm_files


def normalize_file(file, separators=None):
	# type: (Union[Text, PathLike], Optional[Collection[Text]]) -> Text
	"""
	Normalizes the file path to use the POSIX path separator (i.e.,
	``'/'``), and make the paths relative (remove leading ``'/'``).

	*file* (:class:`str` or :class:`pathlib.PurePath`) is the file path.

	*separators* (:class:`~collections.abc.Collection` of :class:`str`; or
	:data:`None`) optionally contains the path separators to normalize.
	This does not need to include the POSIX path separator (``'/'``), but
	including it will not affect the results. Default is :data:`None` for
	:data:`NORMALIZE_PATH_SEPS`. To prevent normalization, pass an empty
	container (e.g., an empty tuple ``()``).

	Returns the normalized file path (:class:`str`).
	"""
	# Normalize path separators.
	if separators is None:
		separators = NORMALIZE_PATH_SEPS

	# Convert path object to string.
	norm_file = str(file)

	for sep in separators:
		norm_file = norm_file.replace(sep, posixpath.sep)

	if norm_file.startswith('/'):
		# Make path relative.
		norm_file = norm_file[1:]

	elif norm_file.startswith('./'):
		# Remove current directory prefix.
		norm_file = norm_file[2:]

	return norm_file


def normalize_files(files, separators=None):
	# type: (Iterable[Union[str, PathLike]], Optional[Collection[Text]]) -> Dict[Text, List[Union[str, PathLike]]]
	"""
	Normalizes the file paths to use the POSIX path separator.

	*files* (:class:`~collections.abc.Iterable` of :class:`str` or
	:class:`pathlib.PurePath`) contains the file paths to be normalized.

	*separators* (:class:`~collections.abc.Collection` of :class:`str`; or
	:data:`None`) optionally contains the path separators to normalize.
	See :func:`normalize_file` for more information.

	Returns a :class:`dict` mapping the each normalized file path
	(:class:`str`) to the original file paths (:class:`list` of
	:class:`str` or :class:`pathlib.PurePath`).
	"""
	norm_files = {}
	for path in files:
		norm_file = normalize_file(path, separators=separators)
		if norm_file in norm_files:
			norm_files[norm_file].append(path)
		else:
			norm_files[norm_file] = [path]

	return norm_files


def register_pattern(name, pattern_factory, override=None):
	# type: (Text, Callable[[AnyStr], Pattern], Optional[bool]) -> None
	"""
	Registers the specified pattern factory.

	*name* (:class:`str`) is the name to register the pattern factory
	under.

	*pattern_factory* (:class:`~collections.abc.Callable`) is used to
	compile patterns. It must accept an uncompiled pattern (:class:`str`)
	and return the compiled pattern (:class:`.Pattern`).

	*override* (:class:`bool` or :data:`None`) optionally is whether to
	allow overriding an already registered pattern under the same name
	(:data:`True`), instead of raising an :exc:`AlreadyRegisteredError`
	(:data:`False`). Default is :data:`None` for :data:`False`.
	"""
	if not isinstance(name, string_types):
		raise TypeError("name:{!r} is not a string.".format(name))
	if not callable(pattern_factory):
		raise TypeError("pattern_factory:{!r} is not callable.".format(pattern_factory))
	if name in _registered_patterns and not override:
		raise AlreadyRegisteredError(name, _registered_patterns[name])
	_registered_patterns[name] = pattern_factory


class AlreadyRegisteredError(Exception):
	"""
	The :exc:`AlreadyRegisteredError` exception is raised when a pattern
	factory is registered under a name already in use.
	"""

	def __init__(self, name, pattern_factory):
		# type: (Text, Callable[[AnyStr], Pattern]) -> None
		"""
		Initializes the :exc:`AlreadyRegisteredError` instance.

		*name* (:class:`str`) is the name of the registered pattern.

		*pattern_factory* (:class:`~collections.abc.Callable`) is the
		registered pattern factory.
		"""
		super(AlreadyRegisteredError, self).__init__(name, pattern_factory)

	@property
	def message(self):
		# type: () -> Text
		"""
		*message* (:class:`str`) is the error message.
		"""
		return "{name!r} is already registered for pattern factory:{pattern_factory!r}.".format(
			name=self.name,
			pattern_factory=self.pattern_factory,
		)

	@property
	def name(self):
		# type: () -> Text
		"""
		*name* (:class:`str`) is the name of the registered pattern.
		"""
		return self.args[0]

	@property
	def pattern_factory(self):
		# type: () -> Callable[[AnyStr], Pattern]
		"""
		*pattern_factory* (:class:`~collections.abc.Callable`) is the
		registered pattern factory.
		"""
		return self.args[1]


class RecursionError(Exception):
	"""
	The :exc:`RecursionError` exception is raised when recursion is
	detected.
	"""

	def __init__(self, real_path, first_path, second_path):
		# type: (Text, Text, Text) -> None
		"""
		Initializes the :exc:`RecursionError` instance.

		*real_path* (:class:`str`) is the real path that recursion was
		encountered on.

		*first_path* (:class:`str`) is the first path encountered for
		*real_path*.

		*second_path* (:class:`str`) is the second path encountered for
		*real_path*.
		"""
		super(RecursionError, self).__init__(real_path, first_path, second_path)

	@property
	def first_path(self):
		# type: () -> Text
		"""
		*first_path* (:class:`str`) is the first path encountered for
		:attr:`self.real_path <RecursionError.real_path>`.
		"""
		return self.args[1]

	@property
	def message(self):
		# type: () -> Text
		"""
		*message* (:class:`str`) is the error message.
		"""
		return "Real path {real!r} was encountered at {first!r} and then {second!r}.".format(
			real=self.real_path,
			first=self.first_path,
			second=self.second_path,
		)

	@property
	def real_path(self):
		# type: () -> Text
		"""
		*real_path* (:class:`str`) is the real path that recursion was
		encountered on.
		"""
		return self.args[0]

	@property
	def second_path(self):
		# type: () -> Text
		"""
		*second_path* (:class:`str`) is the second path encountered for
		:attr:`self.real_path <RecursionError.real_path>`.
		"""
		return self.args[2]


class MatchDetail(object):
	"""
	The :class:`.MatchDetail` class contains information about
	"""

	#: Make the class dict-less.
	__slots__ = ('patterns',)

	def __init__(self, patterns):
		# type: (Sequence[Pattern]) -> None
		"""
		Initialize the :class:`.MatchDetail` instance.

		*patterns* (:class:`~collections.abc.Sequence` of :class:`~pathspec.pattern.Pattern`)
		contains the patterns that matched the file in the order they were
		encountered.
		"""

		self.patterns = patterns
		"""
		*patterns* (:class:`~collections.abc.Sequence` of :class:`~pathspec.pattern.Pattern`)
		contains the patterns that matched the file in the order they were
		encountered.
		"""


class TreeEntry(object):
	"""
	The :class:`.TreeEntry` class contains information about a file-system
	entry.
	"""

	#: Make the class dict-less.
	__slots__ = ('_lstat', 'name', 'path', '_stat')

	def __init__(self, name, path, lstat, stat):
		# type: (Text, Text, os.stat_result, os.stat_result) -> None
		"""
		Initialize the :class:`.TreeEntry` instance.

		*name* (:class:`str`) is the base name of the entry.

		*path* (:class:`str`) is the relative path of the entry.

		*lstat* (:class:`~os.stat_result`) is the stat result of the direct
		entry.

		*stat* (:class:`~os.stat_result`) is the stat result of the entry,
		potentially linked.
		"""

		self._lstat = lstat
		"""
		*_lstat* (:class:`~os.stat_result`) is the stat result of the direct
		entry.
		"""

		self.name = name
		"""
		*name* (:class:`str`) is the base name of the entry.
		"""

		self.path = path
		"""
		*path* (:class:`str`) is the path of the entry.
		"""

		self._stat = stat
		"""
		*_stat* (:class:`~os.stat_result`) is the stat result of the linked
		entry.
		"""

	def is_dir(self, follow_links=None):
		# type: (Optional[bool]) -> bool
		"""
		Get whether the entry is a directory.

		*follow_links* (:class:`bool` or :data:`None`) is whether to follow
		symbolic links. If this is :data:`True`, a symlink to a directory
		will result in :data:`True`. Default is :data:`None` for :data:`True`.

		Returns whether the entry is a directory (:class:`bool`).
		"""
		if follow_links is None:
			follow_links = True

		node_stat = self._stat if follow_links else self._lstat
		return stat.S_ISDIR(node_stat.st_mode)

	def is_file(self, follow_links=None):
		# type: (Optional[bool]) -> bool
		"""
		Get whether the entry is a regular file.

		*follow_links* (:class:`bool` or :data:`None`) is whether to follow
		symbolic links. If this is :data:`True`, a symlink to a regular file
		will result in :data:`True`. Default is :data:`None` for :data:`True`.

		Returns whether the entry is a regular file (:class:`bool`).
		"""
		if follow_links is None:
			follow_links = True

		node_stat = self._stat if follow_links else self._lstat
		return stat.S_ISREG(node_stat.st_mode)

	def is_symlink(self):
		# type: () -> bool
		"""
		Returns whether the entry is a symbolic link (:class:`bool`).
		"""
		return stat.S_ISLNK(self._lstat.st_mode)

	def stat(self, follow_links=None):
		# type: (Optional[bool]) -> os.stat_result
		"""
		Get the cached stat result for the entry.

		*follow_links* (:class:`bool` or :data:`None`) is whether to follow
		symbolic links. If this is :data:`True`, the stat result of the
		linked file will be returned. Default is :data:`None` for :data:`True`.

		Returns that stat result (:class:`~os.stat_result`).
		"""
		if follow_links is None:
			follow_links = True

		return self._stat if follow_links else self._lstat
