# encoding: utf-8
"""
This module provides the base definition for patterns.
"""

import re
try:
	from typing import (
		AnyStr,
		Iterable,
		Iterator,
		Optional,
		Pattern as RegexHint,
		Text,
		Tuple,
		Union)
except ImportError:
	pass

from .compat import unicode


class Pattern(object):
	"""
	The :class:`Pattern` class is the abstract definition of a pattern.
	"""

	# Make the class dict-less.
	__slots__ = ('include',)

	def __init__(self, include):
		# type: (Optional[bool]) -> None
		"""
		Initializes the :class:`Pattern` instance.

		*include* (:class:`bool` or :data:`None`) is whether the matched
		files should be included (:data:`True`), excluded (:data:`False`),
		or is a null-operation (:data:`None`).
		"""

		self.include = include
		"""
		*include* (:class:`bool` or :data:`None`) is whether the matched
		files should be included (:data:`True`), excluded (:data:`False`),
		or is a null-operation (:data:`None`).
		"""

	def match(self, files):
		# type: (Iterable[Text]) -> Iterator[Text]
		"""
		Matches this pattern against the specified files.

		*files* (:class:`~collections.abc.Iterable` of :class:`str`) contains
		each file relative to the root directory (e.g., ``"relative/path/to/file"``).

		Returns an :class:`~collections.abc.Iterable` yielding each matched
		file path (:class:`str`).
		"""
		raise NotImplementedError("{}.{} must override match().".format(self.__class__.__module__, self.__class__.__name__))


class RegexPattern(Pattern):
	"""
	The :class:`RegexPattern` class is an implementation of a pattern
	using regular expressions.
	"""

	# Make the class dict-less.
	__slots__ = ('regex',)

	def __init__(self, pattern, include=None):
		# type: (Union[AnyStr, RegexHint], Optional[bool]) -> None
		"""
		Initializes the :class:`RegexPattern` instance.

		*pattern* (:class:`unicode`, :class:`bytes`, :class:`re.RegexObject`,
		or :data:`None`) is the pattern to compile into a regular
		expression.

		*include* (:class:`bool` or :data:`None`) must be :data:`None`
		unless *pattern* is a precompiled regular expression (:class:`re.RegexObject`)
		in which case it is whether matched files should be included
		(:data:`True`), excluded (:data:`False`), or is a null operation
		(:data:`None`).

			.. NOTE:: Subclasses do not need to support the *include*
			   parameter.
		"""

		self.regex = None
		"""
		*regex* (:class:`re.RegexObject`) is the regular expression for the
		pattern.
		"""

		if isinstance(pattern, (unicode, bytes)):
			assert include is None, "include:{!r} must be null when pattern:{!r} is a string.".format(include, pattern)
			regex, include = self.pattern_to_regex(pattern)
			# NOTE: Make sure to allow a null regular expression to be
			# returned for a null-operation.
			if include is not None:
				regex = re.compile(regex)

		elif pattern is not None and hasattr(pattern, 'match'):
			# Assume pattern is a precompiled regular expression.
			# - NOTE: Used specified *include*.
			regex = pattern

		elif pattern is None:
			# NOTE: Make sure to allow a null pattern to be passed for a
			# null-operation.
			assert include is None, "include:{!r} must be null when pattern:{!r} is null.".format(include, pattern)

		else:
			raise TypeError("pattern:{!r} is not a string, RegexObject, or None.".format(pattern))

		super(RegexPattern, self).__init__(include)
		self.regex = regex

	def __eq__(self, other):
		# type: (RegexPattern) -> bool
		"""
		Tests the equality of this regex pattern with *other* (:class:`RegexPattern`)
		by comparing their :attr:`~Pattern.include` and :attr:`~RegexPattern.regex`
		attributes.
		"""
		if isinstance(other, RegexPattern):
			return self.include == other.include and self.regex == other.regex
		else:
			return NotImplemented

	def match(self, files):
		# type: (Iterable[Text]) -> Iterable[Text]
		"""
		Matches this pattern against the specified files.

		*files* (:class:`~collections.abc.Iterable` of :class:`str`)
		contains each file relative to the root directory (e.g., "relative/path/to/file").

		Returns an :class:`~collections.abc.Iterable` yielding each matched
		file path (:class:`str`).
		"""
		if self.include is not None:
			for path in files:
				if self.regex.match(path) is not None:
					yield path

	@classmethod
	def pattern_to_regex(cls, pattern):
		# type: (Text) -> Tuple[Text, bool]
		"""
		Convert the pattern into an uncompiled regular expression.

		*pattern* (:class:`str`) is the pattern to convert into a regular
		expression.

		Returns the uncompiled regular expression (:class:`str` or :data:`None`),
		and whether matched files should be included (:data:`True`),
		excluded (:data:`False`), or is a null-operation (:data:`None`).

			.. NOTE:: The default implementation simply returns *pattern* and
			   :data:`True`.
		"""
		return pattern, True
