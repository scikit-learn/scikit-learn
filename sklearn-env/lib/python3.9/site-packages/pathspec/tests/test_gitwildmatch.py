# encoding: utf-8
"""
This script tests ``GitWildMatchPattern``.
"""
from __future__ import unicode_literals

import re
import sys
import unittest

import pathspec.patterns.gitwildmatch
import pathspec.util
from pathspec.patterns.gitwildmatch import GitWildMatchPattern, GitWildMatchPatternError

if sys.version_info[0] >= 3:
	unichr = chr


class GitWildMatchTest(unittest.TestCase):
	"""
	The ``GitWildMatchTest`` class tests the ``GitWildMatchPattern``
	implementation.
	"""

	def test_00_empty(self):
		"""
		Tests an empty pattern.
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('')
		self.assertIsNone(include)
		self.assertIsNone(regex)

	def test_01_absolute(self):
		"""
		Tests an absolute path pattern.

		This should match:

			an/absolute/file/path
			an/absolute/file/path/foo

		This should NOT match:

			foo/an/absolute/file/path
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('/an/absolute/file/path')
		self.assertTrue(include)
		self.assertEqual(regex, '^an/absolute/file/path(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'an/absolute/file/path',
			'an/absolute/file/path/foo',
			'foo/an/absolute/file/path',
		]))
		self.assertEqual(results, {
			'an/absolute/file/path',
			'an/absolute/file/path/foo',
		})

	def test_01_absolute_root(self):
		"""
		Tests a single root absolute path pattern.

		This should NOT match any file (according to git check-ignore
		(v2.4.1)).
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('/')
		self.assertIsNone(include)
		self.assertIsNone(regex)

	def test_01_relative(self):
		"""
		Tests a relative path pattern.

		This should match:

			spam
			spam/
			foo/spam
			spam/foo
			foo/spam/bar
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('spam')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?spam(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'spam',
			'spam/',
			'foo/spam',
			'spam/foo',
			'foo/spam/bar',
		]))
		self.assertEqual(results, {
			'spam',
			'spam/',
			'foo/spam',
			'spam/foo',
			'foo/spam/bar',
		})

	def test_01_relative_nested(self):
		"""
		Tests a relative nested path pattern.

		This should match:

			foo/spam
			foo/spam/bar

		This should **not** match (according to git check-ignore (v2.4.1)):

			bar/foo/spam
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('foo/spam')
		self.assertTrue(include)
		self.assertEqual(regex, '^foo/spam(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'foo/spam',
			'foo/spam/bar',
			'bar/foo/spam',
		]))
		self.assertEqual(results, {
			'foo/spam',
			'foo/spam/bar',
		})

	def test_02_comment(self):
		"""
		Tests a comment pattern.
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('# Cork soakers.')
		self.assertIsNone(include)
		self.assertIsNone(regex)

	def test_02_ignore(self):
		"""
		Tests an exclude pattern.

		This should NOT match (according to git check-ignore (v2.4.1)):

			temp/foo
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('!temp')
		self.assertIsNotNone(include)
		self.assertFalse(include)
		self.assertEqual(regex, '^(?:.+/)?temp$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match(['temp/foo']))
		self.assertEqual(results, set())

	def test_03_child_double_asterisk(self):
		"""
		Tests a directory name with a double-asterisk child
		directory.

		This should match:

			spam/bar

		This should **not** match (according to git check-ignore (v2.4.1)):

			foo/spam/bar
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('spam/**')
		self.assertTrue(include)
		self.assertEqual(regex, '^spam/.*$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'spam/bar',
			'foo/spam/bar',
		]))
		self.assertEqual(results, {'spam/bar'})

	def test_03_inner_double_asterisk(self):
		"""
		Tests a path with an inner double-asterisk directory.

		This should match:

			left/right
			left/bar/right
			left/foo/bar/right
			left/bar/right/foo

		This should **not** match (according to git check-ignore (v2.4.1)):

			foo/left/bar/right
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('left/**/right')
		self.assertTrue(include)
		self.assertEqual(regex, '^left(?:/.+)?/right(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'left/right',
			'left/bar/right',
			'left/foo/bar/right',
			'left/bar/right/foo',
			'foo/left/bar/right',
		]))
		self.assertEqual(results, {
			'left/right',
			'left/bar/right',
			'left/foo/bar/right',
			'left/bar/right/foo',
		})

	def test_03_only_double_asterisk(self):
		"""
		Tests a double-asterisk pattern which matches everything.
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('**')
		self.assertTrue(include)
		self.assertEqual(regex, '^.+$')

	def test_03_parent_double_asterisk(self):
		"""
		Tests a file name with a double-asterisk parent directory.

		This should match:

			spam
			foo/spam
			foo/spam/bar
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('**/spam')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?spam(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'spam',
			'foo/spam',
			'foo/spam/bar',
		]))
		self.assertEqual(results, {
			'spam',
			'foo/spam',
			'foo/spam/bar',
		})

	def test_03_duplicate_leading_double_asterisk_edge_case(self):
		"""
		Regression test for duplicate leading **/ bug.
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('**')
		self.assertTrue(include)
		self.assertEqual(regex, '^.+$')

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/**')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/**/**')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

		regex, include = GitWildMatchPattern.pattern_to_regex('**/api')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?api(?:/.*)?$')

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/**/api')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

		regex, include = GitWildMatchPattern.pattern_to_regex('**/api/')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?api/.*$')

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/api/**')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/**/api/**/**')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

	def test_03_double_asterisk_trailing_slash_edge_case(self):
		"""
		Tests the edge-case **/ pattern.

		This should match everything except individual files in the root directory.
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('**/')
		self.assertTrue(include)
		self.assertEqual(regex, '^.+/.*$')

		equivalent_regex, include = GitWildMatchPattern.pattern_to_regex('**/**/')
		self.assertTrue(include)
		self.assertEqual(equivalent_regex, regex)

	def test_04_infix_wildcard(self):
		"""
		Tests a pattern with an infix wildcard.

		This should match:

			foo--bar
			foo-hello-bar
			a/foo-hello-bar
			foo-hello-bar/b
			a/foo-hello-bar/b
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('foo-*-bar')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?foo\\-[^/]*\\-bar(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'foo--bar',
			'foo-hello-bar',
			'a/foo-hello-bar',
			'foo-hello-bar/b',
			'a/foo-hello-bar/b',
		]))
		self.assertEqual(results, {
			'foo--bar',
			'foo-hello-bar',
			'a/foo-hello-bar',
			'foo-hello-bar/b',
			'a/foo-hello-bar/b',
		})

	def test_04_postfix_wildcard(self):
		"""
		Tests a pattern with a postfix wildcard.

		This should match:

			~temp-
			~temp-foo
			~temp-foo/bar
			foo/~temp-bar
			foo/~temp-bar/baz
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('~temp-*')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?\\~temp\\-[^/]*(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'~temp-',
			'~temp-foo',
			'~temp-foo/bar',
			'foo/~temp-bar',
			'foo/~temp-bar/baz',
		]))
		self.assertEqual(results, {
			'~temp-',
			'~temp-foo',
			'~temp-foo/bar',
			'foo/~temp-bar',
			'foo/~temp-bar/baz',
		})

	def test_04_prefix_wildcard(self):
		"""
		Tests a pattern with a prefix wildcard.

		This should match:

			bar.py
			bar.py/
			foo/bar.py
			foo/bar.py/baz
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('*.py')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?[^/]*\\.py(?:/.*)?$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'bar.py',
			'bar.py/',
			'foo/bar.py',
			'foo/bar.py/baz',
		]))
		self.assertEqual(results, {
			'bar.py',
			'bar.py/',
			'foo/bar.py',
			'foo/bar.py/baz',
		})

	def test_05_directory(self):
		"""
		Tests a directory pattern.

		This should match:

			dir/
			foo/dir/
			foo/dir/bar

		This should **not** match:

			dir
		"""
		regex, include = GitWildMatchPattern.pattern_to_regex('dir/')
		self.assertTrue(include)
		self.assertEqual(regex, '^(?:.+/)?dir/.*$')

		pattern = GitWildMatchPattern(re.compile(regex), include)
		results = set(pattern.match([
			'dir/',
			'foo/dir/',
			'foo/dir/bar',
			'dir',
		]))
		self.assertEqual(results, {
			'dir/',
			'foo/dir/',
			'foo/dir/bar',
		})

	def test_06_registered(self):
		"""
		Tests that the pattern is registered.
		"""
		self.assertIs(pathspec.util.lookup_pattern('gitwildmatch'), GitWildMatchPattern)

	def test_06_access_deprecated(self):
		"""
		Tests that the pattern is accessible from the root module using the
		deprecated alias.
		"""
		self.assertTrue(hasattr(pathspec, 'GitIgnorePattern'))
		self.assertTrue(issubclass(pathspec.GitIgnorePattern, GitWildMatchPattern))

	def test_06_registered_deprecated(self):
		"""
		Tests that the pattern is registered under the deprecated alias.
		"""
		self.assertIs(pathspec.util.lookup_pattern('gitignore'), pathspec.GitIgnorePattern)

	def test_07_encode_bytes(self):
		"""
		Test encoding bytes.
		"""
		encoded = "".join(map(unichr, range(0,256))).encode(pathspec.patterns.gitwildmatch._BYTES_ENCODING)
		expected = b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff"
		self.assertEqual(encoded, expected)

	def test_07_decode_bytes(self):
		"""
		Test decoding bytes.
		"""
		decoded = bytes(bytearray(range(0,256))).decode(pathspec.patterns.gitwildmatch._BYTES_ENCODING)
		expected = "\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff"
		self.assertEqual(decoded, expected)

	def test_07_match_bytes_and_bytes(self):
		"""
		Test byte string patterns matching byte string paths.
		"""
		pattern = GitWildMatchPattern(b'*.py')
		results = set(pattern.match([b'a.py']))
		self.assertEqual(results, {b'a.py'})

	def test_07_match_bytes_and_bytes_complete(self):
		"""
		Test byte string patterns matching byte string paths.
		"""
		encoded = bytes(bytearray(range(0,256)))
		escaped = b"".join(b"\\" + encoded[i:i+1] for i in range(len(encoded)))
		pattern = GitWildMatchPattern(escaped)
		results = set(pattern.match([encoded]))
		self.assertEqual(results, {encoded})

	@unittest.skipIf(sys.version_info[0] >= 3, "Python 3 is strict")
	def test_07_match_bytes_and_unicode(self):
		"""
		Test byte string patterns matching byte string paths.
		"""
		pattern = GitWildMatchPattern(b'*.py')
		results = set(pattern.match(['a.py']))
		self.assertEqual(results, {'a.py'})

	@unittest.skipIf(sys.version_info[0] == 2, "Python 2 is lenient")
	def test_07_match_bytes_and_unicode_fail(self):
		"""
		Test byte string patterns matching byte string paths.
		"""
		pattern = GitWildMatchPattern(b'*.py')
		with self.assertRaises(TypeError):
			for _ in pattern.match(['a.py']):
				pass

	@unittest.skipIf(sys.version_info[0] >= 3, "Python 3 is strict")
	def test_07_match_unicode_and_bytes(self):
		"""
		Test unicode patterns with byte paths.
		"""
		pattern = GitWildMatchPattern('*.py')
		results = set(pattern.match([b'a.py']))
		self.assertEqual(results, {b'a.py'})

	@unittest.skipIf(sys.version_info[0] == 2, "Python 2 is lenient")
	def test_07_match_unicode_and_bytes_fail(self):
		"""
		Test unicode patterns with byte paths.
		"""
		pattern = GitWildMatchPattern('*.py')
		with self.assertRaises(TypeError):
			for _ in pattern.match([b'a.py']):
				pass

	def test_07_match_unicode_and_unicode(self):
		"""
		Test unicode patterns with unicode paths.
		"""
		pattern = GitWildMatchPattern('*.py')
		results = set(pattern.match(['a.py']))
		self.assertEqual(results, {'a.py'})

	def test_08_escape(self):
		"""
		Test escaping a string with meta-characters
		"""
		fname = "file!with*weird#naming_[1].t?t"
		escaped = r"file\!with\*weird\#naming_\[1\].t\?t"
		result = GitWildMatchPattern.escape(fname)
		self.assertEqual(result, escaped)

	def test_09_single_escape_fail(self):
		"""
		Test an escape on a line by itself.
		"""
		self._check_invalid_pattern("\\")

	def test_09_single_exclamation_mark_fail(self):
		"""
		Test an escape on a line by itself.
		"""
		self._check_invalid_pattern("!")

	def _check_invalid_pattern(self, git_ignore_pattern):
		expected_message_pattern = re.escape(repr(git_ignore_pattern))
		with self.assertRaisesRegexp(GitWildMatchPatternError, expected_message_pattern):
			GitWildMatchPattern(git_ignore_pattern)

