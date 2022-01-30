# encoding: utf-8
"""
This script tests utility functions.
"""

import errno
import os
import os.path
import shutil
import sys
import tempfile
import unittest

from pathspec.util import iter_tree_entries, iter_tree_files, RecursionError, normalize_file


class IterTreeTest(unittest.TestCase):
	"""
	The ``IterTreeTest`` class tests `pathspec.util.iter_tree_files()`.
	"""

	def make_dirs(self, dirs):
		"""
		Create the specified directories.
		"""
		for dir in dirs:
			os.mkdir(os.path.join(self.temp_dir, self.ospath(dir)))

	def make_files(self, files):
		"""
		Create the specified files.
		"""
		for file in files:
			self.mkfile(os.path.join(self.temp_dir, self.ospath(file)))

	def make_links(self, links):
		"""
		Create the specified links.
		"""
		for link, node in links:
			os.symlink(os.path.join(self.temp_dir, self.ospath(node)), os.path.join(self.temp_dir, self.ospath(link)))

	@staticmethod
	def mkfile(file):
		"""
		Creates an empty file.
		"""
		with open(file, 'wb'):
			pass

	@staticmethod
	def ospath(path):
		"""
		Convert the POSIX path to a native OS path.
		"""
		return os.path.join(*path.split('/'))

	def require_realpath(self):
		"""
		Skips the test if `os.path.realpath` does not properly support
		symlinks.
		"""
		if self.broken_realpath:
			raise unittest.SkipTest("`os.path.realpath` is broken.")

	def require_symlink(self):
		"""
		Skips the test if `os.symlink` is not supported.
		"""
		if self.no_symlink:
			raise unittest.SkipTest("`os.symlink` is not supported.")

	def setUp(self):
		"""
		Called before each test.
		"""
		self.temp_dir = tempfile.mkdtemp()

	def tearDown(self):
		"""
		Called after each test.
		"""
		shutil.rmtree(self.temp_dir)

	def test_1_files(self):
		"""
		Tests to make sure all files are found.
		"""
		self.make_dirs([
			'Empty',
			'Dir',
			'Dir/Inner',
		])
		self.make_files([
			'a',
			'b',
			'Dir/c',
			'Dir/d',
			'Dir/Inner/e',
			'Dir/Inner/f',
		])
		results = set(iter_tree_files(self.temp_dir))
		self.assertEqual(results, set(map(self.ospath, [
			'a',
			'b',
			'Dir/c',
			'Dir/d',
			'Dir/Inner/e',
			'Dir/Inner/f',
		])))

	def test_2_0_check_symlink(self):
		"""
		Tests whether links can be created.
		"""
		# NOTE: Windows does not support `os.symlink` for Python 2. Windows
		# Vista and greater supports `os.symlink` for Python 3.2+.
		no_symlink = None
		try:
			file = os.path.join(self.temp_dir, 'file')
			link = os.path.join(self.temp_dir, 'link')
			self.mkfile(file)

			try:
				os.symlink(file, link)
			except (AttributeError, NotImplementedError):
				no_symlink = True
				raise
			no_symlink = False

		finally:
			self.__class__.no_symlink = no_symlink

	def test_2_1_check_realpath(self):
		"""
		Tests whether `os.path.realpath` works properly with symlinks.
		"""
		# NOTE: Windows does not follow symlinks with `os.path.realpath`
		# which is what we use to detect recursion. See <https://bugs.python.org/issue9949>
		# for details.
		broken_realpath = None
		try:
			self.require_symlink()
			file = os.path.join(self.temp_dir, 'file')
			link = os.path.join(self.temp_dir, 'link')
			self.mkfile(file)
			os.symlink(file, link)

			try:
				self.assertEqual(os.path.realpath(file), os.path.realpath(link))
			except AssertionError:
				broken_realpath = True
				raise
			broken_realpath = False

		finally:
			self.__class__.broken_realpath = broken_realpath

	def test_2_2_links(self):
		"""
		Tests to make sure links to directories and files work.
		"""
		self.require_symlink()
		self.make_dirs([
			'Dir',
		])
		self.make_files([
			'a',
			'b',
			'Dir/c',
			'Dir/d',
		])
		self.make_links([
			('ax', 'a'),
			('bx', 'b'),
			('Dir/cx', 'Dir/c'),
			('Dir/dx', 'Dir/d'),
			('DirX', 'Dir'),
		])
		results = set(iter_tree_files(self.temp_dir))
		self.assertEqual(results, set(map(self.ospath, [
			'a',
			'ax',
			'b',
			'bx',
			'Dir/c',
			'Dir/cx',
			'Dir/d',
			'Dir/dx',
			'DirX/c',
			'DirX/cx',
			'DirX/d',
			'DirX/dx',
		])))

	def test_2_3_sideways_links(self):
		"""
		Tests to make sure the same directory can be encountered multiple
		times via links.
		"""
		self.require_symlink()
		self.make_dirs([
			'Dir',
			'Dir/Target',
		])
		self.make_files([
			'Dir/Target/file',
		])
		self.make_links([
			('Ax', 'Dir'),
			('Bx', 'Dir'),
			('Cx', 'Dir/Target'),
			('Dx', 'Dir/Target'),
			('Dir/Ex', 'Dir/Target'),
			('Dir/Fx', 'Dir/Target'),
		])
		results = set(iter_tree_files(self.temp_dir))
		self.assertEqual(results, set(map(self.ospath, [
			'Ax/Ex/file',
			'Ax/Fx/file',
			'Ax/Target/file',
			'Bx/Ex/file',
			'Bx/Fx/file',
			'Bx/Target/file',
			'Cx/file',
			'Dx/file',
			'Dir/Ex/file',
			'Dir/Fx/file',
			'Dir/Target/file',
		])))

	def test_2_4_recursive_links(self):
		"""
		Tests detection of recursive links.
		"""
		self.require_symlink()
		self.require_realpath()
		self.make_dirs([
			'Dir',
		])
		self.make_files([
			'Dir/file',
		])
		self.make_links([
			('Dir/Self', 'Dir'),
		])
		with self.assertRaises(RecursionError) as context:
			set(iter_tree_files(self.temp_dir))
		self.assertEqual(context.exception.first_path, 'Dir')
		self.assertEqual(context.exception.second_path, self.ospath('Dir/Self'))

	def test_2_5_recursive_circular_links(self):
		"""
		Tests detection of recursion through circular links.
		"""
		self.require_symlink()
		self.require_realpath()
		self.make_dirs([
			'A',
			'B',
			'C',
		])
		self.make_files([
			'A/d',
			'B/e',
			'C/f',
		])
		self.make_links([
			('A/Bx', 'B'),
			('B/Cx', 'C'),
			('C/Ax', 'A'),
		])
		with self.assertRaises(RecursionError) as context:
			set(iter_tree_files(self.temp_dir))
		self.assertIn(context.exception.first_path, ('A', 'B', 'C'))
		self.assertEqual(context.exception.second_path, {
			'A': self.ospath('A/Bx/Cx/Ax'),
			'B': self.ospath('B/Cx/Ax/Bx'),
			'C': self.ospath('C/Ax/Bx/Cx'),
		}[context.exception.first_path])

	def test_2_6_detect_broken_links(self):
		"""
		Tests that broken links are detected.
		"""
		def reraise(e):
			raise e

		self.require_symlink()
		self.make_links([
			('A', 'DOES_NOT_EXIST'),
		])
		with self.assertRaises(OSError) as context:
			set(iter_tree_files(self.temp_dir, on_error=reraise))
		self.assertEqual(context.exception.errno, errno.ENOENT)

	def test_2_7_ignore_broken_links(self):
		"""
		Tests that broken links are ignored.
		"""
		self.require_symlink()
		self.make_links([
			('A', 'DOES_NOT_EXIST'),
		])
		results = set(iter_tree_files(self.temp_dir))
		self.assertEqual(results, set())

	def test_2_8_no_follow_links(self):
		"""
		Tests to make sure directory links can be ignored.
		"""
		self.require_symlink()
		self.make_dirs([
			'Dir',
		])
		self.make_files([
			'A',
			'B',
			'Dir/C',
			'Dir/D',
		])
		self.make_links([
			('Ax', 'A'),
			('Bx', 'B'),
			('Dir/Cx', 'Dir/C'),
			('Dir/Dx', 'Dir/D'),
			('DirX', 'Dir'),
		])
		results = set(iter_tree_files(self.temp_dir, follow_links=False))
		self.assertEqual(results, set(map(self.ospath, [
			'A',
			'Ax',
			'B',
			'Bx',
			'Dir/C',
			'Dir/Cx',
			'Dir/D',
			'Dir/Dx',
			'DirX',
		])))

	def test_3_entries(self):
		"""
		Tests to make sure all files are found.
		"""
		self.make_dirs([
			'Empty',
			'Dir',
			'Dir/Inner',
		])
		self.make_files([
			'a',
			'b',
			'Dir/c',
			'Dir/d',
			'Dir/Inner/e',
			'Dir/Inner/f',
		])
		results = {entry.path for entry in iter_tree_entries(self.temp_dir)}
		self.assertEqual(results, set(map(self.ospath, [
			'a',
			'b',
			'Dir',
			'Dir/c',
			'Dir/d',
			'Dir/Inner',
			'Dir/Inner/e',
			'Dir/Inner/f',
			'Empty',
		])))

	@unittest.skipIf(sys.version_info < (3, 4), "pathlib entered stdlib in Python 3.4")
	def test_4_normalizing_pathlib_path(self):
		"""
		Tests passing pathlib.Path as argument.
		"""
		from pathlib import Path
		first_spec = normalize_file(Path('a.txt'))
		second_spec = normalize_file('a.txt')
		self.assertEqual(first_spec, second_spec)
