# encoding: utf-8
"""Tests for genutils.path"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import sys
import tempfile

import nose.tools as nt

from ..testing.decorators import skip_if_not_win32, skip_win32
from .. import path
from .. import py3compat
from ..tempdir import TemporaryDirectory


def test_filefind():
    f = tempfile.NamedTemporaryFile()
    t = path.filefind(f.name, '.')


def test_ensure_dir_exists():
    with TemporaryDirectory() as td:
        d = os.path.join(td, u'∂ir')
        path.ensure_dir_exists(d) # create it
        assert os.path.isdir(d)
        path.ensure_dir_exists(d) # no-op
        f = os.path.join(td, u'ƒile')
        open(f, 'w').close() # touch
        with nt.assert_raises(IOError):
            path.ensure_dir_exists(f)


class TestLinkOrCopy(object):
    def setUp(self):
        self.tempdir = TemporaryDirectory()
        self.src = self.dst("src")
        with open(self.src, "w") as f:
            f.write("Hello, world!")

    def tearDown(self):
        self.tempdir.cleanup()

    def dst(self, *args):
        return os.path.join(self.tempdir.name, *args)

    def assert_inode_not_equal(self, a, b):
        nt.assert_not_equals(os.stat(a).st_ino, os.stat(b).st_ino,
                             "%r and %r do reference the same indoes" %(a, b))

    def assert_inode_equal(self, a, b):
        nt.assert_equals(os.stat(a).st_ino, os.stat(b).st_ino,
                         "%r and %r do not reference the same indoes" %(a, b))

    def assert_content_equal(self, a, b):
        with open(a) as a_f:
            with open(b) as b_f:
                nt.assert_equals(a_f.read(), b_f.read())

    @skip_win32
    def test_link_successful(self):
        dst = self.dst("target")
        path.link_or_copy(self.src, dst)
        self.assert_inode_equal(self.src, dst)

    @skip_win32
    def test_link_into_dir(self):
        dst = self.dst("some_dir")
        os.mkdir(dst)
        path.link_or_copy(self.src, dst)
        expected_dst = self.dst("some_dir", os.path.basename(self.src))
        self.assert_inode_equal(self.src, expected_dst)

    @skip_win32
    def test_target_exists(self):
        dst = self.dst("target")
        open(dst, "w").close()
        path.link_or_copy(self.src, dst)
        self.assert_inode_equal(self.src, dst)

    @skip_win32
    def test_no_link(self):
        real_link = os.link
        try:
            del os.link
            dst = self.dst("target")
            path.link_or_copy(self.src, dst)
            self.assert_content_equal(self.src, dst)
            self.assert_inode_not_equal(self.src, dst)
        finally:
            os.link = real_link

    @skip_if_not_win32
    def test_windows(self):
        dst = self.dst("target")
        path.link_or_copy(self.src, dst)
        self.assert_content_equal(self.src, dst)

    def test_link_twice(self):
        # Linking the same file twice shouldn't leave duplicates around.
        # See https://github.com/ipython/ipython/issues/6450
        dst = self.dst('target')
        path.link_or_copy(self.src, dst)
        path.link_or_copy(self.src, dst)
        self.assert_inode_equal(self.src, dst)
        nt.assert_equal(sorted(os.listdir(self.tempdir.name)), ['src', 'target'])
