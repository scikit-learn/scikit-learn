# Copyright (c) Amber Brown, 2015
# See LICENSE for details.


import warnings

from twisted.trial.unittest import TestCase

from towncrier import _git


class TestGit(TestCase):
    def test_empty_remove(self):
        """
        If remove_files gets an empty list, it returns gracefully.
        """
        _git.remove_files([])

    def test_get_default_compare_branch_main(self):
        """
        If there's a remote branch origin/main, prefer it over everything else.
        """
        branch = _git.get_default_compare_branch(["origin/master", "origin/main"])

        self.assertEqual("origin/main", branch)

    def test_get_default_compare_branch_fallback(self):
        """
        If there's origin/master and no main, use it and warn about it.
        """
        with warnings.catch_warnings(record=True) as w:
            branch = _git.get_default_compare_branch(["origin/master", "origin/foo"])

        self.assertEqual("origin/master", branch)
        self.assertTrue(w[0].message.args[0].startswith('Using "origin/master'))
