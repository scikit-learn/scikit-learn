# Copyright (c) towncrier contributors, 2025
# See LICENSE for details.

import os.path
import shutil
import unittest

from pathlib import Path
from subprocess import check_call

from click.testing import CliRunner
from twisted.trial.unittest import TestCase

from towncrier import _hg, _vcs

from .helpers import setup_simple_project, write


hg_available = shutil.which("hg") is not None


def create_project(
    pyproject_path="pyproject.toml", main_branch="main", extra_config=""
):
    setup_simple_project(pyproject_path=pyproject_path, extra_config=extra_config)
    Path("foo/newsfragments/123.feature").write_text("Adds levitation")
    initial_commit(branch=main_branch)
    check_call(["hg", "branch", "otherbranch"])


def initial_commit(branch="main"):
    """
    Create a mercurial repo, configure it and make an initial commit

    There must be uncommitted changes
    """
    check_call(["hg", "init", "."])
    check_call(["hg", "branch", branch])
    commit("Initial Commit")


def commit(message):
    """Commit the repo in the current working directory

    There must be uncommitted changes
    """
    check_call(["hg", "addremove", "."])
    check_call(["hg", "commit", "--user", "Example <test@example.com>", "-m", message])


@unittest.skipUnless(hg_available, "requires 'mercurial' to be installed")
class TestHg(TestCase):
    def test_get_default_compare_branch(self):
        """
        Test the 'get_default_compare_branch' behavior.
        """
        assert _hg.get_default_compare_branch(["main", "default"]) == "default"
        assert _hg.get_default_compare_branch(["main", "a_topic"]) is None

    def test_empty_remove(self):
        """
        If remove_files gets an empty list, it returns gracefully.
        """
        _hg.remove_files([])

    def test_complete_scenario(self):
        """
        Tests all the _hg functions that interact with an actual repository.

        Setting up a project is a little slow, hence the grouping of all the
        tests in one.
        """
        runner = CliRunner()
        with runner.isolated_filesystem():
            create_project("pyproject.toml")

            # make sure _vcs._get_mod properly detects a mercurial repo
            self.assertEqual(_hg, _vcs._get_mod("."))

            write("changes/000.misc.rst", "some change")

            _hg.stage_newsfile(".", "changes/000.misc.rst")

            commit("commit 1")

            branches = sorted(_hg.get_remote_branches("."))
            self.assertEqual(["main", "otherbranch"], branches)

            self.assertEqual(
                [os.path.join("changes", "000.misc.rst")],
                _hg.list_changed_files_compared_to_branch(".", "main", False),
            )

            write("changes/001.misc.rst", "some change")
            _hg.remove_files(["changes/000.misc.rst", "changes/001.misc.rst"])
