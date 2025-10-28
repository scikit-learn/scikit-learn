# Copyright (c) Amber Brown, 2017
# See LICENSE for details.

import os
import os.path

from pathlib import Path
from subprocess import check_call

from click.testing import CliRunner
from twisted.trial.unittest import TestCase

from towncrier.build import _main as towncrier_build
from towncrier.check import _main as towncrier_check

from .helpers import setup_simple_project, with_isolated_runner, write


def create_project(
    pyproject_path="pyproject.toml", main_branch="main", extra_config=""
):
    """
    Create the project files in the main branch that already has a
    news-fragment and then switch to a new in-work branch.
    """
    setup_simple_project(pyproject_path=pyproject_path, extra_config=extra_config)
    Path("foo/newsfragments/123.feature").write_text("Adds levitation")
    initial_commit(branch=main_branch)
    check_call(["git", "checkout", "-b", "otherbranch"])


def commit(message):
    """Stage and commit the repo in the current working directory

    There must be uncommitted changes otherwise git will complain:
    "nothing to commit, working tree clean"
    """
    check_call(["git", "add", "."])
    check_call(["git", "commit", "-m", message])


def stage():
    """Stage a commit to the repo in the current working directory

    There must be uncommitted changes otherwise git will complain:
    "nothing to commit, working tree clean"
    """
    check_call(["git", "add", "."])


def initial_commit(branch="main"):
    """
    Create a git repo, configure it and make an initial commit

    There must be uncommitted changes otherwise git will complain:
    "nothing to commit, working tree clean"
    """
    # --initial-branch is explicitly set to `main` because
    # git has deprecated the default branch name.
    check_call(["git", "init", f"--initial-branch={branch}"])
    # Without ``git config` user.name and user.email `git commit` fails
    # unless the settings are set globally
    check_call(["git", "config", "user.name", "user"])
    check_call(["git", "config", "user.email", "user@example.com"])
    commit("Initial Commit")


class TestChecker(TestCase):
    maxDiff = None

    def test_git_fails(self):
        """
        If git fails to report a comparison, git's output is reported to aid in
        debugging the situation.
        """
        runner = CliRunner()
        with runner.isolated_filesystem():
            create_project("pyproject.toml")

            result = runner.invoke(towncrier_check, ["--compare-with", "hblaugh"])
            self.assertIn("git produced output while failing", result.output)
            self.assertIn("hblaugh", result.output)

    def test_no_changes_made(self):
        self._test_no_changes_made(
            "pyproject.toml", lambda runner, main, argv: runner.invoke(main, argv)
        )

    def test_no_changes_made_config_path(self):
        pyproject = "not-pyproject.toml"
        self._test_no_changes_made(
            pyproject,
            lambda runner, main, argv: runner.invoke(
                main, argv + ["--config", pyproject]
            ),
        )

    def _test_no_changes_made(self, pyproject_path, invoke):
        """
        When no changes are made on a new branch, no checks are performed.
        """
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project(pyproject_path, main_branch="master")

            result = invoke(runner, towncrier_check, ["--compare-with", "master"])

            self.assertEqual(0, result.exit_code, result.output)
            self.assertEqual(
                "On master branch, or no diffs, so no newsfragment required.\n",
                result.output,
            )

    @with_isolated_runner
    def test_fragment_exists(self, runner):
        create_project("pyproject.toml")

        write("foo/somefile.py", "import os")
        commit("add a file")

        fragment_path = Path("foo/newsfragments/1234.feature").absolute()
        write(fragment_path, "Adds gravity back")
        commit("add a newsfragment")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])

        self.assertTrue(
            result.output.endswith("Found:\n1. " + str(fragment_path) + "\n"),
            (result.output, str(fragment_path)),
        )
        self.assertEqual(0, result.exit_code, result.output)

    @with_isolated_runner
    def test_fragment_exists_hidden(self, runner):
        """
        Location of fragments can be configured using tool.towncrier.directory.
        """
        create_project("pyproject.toml", extra_config="directory = 'deep/fragz'\n")

        write("foo/bar/somefile.py", "import os")
        commit("add a file")

        fragment_path = Path("deep/fragz/1234.feature").absolute()
        write(fragment_path, "Adds gravity back")
        commit("add a newsfragment")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])

        self.assertTrue(
            result.output.endswith("Found:\n1. " + str(fragment_path) + "\n"),
            (result.output, str(fragment_path)),
        )
        self.assertEqual(0, result.exit_code, result.output)

    def test_fragment_missing(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project("pyproject.toml", main_branch="master")

            file_path = "foo/somefile.py"
            with open(file_path, "w") as f:
                f.write("import os")

            check_call(["git", "add", "foo/somefile.py"])
            check_call(["git", "commit", "-m", "add a file"])

            result = runner.invoke(towncrier_check, ["--compare-with", "master"])

            self.assertEqual(1, result.exit_code)
            self.assertTrue(
                result.output.endswith("No new newsfragments found on this branch.\n")
            )

    def test_fragment_exists_but_not_in_check(self):
        """A fragment that exists but is marked as check=False is ignored by the check."""
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project(
                "pyproject.toml",
                main_branch="master",
                extra_config="[[tool.towncrier.type]]\n"
                'directory = "feature"\n'
                'name = "Features"\n'
                "showcontent = true\n"
                "[[tool.towncrier.type]]\n"
                'directory = "sut"\n'
                'name = "System Under Test"\n'
                "showcontent = true\n"
                "check=false\n",
            )

            file_path = "foo/somefile.py"
            write(file_path, "import os")

            fragment_path = Path("foo/newsfragments/1234.sut").absolute()
            write(fragment_path, "Adds gravity back")
            commit("add a file and a newsfragment")

            result = runner.invoke(towncrier_check, ["--compare-with", "master"])

            self.assertEqual(1, result.exit_code)
            self.assertTrue(
                result.output.endswith(
                    "Found newsfragments of unchecked types in the branch:\n1. "
                    + str(fragment_path)
                    + "\n"
                ),
                (result.output, str(fragment_path)),
            )

    def test_fragment_exists_and_staged(self):
        """A fragment exists and is added in staging. Pass only if staging on the command line"""
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project(
                "pyproject.toml",
                main_branch="master",
                extra_config="[[tool.towncrier.type]]\n"
                'directory = "feature"\n'
                'name = "Features"\n'
                "showcontent = true\n"
                "[[tool.towncrier.type]]\n"
                'directory = "sut"\n'
                'name = "System Under Test"\n'
                "showcontent = true\n",
            )

            file_path = "foo/somefile.py"
            write(file_path, "import os")

            commit("add some files for test initialization")

            fragment_path = Path("foo/newsfragments/1234.feature").absolute()
            write(fragment_path, "Adds gravity back")
            stage()

            result = runner.invoke(towncrier_check, ["--compare-with", "master"])

            self.assertEqual(1, result.exit_code)
            result = runner.invoke(
                towncrier_check, ["--staged", "--compare-with", "master"]
            )
            self.assertEqual(0, result.exit_code)

    def test_fragment_exists_and_in_check(self):
        """
        A fragment that exists but is not marked as check=False is
        not ignored by the check, even if other categories are marked as check=False.
        """
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project(
                "pyproject.toml",
                main_branch="master",
                extra_config="[[tool.towncrier.type]]\n"
                'directory = "feature"\n'
                'name = "Features"\n'
                "showcontent = true\n"
                "[[tool.towncrier.type]]\n"
                'directory = "sut"\n'
                'name = "System Under Test"\n'
                "showcontent = true\n"
                "check=false\n",
            )

            file_path = "foo/somefile.py"
            write(file_path, "import os")

            fragment_path = Path("foo/newsfragments/1234.feature").absolute()
            write(fragment_path, "Adds gravity back")
            commit("add a file and a newsfragment")

            result = runner.invoke(towncrier_check, ["--compare-with", "master"])

            self.assertEqual(0, result.exit_code)
            self.assertTrue(
                result.output.endswith("Found:\n1. " + str(fragment_path) + "\n"),
                (result.output, str(fragment_path)),
            )

    def test_none_stdout_encoding_works(self):
        """
        No failure when output is piped causing None encoding for stdout.
        """
        try:
            runner = CliRunner(mix_stderr=False)
        except TypeError:
            # Fallback for older Click versions (or unexpected signature)
            print("TypeError with mix_stderr=False, falling back to echo_stdin=True")
            runner = CliRunner(echo_stdin=True)

        with runner.isolated_filesystem():
            create_project("pyproject.toml", main_branch="master")

            fragment_path = "foo/newsfragments/1234.feature"
            with open(fragment_path, "w") as f:
                f.write("Adds gravity back")

            check_call(["git", "add", fragment_path])
            check_call(["git", "commit", "-m", "add a newsfragment"])

            result = runner.invoke(towncrier_check, ["--compare-with", "master"])

        self.assertEqual(0, result.exit_code)
        self.assertEqual(0, len(result.stderr))

    def test_first_release(self):
        """
        The checks should be skipped on a branch that creates the news file.

        If the checks are not skipped in this case, towncrier check would fail
        for the first release that has a changelog.
        """
        runner = CliRunner()

        with runner.isolated_filesystem():
            # Arrange
            create_project()
            # Before any release, the NEWS file might no exist.
            self.assertNotIn("NEWS.rst", os.listdir("."))

            runner.invoke(towncrier_build, ["--yes", "--version", "1.0"])
            commit("Prepare a release")
            # When missing,
            # the news file is automatically created with a new release.
            self.assertIn("NEWS.rst", os.listdir("."))

            # Act
            result = runner.invoke(towncrier_check, ["--compare-with", "main"])

            # Assert
            self.assertEqual(0, result.exit_code, (result, result.output))
            self.assertIn("Checks SKIPPED: news file changes detected", result.output)

    def test_release_branch(self):
        """
        The checks for missing news fragments are skipped on a branch that
        modifies the news file.
        This is a hint that we are on a release branch
        and at release time is expected no not have news-fragment files.
        """
        runner = CliRunner()

        with runner.isolated_filesystem():
            # Arrange
            create_project()

            # Do a first release without any checks.
            # And merge the release branch back into the main branch.
            runner.invoke(towncrier_build, ["--yes", "--version", "1.0"])
            commit("First release")
            # The news file is now created.
            self.assertIn("NEWS.rst", os.listdir("."))
            check_call(["git", "checkout", "main"])
            check_call(
                ["git", "merge", "otherbranch", "-m", "Sync release in main branch."]
            )

            # We have a new feature branch that has a news fragment that
            # will be merged to the main branch.
            check_call(["git", "checkout", "-b", "new-feature-branch"])
            write("foo/newsfragments/456.feature", "Foo the bar")
            commit("A feature in the second release.")
            check_call(["git", "checkout", "main"])
            check_call(
                [
                    "git",
                    "merge",
                    "new-feature-branch",
                    "-m",
                    "Merge new-feature-branch.",
                ]
            )

            # We now have the new release branch.
            check_call(["git", "checkout", "-b", "next-release"])
            runner.invoke(towncrier_build, ["--yes", "--version", "2.0"])
            commit("Second release")

            # Act
            result = runner.invoke(towncrier_check, ["--compare-with", "main"])

            # Assert
            self.assertEqual(0, result.exit_code, (result, result.output))
            self.assertIn("Checks SKIPPED: news file changes detected", result.output)

    def test_get_default_compare_branch_missing(self):
        """
        If there's no recognized remote origin, exit with an error.
        """
        runner = CliRunner()

        with runner.isolated_filesystem():
            create_project()

            result = runner.invoke(towncrier_check)

        self.assertEqual(1, result.exit_code)
        self.assertEqual("Could not detect default branch. Aborting.\n", result.output)

    @with_isolated_runner
    def test_in_different_dir_with_nondefault_newsfragments_directory(self, runner):
        """
        It can check the fragments located in a sub-directory
        that is specified using the `--dir` CLI argument.
        """
        main_branch = "main"
        Path("pyproject.toml").write_text(
            # Important to customize `config.directory` because the default
            # already supports this scenario.
            "[tool.towncrier]\n"
            + 'directory = "changelog.d"\n'
        )
        subproject1 = Path("foo")
        (subproject1 / "foo").mkdir(parents=True)
        (subproject1 / "foo/__init__.py").write_text("")
        (subproject1 / "changelog.d").mkdir(parents=True)
        (subproject1 / "changelog.d/123.feature").write_text("Adds levitation")
        initial_commit(branch=main_branch)
        check_call(["git", "checkout", "-b", "otherbranch"])

        # We add a code change but forget to add a news fragment.
        write(subproject1 / "foo/somefile.py", "import os")
        commit("add a file")
        result = runner.invoke(
            towncrier_check,
            (
                "--config",
                "pyproject.toml",
                "--dir",
                str(subproject1),
                "--compare-with",
                "main",
            ),
        )

        self.assertEqual(1, result.exit_code)
        self.assertTrue(
            result.output.endswith("No new newsfragments found on this branch.\n")
        )

        # We add the news fragment.
        fragment_path = (subproject1 / "changelog.d/124.feature").absolute()
        write(fragment_path, "Adds gravity back")
        commit("add a newsfragment")
        result = runner.invoke(
            towncrier_check,
            ("--config", "pyproject.toml", "--dir", "foo", "--compare-with", "main"),
        )

        self.assertEqual(0, result.exit_code, result.output)
        self.assertTrue(
            result.output.endswith("Found:\n1. " + str(fragment_path) + "\n"),
            (result.output, str(fragment_path)),
        )

        # We add a change in a different subproject without a news fragment.
        # Checking subproject1 should pass.
        subproject2 = Path("bar")
        (subproject2 / "bar").mkdir(parents=True)
        (subproject2 / "changelog.d").mkdir(parents=True)
        write(subproject2 / "bar/somefile.py", "import os")
        commit("add a file")
        result = runner.invoke(
            towncrier_check,
            (
                "--config",
                "pyproject.toml",
                "--dir",
                subproject1,
                "--compare-with",
                "main",
            ),
        )

        self.assertEqual(0, result.exit_code, result.output)
        self.assertTrue(
            result.output.endswith("Found:\n1. " + str(fragment_path) + "\n"),
            (result.output, str(fragment_path)),
        )

        # Checking subproject2 should result in an error.
        result = runner.invoke(
            towncrier_check,
            (
                "--config",
                "pyproject.toml",
                "--dir",
                subproject2,
                "--compare-with",
                "main",
            ),
        )
        self.assertEqual(1, result.exit_code)
        self.assertTrue(
            result.output.endswith("No new newsfragments found on this branch.\n")
        )

    @with_isolated_runner
    def test_ignored_files(self, runner):
        """
        When `ignore` is set in config, files with those names are ignored.
        Configuration supports wildcard matching with `fnmatch`.
        """
        create_project(
            "pyproject.toml",
            extra_config='ignore = ["template.jinja", "star_wildcard*"]',
        )

        write(
            "foo/newsfragments/124.feature",
            "This fragment has valid name (control case)",
        )
        write("foo/newsfragments/template.jinja", "This is manually ignored")
        write("foo/newsfragments/.gitignore", "gitignore is automatically ignored")
        write("foo/newsfragments/star_wildcard_foo", "Manually ignored with * wildcard")
        commit("add stuff")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])
        self.assertEqual(0, result.exit_code, result.output)

    @with_isolated_runner
    def test_invalid_fragment_name(self, runner):
        """
        Fails if a news fragment has an invalid name, even if `ignore` is not set in
        the config.
        """
        create_project("pyproject.toml")

        write(
            "foo/newsfragments/124.feature",
            "This fragment has valid name (control case)",
        )
        write(
            "foo/newsfragments/feature.125",
            "This has issue and category wrong way round",
        )
        write(
            "NEWS.rst",
            "Modification of news file should not skip check of invalid names",
        )
        commit("add stuff")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])
        self.assertEqual(1, result.exit_code, result.output)
        self.assertIn("Invalid news fragment name: feature.125", result.output)

    @with_isolated_runner
    def test_issue_pattern(self, runner):
        """
        Fails if an issue name goes against the configured pattern.
        """
        create_project(
            "pyproject.toml",
            extra_config='issue_pattern = "\\\\d+"',
        )
        write(
            "foo/newsfragments/123.feature",
            "This fragment has a valid name",
        )
        write(
            "foo/newsfragments/+abcdefg.feature",
            "This fragment has a valid name (orphan fragment)",
        )
        commit("add stuff")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])
        self.assertEqual(0, result.exit_code, result.output)

    @with_isolated_runner
    def test_issue_pattern_invalid_with_suffix(self, runner):
        """
        Fails if an issue name goes against the configured pattern.
        """
        create_project(
            "pyproject.toml",
            extra_config='issue_pattern = "\\\\d+"',
        )
        write(
            "foo/newsfragments/AAA.BBB.feature.md",
            "This fragment has an invalid name (should be digits only)",
        )
        commit("add stuff")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])
        self.assertEqual(1, result.exit_code, result.output)
        self.assertIn(
            "Error: Issue name 'AAA.BBB' does not match the configured pattern, '\\d+'",
            result.output,
        )

    @with_isolated_runner
    def test_issue_pattern_invalid(self, runner):
        """
        Fails if an issue name goes against the configured pattern.
        """
        create_project(
            "pyproject.toml",
            extra_config='issue_pattern = "\\\\d+"',
        )
        write(
            "foo/newsfragments/AAA.BBB.feature",
            "This fragment has an invalid name (should be digits only)",
        )
        commit("add stuff")

        result = runner.invoke(towncrier_check, ["--compare-with", "main"])
        self.assertEqual(1, result.exit_code, result.output)
        self.assertIn(
            "Error: Issue name 'AAA.BBB' does not match the configured pattern, '\\d+'",
            result.output,
        )
