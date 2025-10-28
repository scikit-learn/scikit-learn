# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

import os
import string

from pathlib import Path
from textwrap import dedent
from unittest import mock

from click.testing import CliRunner
from twisted.trial.unittest import TestCase

from ..create import DEFAULT_CONTENT, _main
from .helpers import setup_simple_project, with_isolated_runner


class TestCli(TestCase):
    maxDiff = None

    def _test_success(
        self,
        content=None,
        config=None,
        mkdir=True,
        additional_args=None,
        eof_newline=True,
    ):
        runner = CliRunner()

        with runner.isolated_filesystem():
            setup_simple_project(config=config, mkdir_newsfragments=mkdir)

            args = ["123.feature.rst"]
            if content is None:
                content = [DEFAULT_CONTENT]
            if additional_args is not None:
                args.extend(additional_args)
            result = runner.invoke(_main, args)

            self.assertEqual(["123.feature.rst"], os.listdir("foo/newsfragments"))

            if eof_newline:
                content.append("")
            with open("foo/newsfragments/123.feature.rst") as fh:
                self.assertEqual("\n".join(content), fh.read())

        self.assertEqual(0, result.exit_code)

    def test_basics(self):
        """Ensure file created where output directory already exists."""
        self._test_success(mkdir=True)

    def test_directory_created(self):
        """Ensure both file and output directory created if necessary."""
        self._test_success(mkdir=False)

    def test_edit_without_comments(self):
        """Create file with dynamic content."""
        content = ["This is line 1", "This is line 2"]
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "\n".join(content)
            self._test_success(content=content, additional_args=["--edit"])
            mock_edit.assert_called_once_with(
                "\n# Please write your news content. Lines starting "
                "with '#' will be ignored, and\n# an empty message aborts.\n",
                extension=".rst",
            )

    def test_edit_with_comment(self):
        """Create file editly with ignored line."""
        content = ["This is line 1", "This is line 2"]
        comment = "# I am ignored"
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "\n".join(content[:1] + [comment] + content[1:])
            self._test_success(content=content, additional_args=["--edit"])

    def test_edit_abort(self):
        """Create file editly and abort."""
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = None

            runner = CliRunner()

            with runner.isolated_filesystem():
                setup_simple_project(config=None, mkdir_newsfragments=True)
                result = runner.invoke(_main, ["123.feature.rst", "--edit"])
                self.assertEqual([], os.listdir("foo/newsfragments"))
                self.assertEqual(1, result.exit_code)

    def test_edit_markdown_extension(self):
        """
        The temporary file extension used when editing is ``.md`` if the main filename
        also uses that extension.
        """

        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "This is line 1"
            self._test_success(
                content=["This is line 1"],
                config=dedent(
                    """\
                    [tool.towncrier]
                    package = "foo"
                    filename = "README.md"
                    """
                ),
                additional_args=["--edit"],
            )
            mock_edit.assert_called_once_with(
                "\n# Please write your news content. Lines starting "
                "with '#' will be ignored, and\n# an empty message aborts.\n",
                extension=".md",
            )

    def test_edit_unknown_extension(self):
        """
        The temporary file extension used when editing is ``.txt`` if it the main
        filename isn't ``.rst`` or ``.md``.
        """

        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "This is line 1"
            self._test_success(
                content=["This is line 1"],
                config=dedent(
                    """\
                    [tool.towncrier]
                    package = "foo"
                    filename = "README.FIRST"
                    """
                ),
                additional_args=["--edit"],
            )
            mock_edit.assert_called_once_with(
                "\n# Please write your news content. Lines starting "
                "with '#' will be ignored, and\n# an empty message aborts.\n",
                extension=".txt",
            )

    def test_content(self):
        """
        When creating a new fragment the content can be passed as a
        command line argument.
        The text editor is not invoked.
        """
        content_line = "This is a content"
        self._test_success(content=[content_line], additional_args=["-c", content_line])

    def test_content_without_eof_newline(self):
        """
        When creating a new fragment the content can be passed as a command line
        argument. The text editor is not invoked, and no eof newline is added if the
        config option is set.
        """
        config = dedent(
            """\
            [tool.towncrier]
            package = "foo"
            create_eof_newline = false
            """
        )
        content_line = "This is a content"
        self._test_success(
            content=[content_line],
            additional_args=["-c", content_line],
            config=config,
            eof_newline=False,
        )

    def test_message_and_edit(self):
        """
        When creating a new message, a initial content can be passed via
        the command line and continue modifying the content by invoking the
        text editor.
        """
        content_line = "This is a content line"
        edit_content = ["This is line 1", "This is line 2"]
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "\n".join(edit_content)
            self._test_success(
                content=edit_content, additional_args=["-c", content_line, "--edit"]
            )
            mock_edit.assert_called_once_with(
                f"{content_line}\n\n# Please write your news content. Lines starting "
                "with '#' will be ignored, and\n# an empty message aborts.\n",
                extension=".rst",
            )

    def test_different_directory(self):
        """Ensure non-standard directories are used."""
        runner = CliRunner()
        config = dedent(
            """\
            [tool.towncrier]
            directory = "releasenotes"
            """
        )

        with runner.isolated_filesystem():
            setup_simple_project(config=config, mkdir_newsfragments=False)
            os.mkdir("releasenotes")

            result = runner.invoke(_main, ["123.feature.rst"])

            self.assertEqual(["123.feature.rst"], os.listdir("releasenotes"))

        self.assertEqual(0, result.exit_code)

    def test_invalid_section(self):
        """Ensure creating a path without a valid section is rejected."""
        runner = CliRunner()

        with runner.isolated_filesystem():
            setup_simple_project()

            self.assertEqual([], os.listdir("foo/newsfragments"))

            result = runner.invoke(_main, ["123.foobar.rst"])

            self.assertEqual([], os.listdir("foo/newsfragments"))

        self.assertEqual(type(result.exception), SystemExit, result.exception)
        self.assertIn(
            "Expected filename '123.foobar.rst' to be of format", result.output
        )

    @with_isolated_runner
    def test_custom_extension(self, runner: CliRunner):
        """Ensure we can still create fragments with custom extensions."""
        setup_simple_project()
        frag_path = Path("foo", "newsfragments")

        result = runner.invoke(_main, ["123.feature.txt"])
        self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        # No '.rst' extension added.
        self.assertEqual(fragments, ["123.feature.txt"])

    @with_isolated_runner
    def test_md_filename_extension(self, runner: CliRunner):
        """Ensure changelog filename extension is used if .md"""
        setup_simple_project(extra_config='filename = "changes.md"')
        frag_path = Path("foo", "newsfragments")

        result = runner.invoke(_main, ["123.feature"])
        self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        # No '.rst' extension added.
        self.assertEqual(fragments, ["123.feature.md"])

    @with_isolated_runner
    def test_no_filename_extension(self, runner: CliRunner):
        """
        When the NEWS filename has no extension, new fragments are will not have an
        extension added.
        """
        # The name of the file where towncrier will generate
        # the final release notes is named `RELEASE_NOTES`
        # for this test (with no file extension).
        setup_simple_project(extra_config='filename = "RELEASE_NOTES"')
        frag_path = Path("foo", "newsfragments")

        result = runner.invoke(_main, ["123.feature"])
        self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        # No '.rst' extension added.
        self.assertEqual(fragments, ["123.feature"])

    @with_isolated_runner
    def test_file_exists(self, runner: CliRunner):
        """Ensure we don't overwrite existing files."""
        setup_simple_project()
        frag_path = Path("foo", "newsfragments")

        for _ in range(3):
            result = runner.invoke(_main, ["123.feature"])
            self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(
            sorted(fragments),
            [
                "123.feature.1.rst",
                "123.feature.2.rst",
                "123.feature.rst",
            ],
        )

    @with_isolated_runner
    def test_file_exists_no_ext(self, runner: CliRunner):
        """
        Ensure we don't overwrite existing files with when not adding filename
        extensions.
        """

        setup_simple_project(extra_config="create_add_extension = false")
        frag_path = Path("foo", "newsfragments")

        for _ in range(3):
            result = runner.invoke(_main, ["123.feature"])
            self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(
            sorted(fragments),
            [
                "123.feature",
                "123.feature.1",
                "123.feature.2",
            ],
        )

    @with_isolated_runner
    def test_file_exists_with_ext(self, runner: CliRunner):
        """
        Ensure we don't overwrite existing files when using an extension after the
        fragment type.
        """
        setup_simple_project()
        frag_path = Path("foo", "newsfragments")

        for _ in range(3):
            result = runner.invoke(_main, ["123.feature.rst"])
            self.assertEqual(result.exit_code, 0, result.output)

        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(
            sorted(fragments),
            [
                "123.feature.1.rst",
                "123.feature.2.rst",
                "123.feature.rst",
            ],
        )

    @with_isolated_runner
    def test_without_filename(self, runner: CliRunner):
        """
        When no filename is provided, the user is prompted for one.
        """
        setup_simple_project()

        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "Edited content"
            result = runner.invoke(_main, input="123\nfeature\n")
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_called_once()
        expected = os.path.join(os.getcwd(), "foo", "newsfragments", "123.feature.rst")
        self.assertEqual(
            result.output,
            f"""Issue number (`+` if none): 123
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}
""",
        )
        with open(expected) as f:
            self.assertEqual(f.read(), "Edited content\n")

    @with_isolated_runner
    def test_without_filename_orphan(self, runner: CliRunner):
        """
        The user can create an orphan fragment from the interactive prompt.
        """
        setup_simple_project()

        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "Orphan content"
            result = runner.invoke(_main, input="+\nfeature\n")
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_called_once()
        expected = os.path.join(os.getcwd(), "foo", "newsfragments", "+")
        self.assertTrue(
            result.output.startswith(
                f"""Issue number (`+` if none): +
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}"""
            ),
            result.output,
        )
        # Check that the file was created with a random name
        created_line = result.output.strip().rsplit("\n", 1)[-1]
        # Get file names in the newsfragments directory.
        files = os.listdir(os.path.join(os.getcwd(), "foo", "newsfragments"))
        # Check that the file name is in the created line.
        created_fragment = created_line.split(" ")[-1]
        self.assertIn(Path(created_fragment).name, files)
        with open(created_fragment) as f:
            self.assertEqual(f.read(), "Orphan content\n")

    @with_isolated_runner
    def test_without_filename_no_orphan_config(self, runner: CliRunner):
        """
        If an empty orphan prefix is set, orphan creation is turned off from interactive
        prompt.
        """
        setup_simple_project(extra_config='orphan_prefix = ""')

        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "Edited content"
            result = runner.invoke(_main, input="+\nfeature\n")
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_called_once()
        expected = os.path.join(os.getcwd(), "foo", "newsfragments", "+.feature.rst")
        self.assertEqual(
            result.output,
            f"""Issue number: +
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}
""",
        )
        with open(expected) as f:
            self.assertEqual(f.read(), "Edited content\n")

    @with_isolated_runner
    def test_sections(self, runner: CliRunner):
        """
        When creating a new fragment, the user can specify the section from the command
        line (and if none is provided, the default section will be used).

        The default section is either the section with a blank path, or else the first
        section defined in the configuration file.
        """
        setup_simple_project(
            extra_config="""
[[tool.towncrier.section]]
name = "Backend"
path = "backend"
[[tool.towncrier.section]]
name = "Frontend"
path = ""
"""
        )
        result = runner.invoke(_main, ["123.feature.rst"])
        self.assertFalse(result.exception, result.output)
        frag_path = Path("foo", "newsfragments")
        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(fragments, ["123.feature.rst"])

        result = runner.invoke(_main, ["123.feature.rst", "--section", "invalid"])
        self.assertTrue(result.exception, result.output)
        self.assertIn(
            "Invalid value for '--section': expected one of 'Backend', 'Frontend'",
            result.output,
        )

        result = runner.invoke(_main, ["123.feature.rst", "--section", "Backend"])
        self.assertFalse(result.exception, result.output)
        frag_path = Path("foo", "backend", "newsfragments")
        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(fragments, ["123.feature.rst"])

    @with_isolated_runner
    def test_sections_without_filename(self, runner: CliRunner):
        """
        When multiple sections exist when the interactive prompt is used, the user is
        prompted to select a section.
        """
        setup_simple_project(
            extra_config="""
[[tool.towncrier.section]]
name = "Backend"
path = ""

[[tool.towncrier.section]]
name = "Frontend"
path = "frontend"
"""
        )
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "Edited content"
            result = runner.invoke(_main, input="2\n123\nfeature\n")
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_called_once()
        expected = os.path.join(
            os.getcwd(), "foo", "frontend", "newsfragments", "123.feature.rst"
        )

        self.assertEqual(
            result.output,
            f"""\
Pick a section:
 1: Backend
 2: Frontend
Section (1, 2) [1]: 2
Issue number (`+` if none): 123
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}
""",
        )

    @with_isolated_runner
    def test_sections_without_filename_with_section_option(self, runner: CliRunner):
        """
        When multiple sections exist and the section is provided via the command line,
        the user isn't prompted to select a section.
        """
        setup_simple_project(
            extra_config="""
[[tool.towncrier.section]]
name = "Backend"
path = ""

[[tool.towncrier.section]]
name = "Frontend"
path = "frontend"
"""
        )
        with mock.patch("click.edit") as mock_edit:
            mock_edit.return_value = "Edited content"
            result = runner.invoke(
                _main, ["--section", "Frontend"], input="123\nfeature\n"
            )
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_called_once()
        expected = os.path.join(
            os.getcwd(), "foo", "frontend", "newsfragments", "123.feature.rst"
        )

        self.assertEqual(
            result.output,
            f"""\
Issue number (`+` if none): 123
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}
""",
        )

    @with_isolated_runner
    def test_sections_all_with_paths(self, runner: CliRunner):
        """
        When all sections have paths, the first is the default.
        """
        setup_simple_project(
            extra_config="""
[[tool.towncrier.section]]
name = "Frontend"
path = "frontend"

[[tool.towncrier.section]]
name = "Backend"
path = "backend"
"""
        )
        result = runner.invoke(_main, ["123.feature.rst"])
        self.assertFalse(result.exception, result.output)
        frag_path = Path("foo", "frontend", "newsfragments")
        fragments = [f.name for f in frag_path.iterdir()]
        self.assertEqual(fragments, ["123.feature.rst"])

    @with_isolated_runner
    def test_without_filename_with_message(self, runner: CliRunner):
        """
        When no filename is provided, the user is prompted for one. If a message is
        provided, the editor isn't opened and the message is used.
        """
        setup_simple_project()

        with mock.patch("click.edit") as mock_edit:
            result = runner.invoke(_main, ["-c", "Fixed this"], input="123\nfeature\n")
            self.assertFalse(result.exception, result.output)
            mock_edit.assert_not_called()
        expected = os.path.join(os.getcwd(), "foo", "newsfragments", "123.feature.rst")
        self.assertEqual(
            result.output,
            f"""Issue number (`+` if none): 123
Fragment type (feature, bugfix, doc, removal, misc): feature
Created news fragment at {expected}
""",
        )
        with open(expected) as f:
            self.assertEqual(f.read(), "Fixed this\n")

    @with_isolated_runner
    def test_create_orphan_fragment(self, runner: CliRunner):
        """
        When a fragment starts with the only the orphan prefix (``+`` by default), the
        create CLI automatically extends the new file's base name to contain a random
        value to avoid commit collisions.
        """
        setup_simple_project()

        frag_path = Path("foo", "newsfragments")
        sub_frag_path = frag_path / "subsection"
        sub_frag_path.mkdir()

        result = runner.invoke(_main, ["+.feature"])
        self.assertEqual(0, result.exit_code)
        result = runner.invoke(
            _main, [str(Path("subsection", "+.feature"))], catch_exceptions=False
        )
        self.assertEqual(0, result.exit_code, result.output)

        fragments = [p for p in frag_path.rglob("*") if p.is_file()]
        self.assertEqual(2, len(fragments))
        change1, change2 = fragments

        self.assertEqual(change1.suffix, ".rst")
        self.assertTrue(change1.stem.startswith("+"))
        self.assertTrue(change1.stem.endswith(".feature"))
        # Length should be '+' character, 8 random hex characters, and ".feature".
        self.assertEqual(len(change1.stem), 1 + 8 + len(".feature"))

        self.assertEqual(change2.suffix, ".rst")
        self.assertTrue(change2.stem.startswith("+"))
        self.assertTrue(change2.stem.endswith(".feature"))
        self.assertEqual(change2.parent, sub_frag_path)
        # Length should be '+' character, 8 random hex characters, and ".feature".
        self.assertEqual(len(change2.stem), 1 + 8 + len(".feature"))

    @with_isolated_runner
    def test_create_orphan_fragment_custom_prefix(self, runner: CliRunner):
        """
        Check that the orphan prefix can be customized.
        """
        setup_simple_project(extra_config='orphan_prefix = "$$$"')

        frag_path = Path("foo", "newsfragments")

        result = runner.invoke(_main, ["$$$.feature"])
        self.assertEqual(0, result.exit_code, result.output)

        fragments = list(frag_path.rglob("*"))
        self.assertEqual(len(fragments), 1)
        change = fragments[0]
        self.assertTrue(change.stem.startswith("$$$"))
        # Length should be '$$$' characters, 8 random hex characters, and ".feature".
        self.assertEqual(len(change.stem), 3 + 8 + len(".feature"))
        # Check the remainder are all hex characters.
        self.assertTrue(
            all(c in string.hexdigits for c in change.stem[3 : -len(".feature")])
        )

    @with_isolated_runner
    def test_in_different_dir_with_nondefault_newsfragments_directory(self, runner):
        """
        When the `--dir` CLI argument is passed,
        it will create a new file in directory that is
        created by combining the `--dir` value
        with the `directory` option from the configuration
        file.
        """
        Path("pyproject.toml").write_text(
            # Important to customize `config.directory` because the default
            # already supports this scenario.
            "[tool.towncrier]\n"
            + 'directory = "changelog.d"\n'
        )
        Path("foo/foo").mkdir(parents=True)
        Path("foo/foo/__init__.py").write_text("")

        result = runner.invoke(
            _main,
            (
                "--config",
                "pyproject.toml",
                "--dir",
                "foo",
                "--content",
                "Adds levitation.",
                "123.feature",
            ),
        )

        self.assertEqual(0, result.exit_code)
        self.assertTrue(Path("foo/changelog.d/123.feature.rst").exists())
