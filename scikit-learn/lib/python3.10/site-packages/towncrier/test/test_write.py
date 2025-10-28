# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

import os

from pathlib import Path
from textwrap import dedent

from click.testing import CliRunner
from twisted.trial.unittest import TestCase

from .._builder import render_fragments, split_fragments
from .._writer import append_to_newsfile
from ..build import _main
from .helpers import read_pkg_resource, write


class WritingTests(TestCase):
    maxDiff = None

    def test_append_at_top(self):
        fragments = {
            "": {
                ("142", "misc", 0): "",
                ("1", "misc", 0): "",
                ("4", "feature", 0): "Stuff!",
                ("4", "feature", 1): "Second Stuff!",
                ("2", "feature", 0): "Foo added.",
                ("72", "feature", 0): "Foo added.",
            },
            "Names": {},
            "Web": {("3", "bugfix", 0): "Web fixed."},
        }

        definitions = {
            "feature": {"name": "Features", "showcontent": True},
            "bugfix": {"name": "Bugfixes", "showcontent": True},
            "misc": {"name": "Misc", "showcontent": False},
        }

        expected_output = """MyProject 1.0 (never)
=====================

Features
--------

- Foo added. (#2, #72)
- Stuff! (#4)
- Second Stuff! (#4)


Misc
----

- #1, #142


Names
-----

No significant changes.


Web
---

Bugfixes
~~~~~~~~

- Web fixed. (#3)


Old text.
"""

        tempdir = self.mktemp()
        os.makedirs(tempdir)

        with open(os.path.join(tempdir, "NEWS.rst"), "w") as f:
            f.write("Old text.\n")

        fragments = split_fragments(fragments, definitions)

        template = read_pkg_resource("templates/default.rst")

        append_to_newsfile(
            tempdir,
            "NEWS.rst",
            ".. towncrier release notes start\n",
            "",
            render_fragments(
                template,
                None,
                fragments,
                definitions,
                ["-", "~"],
                wrap=True,
                versiondata={"name": "MyProject", "version": "1.0", "date": "never"},
            ),
            single_file=True,
        )

        with open(os.path.join(tempdir, "NEWS.rst")) as f:
            output = f.read()

        self.assertEqual(expected_output, output)

    def test_append_at_top_with_hint(self):
        """
        If there is a comment with C{.. towncrier release notes start},
        towncrier will add the version notes after it.
        """
        fragments = {
            "": {
                ("142", "misc", 0): "",
                ("1", "misc", 0): "",
                ("4", "feature", 0): "Stuff!",
                ("2", "feature", 0): "Foo added.",
                ("72", "feature", 0): "Foo added.",
                ("99", "feature", 0): "Foo! " * 100,
            },
            "Names": {},
            "Web": {("3", "bugfix", 0): "Web fixed."},
        }

        definitions = {
            "feature": {"name": "Features", "showcontent": True},
            "bugfix": {"name": "Bugfixes", "showcontent": True},
            "misc": {"name": "Misc", "showcontent": False},
        }

        expected_output = """Hello there! Here is some info.

.. towncrier release notes start

MyProject 1.0 (never)
=====================

Features
--------

- Foo added. (#2, #72)
- Stuff! (#4)
- Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo!
  Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! Foo! (#99)


Misc
----

- #1, #142


Names
-----

No significant changes.


Web
---

Bugfixes
~~~~~~~~

- Web fixed. (#3)


Old text.
"""

        tempdir = self.mktemp()
        write(
            os.path.join(tempdir, "NEWS.rst"),
            contents="""\
                Hello there! Here is some info.

                .. towncrier release notes start
                Old text.
            """,
            dedent=True,
        )

        fragments = split_fragments(fragments, definitions)

        template = read_pkg_resource("templates/default.rst")

        append_to_newsfile(
            tempdir,
            "NEWS.rst",
            ".. towncrier release notes start\n",
            "",
            render_fragments(
                template,
                None,
                fragments,
                definitions,
                ["-", "~"],
                wrap=True,
                versiondata={"name": "MyProject", "version": "1.0", "date": "never"},
            ),
            single_file=True,
        )

        with open(os.path.join(tempdir, "NEWS.rst")) as f:
            output = f.read()

        self.assertEqual(expected_output, output)

    def test_multiple_file_no_start_string(self):
        """
        When no `start_string` is defined, the generated content is added at
        the start of the file.
        """
        tempdir = self.mktemp()
        os.makedirs(tempdir)

        definitions = {}
        fragments = split_fragments(fragments={}, definitions=definitions)

        template = read_pkg_resource("templates/default.rst")

        content = render_fragments(
            template=template,
            issue_format=None,
            fragments=fragments,
            definitions=definitions,
            underlines=["-", "~"],
            wrap=True,
            versiondata={"name": "MyProject", "version": "1.0", "date": "never"},
        )

        append_to_newsfile(
            directory=tempdir,
            filename="NEWS.rst",
            start_string=None,
            top_line="",
            content=content,
            single_file=True,
        )

        with open(os.path.join(tempdir, "NEWS.rst")) as f:
            output = f.read()

        expected_output = dedent(
            """\
            MyProject 1.0 (never)
            =====================
        """
        )

        self.assertEqual(expected_output, output)

    def test_with_title_format_duplicate_version_raise(self):
        """
        When `single_file` enabled as default,
        and fragments of `version` already produced in the newsfile,
        a duplicate `build` will throw a ValueError.
        """
        runner = CliRunner()

        def do_build_once():
            with open("newsfragments/123.feature", "w") as f:
                f.write("Adds levitation")

            result = runner.invoke(
                _main,
                [
                    "--version",
                    "7.8.9",
                    "--name",
                    "foo",
                    "--date",
                    "01-01-2001",
                    "--yes",
                ],
            )
            return result

        # `single_file` default as true
        with runner.isolated_filesystem():
            with open("pyproject.toml", "w") as f:
                f.write(
                    dedent(
                        """
                    [tool.towncrier]
                    title_format="{name} {version} ({project_date})"
                    filename="{version}-notes.rst"
                    """
                    ).lstrip()
                )
            with open("{version}-notes.rst", "w") as f:
                f.write("Release Notes\n\n.. towncrier release notes start\n")
            os.mkdir("newsfragments")

            result = do_build_once()
            self.assertEqual(0, result.exit_code)
            # build again with the same version
            result = do_build_once()
            self.assertNotEqual(0, result.exit_code)
            self.assertIsInstance(result.exception, ValueError)
            self.assertSubstring(
                "already produced newsfiles for this version", result.exception.args[0]
            )

    def test_single_file_false_overwrite_duplicate_version(self):
        """
        When `single_file` disabled, multiple newsfiles generated and
        the content of which get overwritten each time.
        """
        runner = CliRunner()

        def do_build_once():
            with open("newsfragments/123.feature", "w") as f:
                f.write("Adds levitation")

            result = runner.invoke(
                _main,
                [
                    "--version",
                    "7.8.9",
                    "--name",
                    "foo",
                    "--date",
                    "01-01-2001",
                    "--yes",
                ],
            )
            return result

        # single_file = false
        with runner.isolated_filesystem():
            with open("pyproject.toml", "w") as f:
                f.write(
                    dedent(
                        """
                    [tool.towncrier]
                    single_file=false
                    title_format="{name} {version} ({project_date})"
                    filename="{version}-notes.rst"
                    """
                    ).lstrip()
                )
            os.mkdir("newsfragments")

            result = do_build_once()
            self.assertEqual(0, result.exit_code)
            # build again with the same version
            result = do_build_once()
            self.assertEqual(0, result.exit_code)

            notes = list(Path.cwd().glob("*-notes.rst"))
            self.assertEqual(1, len(notes))
            self.assertEqual("7.8.9-notes.rst", notes[0].name)

            with open(notes[0]) as f:
                output = f.read()

        expected_output = dedent(
            """\
            foo 7.8.9 (01-01-2001)
            ======================

            Features
            --------

            - Adds levitation (#123)
            """
        )

        self.assertEqual(expected_output, output)
