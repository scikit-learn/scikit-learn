# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

import os

from textwrap import dedent

from click.testing import CliRunner
from twisted.trial.unittest import TestCase

from .._settings import ConfigError, load_config
from .._shell import cli
from .helpers import with_isolated_runner, write


class TomlSettingsTests(TestCase):
    def mktemp_project(
        self, *, pyproject_toml: str = "", towncrier_toml: str = ""
    ) -> str:
        """
        Create a temporary directory with a pyproject.toml file in it.
        """
        project_dir = self.mktemp()
        os.makedirs(project_dir)

        if pyproject_toml:
            write(
                os.path.join(project_dir, "pyproject.toml"),
                pyproject_toml,
                dedent=True,
            )

        if towncrier_toml:
            write(
                os.path.join(project_dir, "towncrier.toml"),
                towncrier_toml,
                dedent=True,
            )

        return project_dir

    def test_base(self):
        """
        Test a "base config".
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                orphan_prefix = "~"
            """
        )

        config = load_config(project_dir)
        self.assertEqual(config.package, "foobar")
        self.assertEqual(config.package_dir, ".")
        self.assertEqual(config.filename, "NEWS.rst")
        self.assertEqual(config.underlines, ("=", "-", "~"))
        self.assertEqual(config.orphan_prefix, "~")

    def test_markdown(self):
        """
        If the filename references an .md file and the builtin template doesn't have an
        extension, add .md rather than .rst.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                filename = "NEWS.md"
            """
        )

        config = load_config(project_dir)

        self.assertEqual(config.filename, "NEWS.md")

        self.assertEqual(config.template, ("towncrier.templates", "default.md"))

    def test_explicit_template_extension(self):
        """
        If the filename references an .md file and the builtin template has an
        extension, don't change it.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                filename = "NEWS.md"
                template = "towncrier:default.rst"
            """
        )

        config = load_config(project_dir)

        self.assertEqual(config.filename, "NEWS.md")
        self.assertEqual(config.template, ("towncrier.templates", "default.rst"))

    def test_template_extended(self):
        """
        The template can be any package and resource, and although we look for a
        resource's 'templates' package, it could also be in the specified resource
        directly.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                template = "towncrier.templates:default.rst"
            """
        )

        config = load_config(project_dir)

        self.assertEqual(config.template, ("towncrier.templates", "default.rst"))

    def test_incorrect_single_file(self):
        """
        single_file must be a bool.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                single_file = "a"
            """
        )

        with self.assertRaises(ConfigError) as e:
            load_config(project_dir)

        self.assertEqual(e.exception.failing_option, "single_file")

    def test_incorrect_all_bullets(self):
        """
        all_bullets must be a bool.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                all_bullets = "a"
            """
        )

        with self.assertRaises(ConfigError) as e:
            load_config(project_dir)

        self.assertEqual(e.exception.failing_option, "all_bullets")

    def test_mistype_singlefile(self):
        """
        singlefile is not accepted, single_file is.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                singlefile = "a"
            """
        )

        with self.assertRaises(ConfigError) as e:
            load_config(project_dir)

        self.assertEqual(e.exception.failing_option, "singlefile")

    def test_towncrier_toml_preferred(self):
        """
        Towncrier prefers the towncrier.toml for autodetect over pyproject.toml.
        """
        project_dir = self.mktemp_project(
            towncrier_toml="""
                [tool.towncrier]
                package = "a"
            """,
            pyproject_toml="""
                [tool.towncrier]
                package = "b"
            """,
        )

        config = load_config(project_dir)
        self.assertEqual(config.package, "a")

    def test_pyproject_only_pyproject_toml(self):
        """
        Towncrier will fallback to the [project.name] value in pyproject.toml.

        This tests asserts that the minimal configuration is to do *nothing*
        when using a pyproject.toml file.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [project]
                name = "a"
            """,
        )

        config = load_config(project_dir)
        self.assertEqual(config.package, "a")
        self.assertEqual(config.name, "a")

    def test_pyproject_assert_fallback(self):
        """
        This test is an extensive test of the fallback scenarios
        for the `package` and `name` keys in the towncrier section.

        It will fallback to pyproject.toml:name in any case.
        And as such it checks the various fallback mechanisms
        if the fields are not present in the towncrier.toml, nor
        in the pyproject.toml files.

        This both tests when things are *only* in the pyproject.toml
        and default usage of the data in the towncrier.toml file.
        """
        pyproject_toml = dedent(
            """
        [project]
        name = "foo"
        [tool.towncrier]
        """
        )
        towncrier_toml = dedent(
            """
        [tool.towncrier]
        """
        )
        tests = [
            "",
            "name = '{name}'",
            "package = '{package}'",
            "name = '{name}'",
            "package = '{package}'",
        ]

        def factory(name, package):
            def func(test):
                return dedent(test).format(name=name, package=package)

            return func

        for pp_fields in map(factory(name="a", package="b"), tests):
            pp_toml = pyproject_toml + pp_fields
            for tc_fields in map(factory(name="c", package="d"), tests):
                tc_toml = towncrier_toml + tc_fields

                # Create the temporary project
                project_dir = self.mktemp_project(
                    pyproject_toml=pp_toml,
                    towncrier_toml=tc_toml,
                )

                # Read the configuration file.
                config = load_config(project_dir)

                # Now the values depend on where the fallback
                # is.
                # If something is in towncrier.toml, it will be preferred
                # name fallsback to package
                if "package" in tc_fields:
                    package = "d"
                else:
                    package = "foo"
                self.assertEqual(config.package, package)

                if "name" in tc_fields:
                    self.assertEqual(config.name, "c")
                else:
                    # fall-back to package name
                    self.assertEqual(config.name, package)

    @with_isolated_runner
    def test_load_no_config(self, runner: CliRunner):
        """
        Calling the root CLI without an existing configuration file in the base directory,
        will exit with code 1 and an informative message is sent to standard output.
        """
        temp = self.mktemp()
        os.makedirs(temp)

        result = runner.invoke(cli, ("--dir", temp))

        self.assertEqual(
            result.output,
            f"No configuration file found.\nLooked back from: {os.path.abspath(temp)}\n",
        )
        self.assertEqual(result.exit_code, 1)

    @with_isolated_runner
    def test_load_explicit_missing_config(self, runner: CliRunner):
        """
        Calling the CLI with an incorrect explicit configuration file will exit with
        code 1 and an informative message is sent to standard output.
        """
        config = "not-there.toml"
        result = runner.invoke(cli, ("--config", config))

        self.assertEqual(result.exit_code, 1)
        self.assertEqual(
            result.output,
            f"Configuration file '{os.path.abspath(config)}' not found.\n",
        )

    def test_missing_template(self):
        """
        Towncrier will raise an exception saying when it can't find a template.
        """
        project_dir = self.mktemp_project(
            towncrier_toml="""
                [tool.towncrier]
                template = "foo.rst"
            """
        )

        with self.assertRaises(ConfigError) as e:
            load_config(project_dir)

        self.assertEqual(
            str(e.exception),
            "The template file '{}' does not exist.".format(
                os.path.normpath(os.path.join(project_dir, "foo.rst")),
            ),
        )

    def test_missing_template_in_towncrier(self):
        """
        Towncrier will raise an exception saying when it can't find a template
        from the Towncrier templates.
        """
        project_dir = self.mktemp_project(
            towncrier_toml="""
                [tool.towncrier]
                template = "towncrier:foo"
            """
        )

        with self.assertRaises(ConfigError) as e:
            load_config(project_dir)

        self.assertEqual(
            str(e.exception), "'towncrier' does not have a template named 'foo.rst'."
        )

    def test_custom_types_as_tables_array_deprecated(self):
        """
        Custom fragment categories can be defined inside
        the toml config file using an array of tables
        (a table name in double brackets).

        This functionality is considered deprecated, but we continue
        to support it to keep backward compatibility.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                [[tool.towncrier.type]]
                directory="foo"
                name="Foo"
                showcontent=false

                [[tool.towncrier.type]]
                directory="spam"
                name="Spam"
                showcontent=true

                [[tool.towncrier.type]]
                directory="auto"
                name="Automatic"
                showcontent=true
                check=false
            """
        )
        config = load_config(project_dir)
        expected = [
            (
                "foo",
                {
                    "name": "Foo",
                    "showcontent": False,
                    "check": True,
                },
            ),
            (
                "spam",
                {
                    "name": "Spam",
                    "showcontent": True,
                    "check": True,
                },
            ),
            (
                "auto",
                {
                    "name": "Automatic",
                    "showcontent": True,
                    "check": False,
                },
            ),
        ]
        expected = dict(expected)
        actual = config.types
        self.assertDictEqual(expected, actual)

    def test_custom_types_as_tables(self):
        """
        Custom fragment categories can be defined inside
        the toml config file using tables.
        """
        project_dir = self.mktemp_project(
            pyproject_toml="""
                [tool.towncrier]
                package = "foobar"
                [tool.towncrier.fragment.feat]
                ignored_field="Bazz"
                [tool.towncrier.fragment.fix]
                [tool.towncrier.fragment.chore]
                name = "Other Tasks"
                showcontent = false
                [tool.towncrier.fragment.auto]
                name = "Automatic"
                check = false
            """
        )
        config = load_config(project_dir)
        expected = {
            "chore": {
                "name": "Other Tasks",
                "showcontent": False,
                "check": True,
            },
            "feat": {
                "name": "Feat",
                "showcontent": True,
                "check": True,
            },
            "fix": {
                "name": "Fix",
                "showcontent": True,
                "check": True,
            },
            "auto": {
                "name": "Automatic",
                "showcontent": True,
                "check": False,
            },
        }
        actual = config.types
        self.assertDictEqual(expected, actual)
