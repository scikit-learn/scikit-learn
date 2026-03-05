"""Tests for the minimum dependencies in README.rst and pyproject.toml"""

import os
import re
import tomllib
from collections import defaultdict
from pathlib import Path

import pytest

import sklearn
from sklearn._min_dependencies import dependent_packages
from sklearn.utils.fixes import parse_version

# minimal dependencies and pyproject definitions for testing the pyproject tests

TOY_MIN_DEPENDENCIES_PY_INFO = {
    "joblib": ("1.3.0", "install"),
    "scipy": ("1.10.0", "build, install"),
    "conda-lock": ("3.0.1", "maintenance"),
}

TOY_MATCHING_PYPROJECT_SECTIONS = """
[project]
dependencies = ["joblib>=1.3.0", "scipy>=1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0"]
install = ["joblib>=1.3.0", "scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0"]
"""

TOY_MATCHING_PYPROJECT_SECTIONS_WITH_UPPER_BOUND = """
[project]
dependencies = ["joblib>=1.3.0,<2.0", "scipy>=1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0,<1.19.0"]
install = ["joblib>=1.3.0,<2.0", "scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0,<1.19.0"]
"""

TOY_WRONG_SYMBOL_PYPROJECT_SECTIONS = """
[project]
dependencies = ["scipy<1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0"]
install = ["scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0"]
"""

TOY_MISSING_PACKAGE_PYPROJECT_SECTIONS = """
[project]
dependencies = ["scipy>=1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0"]
install = ["scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0"]
"""

TOY_ADDITIONAL_PACKAGE_PYPROJECT_SECTIONS = """
[project]
dependencies = ["joblib>=1.3.0", "scipy>=1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0", "package_not_in_min_dependencies_py_file>=4.2"]
install = ["joblib>=1.3.0", "scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0"]
"""

TOY_NON_MATCHING_VERSION_PYPROJECT_SECTIONS = """
[project]
dependencies = ["joblib>=1.42.0", "scipy>=1.10.0"]
[project.optional-dependencies]
build = ["scipy>=1.10.0"]
install = ["joblib>=1.3.0", "scipy>=1.10.0"]
maintenance = ["conda-lock==3.0.1"]
[build-system]
requires = ["scipy>=1.10.0"]
"""


def test_min_dependencies_readme():
    # Test that the minimum dependencies in the README.rst file are
    # consistent with the minimum dependencies defined at the file:
    # sklearn/_min_dependencies.py

    pattern = re.compile(
        r"\.\. \|"
        r"([A-Za-z-]+)"
        r"MinVersion\| replace::"
        r"( [0-9]+\.[0-9]+(\.[0-9]+)?)"
    )

    readme_path = Path(sklearn.__file__).parent.parent
    readme_file = readme_path / "README.rst"

    if not os.path.exists(readme_file):
        # Skip the test if the README.rst file is not available.
        # For instance, when installing scikit-learn from wheels
        pytest.skip("The README.rst file is not available.")

    with readme_file.open("r") as f:
        for line in f:
            matched = pattern.match(line)

            if not matched:
                continue

            package, version = matched.group(1), matched.group(2)
            package = package.lower()

            if package in dependent_packages:
                version = parse_version(version)
                min_version = parse_version(dependent_packages[package][0])

                message = (
                    f"{package} has inconsistent minimum versions in README.rst and"
                    f" _min_depencies.py: {version} != {min_version}"
                )
                assert version == min_version, message


def extract_packages_and_pyproject_tags(dependencies):
    min_depencies_tag_to_packages_without_version = defaultdict(list)
    for package, (min_version, tags) in dependencies.items():
        for t in tags.split(", "):
            min_depencies_tag_to_packages_without_version[t].append(package)

    pyproject_section_to_min_dependencies_tag = {
        "build-system.requires": "build",
        "project.dependencies": "install",
    }
    for tag in min_depencies_tag_to_packages_without_version:
        section = f"project.optional-dependencies.{tag}"
        pyproject_section_to_min_dependencies_tag[section] = tag

    return (
        min_depencies_tag_to_packages_without_version,
        pyproject_section_to_min_dependencies_tag,
    )


def check_pyproject_sections(pyproject_toml, min_dependencies):
    packages, pyproject_tags = extract_packages_and_pyproject_tags(min_dependencies)

    for pyproject_section, min_dependencies_tag in pyproject_tags.items():
        # Special situation for numpy: we have numpy>=2 in
        # build-system.requires to make sure we build wheels against numpy>=2.
        # TODO remove this when our minimum supported numpy version is >=2.
        skip_version_check_for = (
            ["numpy"] if pyproject_section == "build-system.requires" else []
        )

        expected_packages = packages[min_dependencies_tag]

        pyproject_section_keys = pyproject_section.split(".")
        info = pyproject_toml
        # iterate through nested keys to get packages and version
        for key in pyproject_section_keys:
            info = info[key]

        pyproject_build_min_versions = {}
        # Assuming pyproject.toml build section has something like "my-package>=2.3.0"
        pattern = r"([\w-]+)\s*[>=]=\s*([\d\w.]+)"
        for requirement in info:
            match = re.search(pattern, requirement)
            if match is None:
                raise NotImplementedError(
                    f"{requirement} does not match expected regex {pattern!r}. "
                    "Only >= and == are supported for version requirements"
                )

            package, version = match.group(1), match.group(2)

            pyproject_build_min_versions[package] = version

        msg = f"Packages in {pyproject_section} differ from _min_depencies.py"

        assert sorted(pyproject_build_min_versions) == sorted(expected_packages), msg

        for package, version in pyproject_build_min_versions.items():
            version = parse_version(version)
            expected_min_version = parse_version(min_dependencies[package][0])
            if package in skip_version_check_for:
                continue

            message = (
                f"{package} has inconsistent minimum versions in pyproject.toml and"
                f" _min_depencies.py: {version} != {expected_min_version}"
            )
            assert version == expected_min_version, message


def test_min_dependencies_pyproject_toml():
    """Check versions in pyproject.toml is consistent with _min_dependencies."""

    root_directory = Path(sklearn.__file__).parent.parent
    pyproject_toml_path = root_directory / "pyproject.toml"

    if not pyproject_toml_path.exists():
        # Skip the test if the pyproject.toml file is not available.
        # For instance, when installing scikit-learn from wheels
        pytest.skip("pyproject.toml is not available.")

    with pyproject_toml_path.open("rb") as f:
        pyproject_toml = tomllib.load(f)

    check_pyproject_sections(pyproject_toml, dependent_packages)


@pytest.mark.parametrize(
    "example_pyproject",
    [
        TOY_MATCHING_PYPROJECT_SECTIONS,
        TOY_MATCHING_PYPROJECT_SECTIONS_WITH_UPPER_BOUND,
    ],
)
def test_check_matching_pyproject_section(example_pyproject):
    """Test the version check for matching packages."""

    pyproject_toml = tomllib.loads(example_pyproject)

    check_pyproject_sections(pyproject_toml, TOY_MIN_DEPENDENCIES_PY_INFO)


@pytest.mark.parametrize(
    "example_non_matching_pyproject, error_msg",
    [
        (
            TOY_WRONG_SYMBOL_PYPROJECT_SECTIONS,
            ".* does not match expected regex .*. "
            "Only >= and == are supported for version requirements",
        ),
        (
            TOY_MISSING_PACKAGE_PYPROJECT_SECTIONS,
            "Packages in .* differ from _min_depencies.py",
        ),
        (
            TOY_ADDITIONAL_PACKAGE_PYPROJECT_SECTIONS,
            "Packages in .* differ from _min_depencies.py",
        ),
        (
            TOY_NON_MATCHING_VERSION_PYPROJECT_SECTIONS,
            ".* has inconsistent minimum versions in pyproject.toml and"
            " _min_depencies.py: .* != .*",
        ),
    ],
)
def test_check_non_matching_pyproject_section(
    example_non_matching_pyproject, error_msg
):
    """Test the version check for non-matching packages and versions."""

    pyproject_toml = tomllib.loads(example_non_matching_pyproject)

    with pytest.raises(Exception, match=error_msg):
        check_pyproject_sections(pyproject_toml, TOY_MIN_DEPENDENCIES_PY_INFO)
