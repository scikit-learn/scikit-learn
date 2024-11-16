"""Tests for the minimum dependencies in README.rst and pyproject.toml"""

import os
import re
from collections import defaultdict
from pathlib import Path

import pytest

import sklearn
from sklearn._min_dependencies import dependent_packages
from sklearn.utils.fixes import parse_version

min_depencies_tag_to_packages_without_version = defaultdict(list)
for package, (min_version, extras) in dependent_packages.items():
    for extra in extras.split(", "):
        min_depencies_tag_to_packages_without_version[extra].append(package)

min_dependencies_tag_to_pyproject_section = {
    "build": "build-system.requires",
    "install": "project.dependencies",
}
for tag in min_depencies_tag_to_packages_without_version:
    min_dependencies_tag_to_pyproject_section[tag] = (
        f"project.optional-dependencies.{tag}"
    )


def test_min_dependencies_readme():
    # Test that the minimum dependencies in the README.rst file are
    # consistent with the minimum dependencies defined at the file:
    # sklearn/_min_dependencies.py

    pattern = re.compile(
        r"(\.\. \|)"
        + r"(([A-Za-z]+\-?)+)"
        + r"(MinVersion\| replace::)"
        + r"( [0-9]+\.[0-9]+(\.[0-9]+)?)"
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

            package, version = matched.group(2), matched.group(5)
            package = package.lower()

            if package in dependent_packages:
                version = parse_version(version)
                min_version = parse_version(dependent_packages[package][0])

                assert version == min_version, f"{package} has a mismatched version"


def check_pyproject_section(
    pyproject_section, min_dependencies_tag, skip_version_check_for=None
):
    # tomllib is available in Python 3.11
    tomllib = pytest.importorskip("tomllib")

    if skip_version_check_for is None:
        skip_version_check_for = []

    expected_packages = min_depencies_tag_to_packages_without_version[
        min_dependencies_tag
    ]

    root_directory = Path(sklearn.__file__).parent.parent
    pyproject_toml_path = root_directory / "pyproject.toml"

    if not pyproject_toml_path.exists():
        # Skip the test if the pyproject.toml file is not available.
        # For instance, when installing scikit-learn from wheels
        pytest.skip("pyproject.toml is not available.")

    with pyproject_toml_path.open("rb") as f:
        pyproject_toml = tomllib.load(f)

    pyproject_section_keys = pyproject_section.split(".")
    info = pyproject_toml
    for key in pyproject_section_keys:
        info = info[key]

    pyproject_build_min_versions = {}
    for requirement in info:
        if ">=" in requirement:
            package, version = requirement.split(">=")
        elif "==" in requirement:
            package, version = requirement.split("==")
        else:
            raise NotImplementedError(
                f"{requirement} not supported yet in this test. "
                "Only >= and == are supported for version requirements"
            )

        pyproject_build_min_versions[package] = version

    assert sorted(pyproject_build_min_versions) == sorted(expected_packages)

    for package, version in pyproject_build_min_versions.items():
        version = parse_version(version)
        expected_min_version = parse_version(dependent_packages[package][0])
        if package in skip_version_check_for:
            continue

        assert version == expected_min_version, f"{package} has a mismatched version"


@pytest.mark.parametrize(
    "min_dependencies_tag, pyproject_section",
    min_dependencies_tag_to_pyproject_section.items(),
)
def test_min_dependencies_pyproject_toml(pyproject_section, min_dependencies_tag):
    """Check versions in pyproject.toml is consistent with _min_dependencies."""
    # NumPy is more complex because build-time (>=1.25) and run-time (>=1.19.5)
    # requirement currently don't match
    skip_version_check_for = ["numpy"] if min_dependencies_tag == "build" else None
    check_pyproject_section(
        pyproject_section,
        min_dependencies_tag,
        skip_version_check_for=skip_version_check_for,
    )
