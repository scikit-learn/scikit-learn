"""Tests for the minimum dependencies in the README.rst file and pyproject.toml."""


import os
import re
import platform
from pathlib import Path

import pytest
import sklearn
from sklearn._min_dependencies import dependent_packages, tag_to_packages
from sklearn.utils.fixes import parse_version


def test_min_dependencies_readme():
    # Test that the minimum dependencies in the README.rst file are
    # consistent with the minimum dependencies defined at the file:
    # sklearn/_min_dependencies.py

    if platform.python_implementation() == "PyPy":
        pytest.skip("PyPy does not always share the same minimum deps")

    pattern = re.compile(
        r"(\.\. \|)"
        + r"(([A-Za-z]+\-?)+)"
        + r"(MinVersion\| replace::)"
        + r"( [0-9]+\.[0-9]+(\.[0-9]+)?)"
    )

    readme_path = Path(sklearn.__path__[0]).parents[0]
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


def test_min_dependencies_pyproject_toml():
    """Verify that pyproject.toml is consistent with _min_dependencies."""
    root_path = Path(sklearn.__path__[0]).parents[0]
    pyproject_toml = root_path / "pyproject.toml"

    if not pyproject_toml.exists():
        # Skip the test if the pyproject.toml file is not available.
        # For instance, when installing scikit-learn from wheels
        pytest.skip("The pyproject.toml does not exist")

    toml_content = pyproject_toml.read_text()

    for tag, constraints in tag_to_packages.items():
        if tag == "maintenance":
            # maintenance does not need to be in pyproject.toml
            continue
        for constraint in constraints:
            assert (
                constraint in toml_content
            ), f"{constraint} should be in pyproject.toml"
