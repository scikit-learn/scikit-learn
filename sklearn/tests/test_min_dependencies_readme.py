"""Tests for the minimum dependencies in the README.rst file."""


import os
import re
from pathlib import Path
from packaging.version import parse

import pytest
import sklearn
from sklearn._build_utils.min_dependencies import dependent_packages


def test_min_dependencies_readme():
    # Test that the minimum dependencies in the README.rst file are
    # consistent with the minimum dependencies defined at the file:
    # sklearn/_build_utils/min_dependencies.py

    pattern = re.compile(r"(\.\. \|)" +
                         r"(([A-Za-z]+\-?)+)" +
                         r"(MinVersion\| replace::)" +
                         r"( [0-9]+\.[0-9]+(\.[0-9]+)?)")

    readme_path = Path(sklearn.__path__[0]).parents[0]
    readme_file = readme_path / "README.rst"

    if not os.path.exists(readme_file):
        # Skip the test if the README.rst file is not available.
        # For instance, when installing scikit-learn from wheels
        pytest.skip("The README.rst file is not available.")

    with readme_file.open("r") as f:
        for line in f:
            if pattern.match(line):
                dependency = pattern.sub(r"\2\5", line)
                (package, version) = dependency.lower().split(" ")

                if package in dependent_packages:
                    version = parse(version)
                    min_version = parse(dependent_packages[package][0])

                    assert version == min_version
