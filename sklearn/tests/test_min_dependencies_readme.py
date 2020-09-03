"""Tests for the minimum dependencies in the README.rst file."""

import os
import re
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

    readme_path = os.path.abspath(os.path.join(sklearn.__path__[0], ".."))
    print(os.listdir(readme_path))
    readme_filename = os.path.join(readme_path, "README.rst")

    # if not os.path.exists(readme_filename):
    # Skip because the README.rst
    # file is not found in Windows
    # pytest.skip("README.rst is not available.")

    with open(readme_filename, "r") as f:
        for line in f:
            if pattern.match(line):
                dependency = pattern.sub(r"\2\5", line)
                (package, version) = dependency.lower().split(" ")

                if package in dependent_packages:
                    version = parse(version)
                    min_version = parse(dependent_packages[package][0])

                    assert version == min_version
