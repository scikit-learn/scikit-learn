"""Tests for the minimum dependencies in the README.rst file."""

import os
import re
from packaging.version import parse

import sklearn
from sklearn._build_utils.min_dependencies import dependent_packages


def test_min_dependencies_readme():
    # Test that the minimum dependencies in the README.rst file
    # are compliant with the minimum dependencies for scikit-learn
    pattern = re.compile(r"(\.\. \|)" +
                         r"(([A-Za-z]+\-?)+)" +
                         r"(MinVersion\| replace::)" +
                         r"( [0-9]+\.[0-9]+(\.[0-9]+)?)")

    with open(os.path.join(sklearn.__path__[0], "..", "README.rst"), "r") as f:
        for line in f:
            if pattern.match(line):
                dependency = pattern.sub(r"\2\5", line)
                (package, version) = dependency.lower().split(" ")

                # Avoid non-defined packages
                if package in dependent_packages:
                    version = parse(version)
                    min_version = parse(dependent_packages[package][0])

                    assert version == min_version
