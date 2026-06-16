#!/bin/bash

set -e

# Verify that scikit-learn's estimator checks can be imported and used when
# pytest is NOT installed (pytest is only a soft dependency for end users).
#
# This script runs in the dedicated `pytest-soft-dependency` pixi environment
# (see pyproject.toml), which deliberately does not include pytest. We therefore
# do not modify the environment here: there is nothing to uninstall.
#
# `sklearn/utils/tests/test_estimator_checks.py` is intentionally written
# without pytest (it relies on `unittest` and a local `raises` helper) so that
# it can be executed with a plain `python -m`.

mkdir -p "$TEST_DIR"
# Run from outside the source tree so that the installed scikit-learn is used.
cd "$TEST_DIR"
python -m sklearn.utils.tests.test_estimator_checks
