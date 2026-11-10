# This script checks that the common tests marked with xfail are actually
# failing.
# Note that in some cases, a test might be marked with xfail because it is
# failing on certain machines, and might not be triggered by this script.

import contextlib
import io

from sklearn.utils._test_common.instance_generator import (
    _get_expected_failed_checks,
    _tested_estimators,
)
from sklearn.utils.estimator_checks import check_estimator

for estimator in _tested_estimators():
    # calling check_estimator w/o passing expected_failed_checks will find
    # all the failing tests in your environment.
    # suppress stdout/stderr while running checks
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        check_results = check_estimator(estimator, on_skip=None, on_fail=None)
    failed_tests = [e for e in check_results if e["status"] == "failed"]
    failed_test_names = set(e["check_name"] for e in failed_tests)
    expected_failed_tests = set(_get_expected_failed_checks(estimator).keys())
    unexpected_failures = failed_test_names - expected_failed_tests
    if unexpected_failures:
        print(f"{estimator.__class__.__name__} failed with unexpected failures:")
        for failure in unexpected_failures:
            print(f"  {failure}")

    expected_but_not_raised = expected_failed_tests - failed_test_names
    if expected_but_not_raised:
        print(f"{estimator.__class__.__name__} did not fail expected failures:")
        for failure in expected_but_not_raised:
            print(f"  {failure}")
