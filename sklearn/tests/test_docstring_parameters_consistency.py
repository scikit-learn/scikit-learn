# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn import metrics
from sklearn.ensemble import StackingClassifier, StackingRegressor
from sklearn.utils._testing import assert_docstring_consistency, skip_if_no_numpydoc

CLASS_DOCSTRING_CONSISTENCY_CASES = [
    {
        "objects": [StackingClassifier, StackingRegressor],
        "include_params": ["cv", "n_jobs", "passthrough", "verbose"],
        "exclude_params": None,
        "include_attrs": True,
        "exclude_attrs": ["final_estimator_"],
        "include_returns": False,
        "exclude_returns": None,
        "descr_regex_patterns": {},
    },
]

FUNCTION_DOCSTRING_CONSISTENCY_CASES = [
    {
        "objects": [
            metrics.precision_recall_fscore_support,
            metrics.f1_score,
            metrics.fbeta_score,
            metrics.precision_score,
            metrics.recall_score,
        ],
        "include_params": True,
        "exclude_params": [
            "zero_division",
        ],
        "include_attrs": False,
        "exclude_attrs": None,
        "include_returns": False,
        "exclude_returns": None,
        "descr_regex_patterns": {
            # Non-greedily matches anything until non-capturing group
            # "Weighted recall is equal to accuracy. " is found, then
            # matches anything until end of line. If not found first group
            # matches everything.
            "average": r"^(.+?)(?:Weighted recall is equal to accuracy\.\s(.*))?$",
        },
    },
]


@pytest.mark.parametrize("case", CLASS_DOCSTRING_CONSISTENCY_CASES)
@skip_if_no_numpydoc
def test_class_docstring_consistency(case):
    """Check docstrings parameters consistency between related classes."""
    assert_docstring_consistency(**case)


@pytest.mark.parametrize("case", FUNCTION_DOCSTRING_CONSISTENCY_CASES)
@skip_if_no_numpydoc
def test_function_docstring_consistency(case):
    """Check docstrings parameters consistency between related functions."""
    assert_docstring_consistency(**case)
