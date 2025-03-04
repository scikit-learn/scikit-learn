# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn import metrics
from sklearn.ensemble import (
    AdaBoostClassifier,
    AdaBoostRegressor,
    StackingClassifier,
    StackingRegressor,
)
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture
from sklearn.utils._testing import assert_docstring_consistency, skip_if_no_numpydoc

CLASS_DOCSTRING_CONSISTENCY_CASES = [
    {
        "objects": [AdaBoostClassifier, AdaBoostRegressor],
        "include_params": True,
        "exclude_params": None,
        "include_attrs": True,
        "exclude_attrs": [
            # type differs
            "estimators_",
            "",
        ],
        "include_returns": False,
        "exclude_returns": None,
        "descr_regex_patterns": {
            # Excludes 2nd sentence if present, matches last sentence until "."
            "estimator": r"^([^.]+\.)(?:\s+[^.]+\.)?\s+([^.]+\.)",
            "learning_rate": (
                r"^(.*?)(?:regressor |classifier )(.*?)(?:regressor|classifier)(.*)$"
            ),
            # Excludes 3rd sentence if present.
            "random_state": r"([^.]+\.[^.]+\.)(?:\s+[^.]+\.)?(\s[^.]+\.[^.]+\.)",
            # Excludes words "Classification "/"Regression "
            "estimator_errors_": r"^(?:Classification |Regression )(.*)$",
        },
    },
    {
        "objects": [GaussianMixture, BayesianGaussianMixture],
        "include_params": True,
        "exclude_params": None,
        "include_attrs": True,
        "exclude_attrs": ["lower_bound_"],
        "include_returns": False,
        "exclude_returns": None,
        "descr_regex_patterns": {
            # Match first sentence only
            "n_components": r"^[^.?!]*[.]",
            # Excludes words "precisions" or "covariances"
            "init_params": r"^(.*?)(?:precisions.|covariances.)(.*)$",
            # Excludes words "inference " or "EM "
            "converged_": r"^(.*?)(?:inference |EM )(.*)$",
            # Excludes words "inference " or "EM "
            "n_iter_": r"^(.*?)(?:inference |EM )(.*)$",
        },
    },
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
            # Matches everything, excluding "Weighted recall is equal to accuracy. "
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
