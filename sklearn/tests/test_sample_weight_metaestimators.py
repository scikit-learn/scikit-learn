"""Tests that composed meta-estimators honor integer ``sample_weight``.

The equivalence test compares top-level predictions of real meta-estimators
fitted with integer weights against the same estimators fitted on repeated
rows.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn import config_context
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import make_scorer, mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler
from sklearn.utils.estimator_checks import check_sample_weight_equivalence_on_dense_data


def _pipeline_logistic():
    return Pipeline(
        [
            ("scale", StandardScaler()),
            ("model", LogisticRegression(solver="lbfgs")),
        ]
    )


def _grid_search_ridge():
    return GridSearchCV(
        Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", Ridge()),
            ]
        ),
        param_grid={"model__alpha": [0.1, 1.0, 10.0]},
        scoring=make_scorer(mean_squared_error, greater_is_better=False),
    )


def _calibrated_logistic():
    return CalibratedClassifierCV(
        LogisticRegression(solver="lbfgs"),
        method="sigmoid",
    )


def _grid_search_spline_union_logistic():
    # Same composition as the meta-estimator in
    # https://github.com/scikit-learn/scikit-learn/pull/30887
    return GridSearchCV(
        Pipeline(
            [
                (
                    "union",
                    FeatureUnion(
                        [
                            ("scaler", StandardScaler()),
                            ("spline", SplineTransformer(knots="quantile")),
                        ]
                    ),
                ),
                ("logistic", LogisticRegression(solver="lbfgs")),
            ]
        ),
        param_grid={
            "union__spline__n_knots": [5, 6],
            "logistic__C": [0.1, 1.0],
        },
    )


@pytest.mark.parametrize(
    "estimator",
    [
        _pipeline_logistic(),
        _grid_search_ridge(),
        _calibrated_logistic(),
        _grid_search_spline_union_logistic(),
    ],
    ids=[
        "pipeline-logistic",
        "gridsearch-pipeline-ridge",
        "calibrated-logistic",
        "gridsearch-spline-union-logistic",
    ],
)
@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_integer_sample_weight_equivalence(estimator):
    """Integer weights and repeated rows agree on the meta-estimator's outputs.

    Auto-requests are enabled, so the test does not call ``set_*_request``. The
    outputs should match only when routing delivers ``sample_weight`` to every
    consuming sub-estimator and scorer, and each of those pieces itself treats
    integer weights as repeated rows.
    """
    check_sample_weight_equivalence_on_dense_data(
        estimator.__class__.__name__, estimator
    )
