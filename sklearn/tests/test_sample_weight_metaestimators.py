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
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
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


EQUIVALENCE_ESTIMATORS = [
    pytest.param(_pipeline_logistic(), id="pipeline-logistic"),
    pytest.param(_grid_search_ridge(), id="gridsearch-pipeline-ridge"),
    pytest.param(_calibrated_logistic(), id="calibrated-logistic"),
]


@pytest.mark.parametrize("estimator", EQUIVALENCE_ESTIMATORS)
@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_integer_sample_weight_equivalence(estimator):
    """Integer weights and repeated rows agree on the meta-estimator's outputs.

    ``sample_weight`` is not passed to ``set_*_request``. Routing and
    auto-requests have to deliver it to every consuming sub-estimator for the
    fitted predictions to match.
    """
    check_sample_weight_equivalence_on_dense_data(
        estimator.__class__.__name__, estimator
    )
