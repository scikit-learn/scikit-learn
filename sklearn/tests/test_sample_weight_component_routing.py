"""Check that auto-requested ``sample_weight`` reaches each sub-component.

The meta-estimator here is only a router. Its children mimic the pieces used
in composed estimators: a transformer, a meta-estimator that both consumes and
routes ``sample_weight``, a regressor, and a scorer. No ``set_*_request`` call
is used, so a method receives ``sample_weight`` only when its signature accepts
it.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from functools import partial

import numpy as np

from sklearn import config_context
from sklearn.base import (
    BaseEstimator,
    MetaEstimatorMixin,
    RegressorMixin,
    TransformerMixin,
    clone,
)
from sklearn.metrics import mean_squared_error
from sklearn.metrics._scorer import _Scorer
from sklearn.tests.metadata_routing_common import (
    check_recorded_metadata,
    record_metadata,
)
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    get_routing_for_object,
    process_routing,
)


class RecordingTransformer(TransformerMixin, BaseEstimator):
    """Transformer that consumes ``sample_weight`` in ``fit`` only.

    Mirrors ``StandardScaler``: ``transform`` does not accept ``sample_weight``.
    """

    def fit(self, X, y=None, sample_weight=None):
        record_metadata(self, sample_weight=sample_weight)
        self.fitted_ = True
        return self

    def transform(self, X):
        return X


class RecordingRegressor(RegressorMixin, BaseEstimator):
    """Regressor that consumes ``sample_weight`` in ``fit`` only.

    Mirrors ``Ridge``: ``predict`` does not accept ``sample_weight``.
    """

    def fit(self, X, y, sample_weight=None):
        record_metadata(self, sample_weight=sample_weight)
        self.fitted_ = True
        return self

    def predict(self, X):
        return np.zeros(len(X))


class RecordingConsumerRouter(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    """Meta-estimator that consumes ``sample_weight`` and routes it onward.

    Mirrors meta-estimators such as ``CalibratedClassifierCV``: ``fit`` uses
    ``sample_weight`` itself and forwards it to ``estimator.fit``.
    """

    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y, sample_weight=None, **fit_params):
        record_metadata(self, sample_weight=sample_weight)
        params = process_routing(self, "fit", sample_weight=sample_weight, **fit_params)
        self.estimator_ = clone(self.estimator).fit(X, y, **params.estimator.fit)
        return self

    def predict(self, X):
        return self.estimator_.predict(X)

    def get_metadata_routing(self):
        return (
            MetadataRouter(owner=self)
            .add_self_request(self)
            .add(
                estimator=self.estimator,
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
        )


def consuming_metric(y_true, y_pred, sample_weight=None, scorer=None):
    """Metric with an explicit ``sample_weight`` parameter.

    The function name is what ``record_metadata`` uses as the callee.
    """
    record_metadata(scorer, sample_weight=sample_weight)
    return mean_squared_error(y_true, y_pred, sample_weight=sample_weight)


class RecordingScorer(_Scorer):
    """Scorer that auto-requests ``sample_weight``, like ``make_scorer``."""

    def __init__(self):
        super().__init__(
            score_func=partial(consuming_metric, scorer=self),
            sign=1,
            kwargs={},
            response_method="predict",
        )


class ComposedMeta(MetaEstimatorMixin, RegressorMixin, BaseEstimator):
    """Router around a transformer, a consuming meta-estimator, and a scorer.

    ``fit`` does not consume ``sample_weight`` itself. It only forwards the
    routed values to each child.
    """

    def __init__(self, transformer, estimator, scorer):
        self.transformer = transformer
        self.estimator = estimator
        self.scorer = scorer

    def fit(self, X, y, sample_weight=None):
        routed = process_routing(self, "fit", sample_weight=sample_weight)
        self.transformer_ = clone(self.transformer).fit(X, y, **routed.transformer.fit)
        Xt = self.transformer_.transform(X)
        self.estimator_ = clone(self.estimator).fit(Xt, y, **routed.estimator.fit)
        self.score_ = self.scorer(self.estimator_, Xt, y, **routed.scorer.score)
        return self

    def get_metadata_routing(self):
        return (
            MetadataRouter(owner=self)
            .add(
                transformer=self.transformer,
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
            .add(
                estimator=self.estimator,
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
            .add(
                scorer=self.scorer,
                method_mapping=MethodMapping().add(caller="fit", callee="score"),
            )
        )


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_sample_weight_routed_to_each_component():
    """The same ``sample_weight`` array reaches every consuming component.

    Auto-requests are enabled, so the test does not call ``set_*_request``.
    ``transform`` and ``predict`` do not accept ``sample_weight`` and must not
    gain a request for it.
    """
    transformer = RecordingTransformer()
    regressor = RecordingRegressor()
    router = RecordingConsumerRouter(estimator=regressor)
    scorer = RecordingScorer()
    sample_weight = np.array([1, 0, 2, 1, 3, 0, 4, 2])
    X = np.arange(sample_weight.shape[0] * 3, dtype=float).reshape(-1, 3)
    y = np.arange(sample_weight.shape[0], dtype=float)

    transformer_requests = get_routing_for_object(transformer)
    assert transformer_requests.fit.requests["sample_weight"] is True
    assert "sample_weight" not in transformer_requests.transform.requests

    regressor_requests = get_routing_for_object(regressor)
    assert regressor_requests.fit.requests["sample_weight"] is True
    assert "sample_weight" not in regressor_requests.predict.requests

    scorer_requests = get_routing_for_object(scorer)
    assert scorer_requests.score.requests["sample_weight"] is True

    meta = ComposedMeta(transformer=transformer, estimator=router, scorer=scorer)
    meta.fit(X, y, sample_weight=sample_weight)

    check_recorded_metadata(
        meta.transformer_,
        method="fit",
        parent="fit",
        sample_weight=sample_weight,
    )
    check_recorded_metadata(
        meta.estimator_,
        method="fit",
        parent="fit",
        sample_weight=sample_weight,
    )
    check_recorded_metadata(
        meta.estimator_.estimator_,
        method="fit",
        parent="fit",
        sample_weight=sample_weight,
    )
    check_recorded_metadata(
        meta.scorer,
        method="consuming_metric",
        parent="fit",
        sample_weight=sample_weight,
    )
