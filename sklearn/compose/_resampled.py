from ..base import BaseEstimator, MetaEstimatorMixin, clone
from ..utils.metaestimators import if_delegate_has_method


class ResampledTrainer(MetaEstimatorMixin, BaseEstimator):
    """Composition of a resampler and a predictor/transfomer

    Read more in the :ref:`User Guide <pipeline_resamplers>`.

    Parameters
    ----------
    resampler : Estimator supporting fit_resample
    predictor : Estimator

    Attributes
    ----------
    resampler_ : Estimator
        Fitted clone of `resampler`.

    predictor_ : Estimator
        Fitted clone of `predictor`.

    Examples
    --------
    >>> from sklearn.base import BaseEstimator
    >>> from sklearn.compose import ResampledTrainer
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.linear_model import LogisticRegression
    >>>
    >>> class HalfSampler(BaseEstimator):
    ...     "Train with every second sample"
    ...     def fit_resample(self, X, y, **kw):
    ...         return X[::2], y[::2]
    >>>
    >>> est = ResampledTrainer(HalfSampler(), LogisticRegression())
    >>> est.fit(X, y)
    >>> est.predict(X[:2])
    """

    def __init__(self, resampler, predictor):
        self.resampler = resampler
        self.predictor = predictor

    # TODO: tags?

    def fit(self, X, y=None, **kw):
        self.resampler_ = clone(self.resampler)
        ret = self.resampler_.fit_resample(X, y, **kw)
        if len(ret) == 2:
            kw = {}
            X, y = ret
        else:
            X, y, kw = ret
        self.predictor_ = clone(self.predictor).fit(X, y, **kw)
        return self

    @if_delegate_has_method(delegate='predictor_')
    def predict(self, X, **predict_params):
        return self.predictor_.predict(X, **predict_params)

    @if_delegate_has_method(delegate='predictor_')
    def predict_proba(self, X):
        return self.predictor_.predict_proba(X)

    @if_delegate_has_method(delegate='predictor_')
    def predict_log_proba(self, X):
        return self.predictor_.predict_log_proba(X)

    @if_delegate_has_method(delegate='predictor_')
    def decision_function(self, X):
        return self.predictor_.decision_function(X)

    @if_delegate_has_method(delegate='predictor_')
    def score(self, X, y, **kw):
        return self.predictor_.score(X, y, **kw)

    @property
    def fit_transform(self):
        transform = self.predictor_.transform

        def fit_transform(X, y, **kwargs):
            self.fit(X, y, **kwargs)
            return transform(X)

    @property
    def fit_predict(self):
        predict = self.predictor_.predict

        def fit_predict(X, y, **kwargs):
            self.fit(X, y, **kwargs)
            return predict(X)

    @property
    def _estimator_type(self):
        return self.predictor._estimator_type

    @property
    def classes_(self):
        return self.predictor_.classes_
