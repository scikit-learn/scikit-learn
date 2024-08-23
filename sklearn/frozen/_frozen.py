# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..base import BaseEstimator
from ..exceptions import NotFittedError
from ..utils._available_if import available_if
from ..utils.validation import check_is_fitted


def _estimator_has(attr):
    """Check that final_estimator has `attr`.

    Used together with `available_if`.
    """

    def check(self):
        # raise original `AttributeError` if `attr` does not exist
        getattr(self.estimator, attr)
        return True

    return check


def _estimator_fitted():
    """Check that final_estimator has `attr`.

    Used together with `available_if`.
    """

    def check(self):
        # raise original `AttributeError` if `attr` does not exist
        try:
            check_is_fitted(self)
            return True
        except NotFittedError:
            return False

    return check


class FrozenEstimator(BaseEstimator):
    """Frozen estimator.

    This meta-estimator takes an estimator and freezes it, in the sense that calling
    `.fit` on it has no effect.

    This is particularly useful when you have a fit or a pre-trained model as a
    transformer in a pipeline, and you'd like `pipeline.fit` to have no effect on this
    step.

    Parameters
    ----------
    estimator : estimator
        The estimator which is to be kept frozen.

    See Also
    --------
    None: No siimlar entry in the scikit-learn documentation.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.frozen import FrozenEstimator
    >>> from sklearn.linear_model import LogisticRegression
    >>> X, y = make_classification(random_state=0)
    >>> clf = LogisticRegression(random_state=0).fit(X, y)
    >>> frozen_clf = FrozenEstimator(clf)
    >>> frozen_clf.fit(X, y)  # No-op
    >>> frozen_clf.predict(X)  # Predictions from `clf.predict`
    """

    _required_parameters = ["estimator"]

    def __init__(self, estimator):
        self.estimator = estimator

    @available_if(_estimator_has("__getitem__"))
    def __getitem__(self, *args, **kwargs):
        return self.estimator.__getitem__(*args, **kwargs)

    def __getattr__(self, name):
        # `estimator`'s attributes are now accessible
        return getattr(self.estimator, name)

    def __sklearn_clone__(self):
        return self

    def fit(self, X, y, *args, **kwargs):
        """No-op.

        As a frozen estimator, calling fit has no effect.

        Parameters
        ----------
        X : object
            Ignored.

        y : object
            Ignored.

        *args : tuple
            Additional positional arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        **kwargs : dict
            Additional keyword arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        check_is_fitted(self.estimator)
        return self

    @available_if(_estimator_has("fit_transform") and _estimator_fitted())
    def fit_transform(self, X, y=None, *args, **kwargs):
        """Call `estimator.transform`.

        As a frozen estimator, the call to `fit` is skipped.

        Parameters
        ----------
        X : object
            Input data.

        y : object
            Ignored.

        *args : tuple
            Additional positional arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        **kwargs : dict
            Additional keyword arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        Returns
        -------
        X_transformed : object
            Transformed data.
        """
        # fit_transform only transforms the data
        return self.estimator.transform(X, *args, **kwargs)

    @available_if(_estimator_has("fit_predict") and _estimator_fitted())
    def fit_predict(self, X, y=None, *args, **kwargs):
        """Call `estimator.predict`.

        As a frozen estimator, the call to `fit` is skipped.

        Parameters
        ----------
        X : object
            Input data.

        y : object
            Ignored.

        *args : tuple
            Additional positional arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        **kwargs : dict
            Additional keyword arguments. Ignored, but present for API compatibility
            with `self.estimator`.

        Returns
        -------
        y_pred : object
            Predictions from `self.estimator.predict`.
        """
        # fit_transform only transforms the data
        return self.estimator.predict(X, *args, **kwargs)

    def _more_tags(self):
        return {
            "_skip_test": True,
        }
