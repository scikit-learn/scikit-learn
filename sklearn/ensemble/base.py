from ..base import BaseEstimator


class BaseEnsemble(BaseEstimator):
    """
    Base class for all ensemble classes
    """

    def __init__(self, estimator, n_estimators, **params):
        if not isinstance(estimator, BaseEstimator):
            raise TypeError("estimator must be a subclass of BaseEstimator")
        estimator.set_params(**params)
        self.estimators = [clone(estimator) for i in xrange(n_estimators)]

    def __len__(self):
        """Returns the number of estimators in the ensemble"""
        return len(self.estimators)

    def __getitem__(self, index):
        """Returns the index'th estimator in the ensemble"""
        return self.estimators[index]
