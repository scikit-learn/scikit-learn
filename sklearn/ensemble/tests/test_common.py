import pytest

from sklearn.base import clone

from sklearn.datasets import make_classification
from sklearn.datasets import make_regression

from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.svm import LinearSVC, LinearSVR
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from sklearn.ensemble import StackingClassifier, StackingRegressor
from sklearn.ensemble import VotingClassifier, VotingRegressor


@pytest.mark.parametrize(
    "X, y, estimator",
    [(*make_classification(n_samples=10),
      StackingClassifier(estimators=[('lr', LogisticRegression()),
                                     ('svm', LinearSVC()),
                                     ('rf', RandomForestClassifier())])),
     (*make_classification(n_samples=10),
      VotingClassifier(estimators=[('lr', LogisticRegression()),
                                   ('svm', LinearSVC()),
                                   ('rf', RandomForestClassifier())])),
     (*make_regression(n_samples=10),
      StackingRegressor(estimators=[('lr', LinearRegression()),
                                    ('svm', LinearSVR()),
                                    ('rf', RandomForestRegressor())])),
     (*make_regression(n_samples=10),
      VotingRegressor(estimators=[('lr', LinearRegression()),
                                  ('svm', LinearSVR()),
                                  ('rf', RandomForestRegressor())]))]
)
def test_ensemble_heterogeneous(X, y, estimator):
    estimator.fit(X, y)
    assert len(estimator.named_estimators) == 3
    assert len(estimator.named_estimators_) == 3

    estimator_dropped = clone(estimator)
    estimator_dropped.set_params(svm='drop')
    estimator_dropped.fit(X, y)
    # be sure that we are not reporting the dropped estimator once fitted
    assert len(estimator_dropped.named_estimators) == 3
    assert estimator_dropped.named_estimators.svm == 'drop'
    assert len(estimator_dropped.named_estimators_) == 2
    assert (sorted(list(estimator_dropped.named_estimators_.keys())) ==
            sorted(['lr', 'rf']))
    for sub_est in estimator_dropped.named_estimators_:
        # check that the correspondence is correct
        assert not isinstance(sub_est, type(estimator.named_estimators.svm))
