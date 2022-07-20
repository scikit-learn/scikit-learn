"""Test Multioutput."""


def test_multioutput():
    from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
    from sklearn.datasets import make_regression, make_multilabel_classification
    from sklearn.multioutput import RegressorChain, ClassifierChain, MultiOutputRegressor, MultiOutputClassifier
    import numpy as np

    X, y = make_regression(n_targets=2)
    X[0, 0] = np.nan

    cl = RegressorChain(HistGradientBoostingRegressor())
    assert cl._get_tags()['allow_nan']
    cl.fit(X, y)
    assert cl.predict(X).shape[1] == 2

    cl = MultiOutputRegressor(HistGradientBoostingRegressor())
    assert cl._get_tags()['allow_nan']
    cl.fit(X, y)
    assert cl.predict(X).shape[1] == 2

    X, y = make_multilabel_classification(n_classes=2)
    X[0, 0] = np.nan

    cl = ClassifierChain(HistGradientBoostingClassifier())
    assert cl._get_tags()['allow_nan']
    cl.fit(X, y)
    assert cl.predict(X).shape[1] == 2

    cl = MultiOutputClassifier(HistGradientBoostingClassifier())
    assert cl._get_tags()['allow_nan']
    cl.fit(X, y)
    assert cl.predict(X).shape[1] == 2

test_multioutput()