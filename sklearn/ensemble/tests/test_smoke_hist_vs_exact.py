# sklearn/ensemble/tests/test_smoke_hist_vs_exact.py
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier

def test_smoke_random_forest_small():
    X, y = make_classification(n_samples=200, n_features=10, random_state=0)
    clf = RandomForestClassifier(n_estimators=5, random_state=0)
    clf.fit(X, y)
    preds = clf.predict(X[:10])
    assert len(preds) == 10
