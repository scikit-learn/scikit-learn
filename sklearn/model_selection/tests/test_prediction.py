from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

from sklearn.model_selection import CutOffClassifier


def test_xxx():
    X, y = make_classification(random_state=0)
    clf = CutOffClassifier(LogisticRegression())
    clf.fit(X, y)

    assert clf.predict(X) is not None
