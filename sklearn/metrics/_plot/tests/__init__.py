import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.datasets import load_iris
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split


@pytest.fixture(scope="module")
def iris():
    return load_iris()


@pytest.fixture(scope="module")
def clf_iris(iris):
    X_train, _, y_train, _ = train_test_split(iris.data, iris.target,
                                              random_state=42)
    clf = SVC(kernel='linear', C=0.01)
    return clf.fit(X_train, y_train)


@pytest.mark.parametrize("normalize", [True, False])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_plot_confusion_matrix(pyplot, iris, clf_iris, with_sample_weight,
                               normalize):
    _, X_test, _, y_test = train_test_split(iris.data, iris.target,
                                            random_state=42)
    if with_sample_weight:
        rng = np.random.RandomState(42)
        sample_weight = rng.randint(1, 4, size=(X_test.shape[0]))
    else:
        sample_weight = None

    y_pred = clf_iris.predict(X_test)
    cm = confusion_matrix(y_test, y_pred, labels=clf_iris.classes_,
                          sample_weight=sample_weight)
    if normalize:
        cm = cm.astype(np.float) / cm.sum(axies=1)[:, None]

    disp = plot_confusion_matrix(iris, y_test, sample_weight=sample_weight,
                                  normalize=normalize)

    assert_allclose(disp.confusion_matrix, cm)
    assert np.all(disp.target_names == clf_iris.classes_)
