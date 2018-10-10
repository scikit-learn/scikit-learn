import numpy as np
import pytest
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns_message
from sklearn.exceptions import NotFittedError

from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn.dummy import DummyClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_iris
from sklearn.preprocessing import LabelBinarizer

# Author: Oliver Rausch <rauscho@ethz.ch>
# License: BSD 3 clause


# load the iris dataset and randomly permute it
iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                    iris.target,
                                                    random_state=0)

n_missing_samples = 50

y_train_missing_labels = y_train.copy()
y_train_missing_labels[n_missing_samples:] = -1
y_train_missing_labels_strings = y_train_missing_labels.copy()
mapping = {0: 'A', 1: 'B', 2: 'C', -1: -1}
y_train_missing_strings = np.vectorize(mapping.get)(y_train_missing_labels)

lb = LabelBinarizer()
y_train_missing_dummies = lb.fit_transform(y_train_missing_labels)


def test_missing_predict_proba():
    # Check that an error is thrown if predict_proba is not implemented
    base_classifier = SVC(gamma="scale")
    self_training = SelfTrainingClassifier(base_classifier)
    message = "base_classifier (SVC) should implement predict_proba!"
    assert_raise_message(ValueError, message, self_training.fit, X_train,
                         y_train)


def test_none_classifier():
    st = SelfTrainingClassifier(None)
    msg = "base_classifier cannot be None"
    assert_raise_message(ValueError, msg, st.fit, X_train,
                         y_train_missing_labels)


@pytest.mark.parametrize("max_iter, threshold",
                         [(-1, 1.0), (-100, -2), (-10, 10)])
def test_invalid_params(max_iter, threshold):
    # Test negative iterations
    base_classifier = SVC(gamma="scale", probability=True)
    st = SelfTrainingClassifier(base_classifier, max_iter=max_iter)
    message = "max_iter must be >= 0 or None, got"
    assert_raise_message(ValueError, message, st.fit, X_train, y_train)

    base_classifier = SVC(gamma="scale", probability=True)
    st = SelfTrainingClassifier(base_classifier, threshold=threshold)
    message = "threshold must be in [0,1)"
    assert_raise_message(ValueError, message, st.fit, X_train, y_train)


@pytest.mark.parametrize("base_classifier",
                         [DummyClassifier(random_state=0),
                          DecisionTreeClassifier(random_state=0),
                          KNeighborsClassifier(),
                          SVC(gamma="scale", probability=True,
                              random_state=0)])
@pytest.mark.parametrize("max_iter", [1, 50, 100])
@pytest.mark.parametrize("threshold", [0.0, 0.5, 0.9])
def test_classification(base_classifier, max_iter, threshold):
    # Check classification for various parameter settings.
    # Also assert that predictions for strings and numerical labels are equal.
    # Also test for multioutput classification
    st = SelfTrainingClassifier(base_classifier, max_iter=max_iter,
                                threshold=threshold)
    st.fit(X_train, y_train_missing_labels)
    pred = st.predict(X_test)
    proba = st.predict_proba(X_test)

    st_string = SelfTrainingClassifier(base_classifier, max_iter=max_iter,
                                       threshold=threshold)
    st_string.fit(X_train, y_train_missing_strings)
    pred_string = st_string.predict(X_test)
    proba_string = st_string.predict_proba(X_test)

    assert_array_equal(np.vectorize(mapping.get)(pred), pred_string)
    assert_array_equal(proba, proba_string)


def test_none_iter():
    # Check that the all samples were labeled after a 'reasonable' number of
    # iterations.
    st = SelfTrainingClassifier(KNeighborsClassifier(), threshold=.55,
                                max_iter=None)
    st.fit(X_train, y_train_missing_labels)

    assert st.n_iter_ < 10


def test_zero_iterations():
    # Check classification for zero iterations.
    # Fitting a SelfTrainingClassifier with zero iterations should give the
    # same results as fitting a supervised classifier.

    clf1 = SelfTrainingClassifier(KNeighborsClassifier(),
                                  max_iter=0).fit(X_train,
                                                  y_train_missing_labels)

    clf2 = KNeighborsClassifier().fit(
        X_train[:n_missing_samples], y_train[:n_missing_samples])

    assert_array_equal(clf1.predict(X_test), clf2.predict(X_test))


def test_notfitted():
    # Test that predicting without training throws an error
    st = SelfTrainingClassifier(KNeighborsClassifier())
    msg = ("This SelfTrainingClassifier instance is not fitted yet. Call "
           "\'fit\' with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, st.predict, X_train)
    assert_raise_message(NotFittedError, msg, st.predict_proba, X_train)


def test_prefitted_throws_error():
    # Test that passing a pre-fitted classifier and calling predict throws an
    # error
    knn = KNeighborsClassifier()
    knn.fit(X_train, y_train)
    st = SelfTrainingClassifier(knn)
    msg = ("This SelfTrainingClassifier instance is not fitted yet. Call "
           "\'fit\' with appropriate arguments before using this method.")
    assert_raise_message(NotFittedError, msg, st.predict, X_train)


@pytest.mark.parametrize("max_iter", range(1, 5))
def test_y_labeled_iter(max_iter):
    # Check that the amount of datapoints labeled in iteration 0 is equal to
    # the amount of labeled datapoints we passed.
    st = SelfTrainingClassifier(KNeighborsClassifier(), max_iter=max_iter)
    st.fit(X_train, y_train_missing_labels)
    amount_iter_0 = len(st.y_labeled_iter_[st.y_labeled_iter_ == 0])
    assert amount_iter_0 == n_missing_samples
    # Check that the max of the iterations is less than the total amount of
    # iterations
    assert np.max(st.y_labeled_iter_) <= st.n_iter_ <= max_iter


def test_no_unlabeled():
    # Test that training on a fully labeled dataset produces the same results
    # as training the classifier by itself.
    knn = KNeighborsClassifier()
    knn.fit(X_train, y_train)
    st = SelfTrainingClassifier(knn)
    msg = "y contains no unlabeled samples"
    assert_warns_message(RuntimeWarning, msg, st.fit, X_train, y_train)
    assert_array_equal(knn.predict(X_test), st.predict(X_test))
    # Assert that all samples were labeled in iteration 0 (since there were no
    # unlabeled samples).
    assert np.all(st.y_labeled_iter_ == 0)
