import pytest
from sklearn.datasets import make_classification
from sklearn.utils._mocking import CheckingClassifier


X, y = make_classification(n_samples=20, random_state=42)


def test_no_checks():
    # check that none of the method calls raises
    clf = CheckingClassifier().fit(X, y)
    clf.predict(X)
    clf.score(X, y)


def test_check_X_on_fit_success():
    def success(X):
        return True

    # does not raise
    CheckingClassifier(check_X=success).fit(X, y)


def test_check_X_on_fit_fail():
    def fail(X):
        return False

    clf = CheckingClassifier(check_X=fail)
    with pytest.raises(AssertionError):
        clf.fit(X, y)


def test_check_y_on_fit_success():
    def success(y):
        return True

    # does not raise
    CheckingClassifier(check_y=success).fit(X, y)


def test_check_y_on_fit_fail():
    def fail(y):
        return False

    clf = CheckingClassifier(check_y=fail)
    with pytest.raises(AssertionError):
        clf.fit(X, y)


def test_check_X_on_predict_success():
    def success(X):
        return True

    clf = CheckingClassifier(check_X=success).fit(X, y)
    clf.predict(X)  # does not raise


def test_check_X_on_predict_fail():
    def fail(X):
        return False

    clf = CheckingClassifier(check_X=fail)
    with pytest.raises(AssertionError):
        clf.predict(X)


def test_expected_fit_params_success():
    fit_params = {'a': y, 'b': y - 1}
    clf = CheckingClassifier(expected_fit_params=fit_params)
    clf.fit(X, y, **fit_params)


def test_expected_fit_params_fails_missing_key():
    # a and b are expected fit_params but only a is passed
    fit_params = {'a': y, 'b': y - 1}
    clf = CheckingClassifier(expected_fit_params=fit_params)
    with pytest.raises(AssertionError) as exc:
        clf.fit(X, y, a=y)

    msg = "Expected fit parameter(s) ['b'] not seen."
    assert exc.value.args[0] == msg


def test_expected_fit_params_fails_wrong_length():
    # fit_param b has only half the expected length
    fit_params = {'a': y, 'b': y - 1}
    clf = CheckingClassifier(expected_fit_params=fit_params)
    with pytest.raises(AssertionError) as exc:
        clf.fit(X, y, a=y, b=y[::2])

    msg = 'Fit parameter b has length 10; expected 20.'
    assert exc.value.args[0] == msg


def test_foo_param_score_1():
    clf = CheckingClassifier(foo_param=2)
    score = clf.score(X, y)
    assert score == 1.


def test_foo_param_score_0():
    clf = CheckingClassifier(foo_param=0)
    score = clf.score(X, y)
    assert score == 0.
