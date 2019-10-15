from numpy.testing import assert_array_equal
import pytest
from sklearn.datasets import make_classification
from sklearn.utils._mocking import CheckingClassifier


X, y = make_classification(n_samples=20, random_state=42)


def _success(x):
    return True


def _fail(x):
    return False


@pytest.mark.parametrize('kwargs', [
    {},
    {'check_X': _success},
    {'check_y': _success},
    {'check_X': _success, 'check_y': _success},
])
def test_check_on_fit_success(kwargs):
    # does not raise
    CheckingClassifier(**kwargs).fit(X, y)


@pytest.mark.parametrize('kwargs', [
    {'check_X': _fail},
    {'check_y': _fail},
    {'check_X': _success, 'check_y': _fail},
    {'check_X': _fail, 'check_y': _success},
    {'check_X': _fail, 'check_y': _fail},
])
def test_check_on_fit_fail(kwargs):
    clf = CheckingClassifier(**kwargs)
    with pytest.raises(AssertionError):
        clf.fit(X, y)


def test_check_X_on_predict_success():
    clf = CheckingClassifier(check_X=_success).fit(X, y)
    clf.predict(X)  # does not raise


def test_check_X_on_predict_fail():
    clf = CheckingClassifier(check_X=_fail)
    with pytest.raises(AssertionError):
        clf.predict(X)


def test_classes_set_in_fit():
    clf = CheckingClassifier().fit(X, y)
    assert_array_equal(clf.classes_, [0, 1])


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


@pytest.mark.parametrize('foo_param, expected', [
    (2, 1.),
    (0, 0.),
])
def test_foo_param_score(foo_param, expected):
    clf = CheckingClassifier(foo_param=foo_param)
    score = clf.score(X, y)
    assert score == pytest.approx(expected)
