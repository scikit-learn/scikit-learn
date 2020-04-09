# Should be in test_estimator_checks. Putting this here temporarily because I'm
# too lazy to make it work without pytest.

import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import SkipTestWarning
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.estimator_checks import some_strict_check


def test_strict_mode(capsys):

    est = LogisticRegression()
    with pytest.warns(SkipTestWarning, match='is marked as strict'):
        # some_strict_check() was skipped
        check_estimator(est, strict_mode=False)

    # some_strict_check was ran
    check_estimator(est, strict_mode=True)
    assert capsys.readouterr().out == 'IN THE TEST\n'

    # Make sure setting strict mode global var is reverted to its default
    # (True) at the end of check_estimator.
    check_estimator(est, strict_mode=False)
    some_strict_check('name', est)
    assert capsys.readouterr().out == 'IN THE TEST\n'
