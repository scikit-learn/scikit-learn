# Should be in test_estimator_checks. Putting this here temporarily because I'm
# too lazy to make it work without pytest.

import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import SkipTestWarning
from sklearn.utils.estimator_checks import check_estimator


def test_strict_mode(capsys):

    est = LogisticRegression()
    with pytest.warns(SkipTestWarning, match='strict mode is off'):
        check_estimator(est, strict_mode=False)

    expected_out = 'in non-strict part of some_partially_strict_check\n'
    assert capsys.readouterr().out == expected_out

    check_estimator(est, strict_mode=True)
    expected_output = (
        "in some_strict_check\n"
        "in non-strict part of some_partially_strict_check\n"
        "in strict part of some_partially_strict_check\n"
    )
    assert capsys.readouterr().out == expected_output
