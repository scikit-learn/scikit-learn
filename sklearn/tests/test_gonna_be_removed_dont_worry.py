# Should be in test_estimator_checks. Putting this here temporarily because I'm
# too lazy to make it work without pytest.

# import functools
import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import SkipTestWarning
from sklearn.utils.estimator_checks import check_estimator
# from sklearn.utils.estimator_checks import some_strict_check
# from sklearn.utils.estimator_checks import some_partially_strict_check


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


# class LR1(LogisticRegression):
#     def _more_tags(self):
#         return {
#             '_xfail_checks': {
#                 some_strict_check: 'whatever reason',
#                 some_partially_strict_check: 'whatever reason',
#             }
#         }


# class LR2(LogisticRegression):
#     def _more_tags(self):
#         return {
#             '_xfail_checks': {
#                 functools.partial(some_strict_check): 'whatever reason',
#                 functools.partial(some_partially_strict_check): 'whatever reason',  # noqa
#             }
#         }


# @pytest.mark.parametrize('est, expected_output', (
#     (LR1(), ""),
#     # (LR2(), ""),
#     )
# )
# def test_xfail_tag_is_partial_func(est, expected_output, capsys):

#     check_estimator(est, strict_mode=True)
    # assert capsys.readouterr().out == expected_output
