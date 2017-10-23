from sklearn.utils.estimator_checks import check_estimator

from ..sfa import SFA


def test_sfa_is_estimator():
    return check_estimator(SFA)
