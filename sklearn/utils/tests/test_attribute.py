import pytest

from sklearn.utils._attribute_validation import FittedAttribute


class MyObj:
    coef_ = FittedAttribute()
    cv_results_ = FittedAttribute("cv=None")

    def fit(self):
        self.coef_ = 123
        self.cv_results_ = [1, 2, 3]


def test_attribute_error():
    my_obj = MyObj()

    msg = "Call `fit` to define 'coef_'"
    with pytest.raises(AttributeError, match=msg):
        my_obj.coef_

    msg = "Call `fit` with cv=None to define 'cv_results_'"
    with pytest.raises(AttributeError, match=msg):
        my_obj.cv_results_

    my_obj.fit()
    assert my_obj.coef_ == 123
    assert my_obj.cv_results_ == [1, 2, 3]
