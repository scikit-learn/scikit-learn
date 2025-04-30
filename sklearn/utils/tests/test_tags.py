from dataclasses import dataclass, fields

import numpy as np
import pytest

from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
    RegressorMixin,
    TransformerMixin,
)
from sklearn.pipeline import Pipeline
from sklearn.utils import (
    Tags,
    get_tags,
)
from sklearn.utils.estimator_checks import (
    check_estimator_tags_renamed,
    check_valid_tag_types,
)


class NoTagsEstimator:
    pass


class ClassifierEstimator:
    # This is to test whether not inheriting from mixins works.
    _estimator_type = "classifier"


class EmptyTransformer(TransformerMixin, BaseEstimator):
    pass


class EmptyRegressor(RegressorMixin, BaseEstimator):
    pass


# TODO(1.8): Update when implementing __sklearn_tags__ is required
@pytest.mark.filterwarnings(
    "ignore:.*no attribute '__sklearn_tags__'.*:DeprecationWarning"
)
@pytest.mark.parametrize(
    "estimator, value",
    [
        [NoTagsEstimator(), False],
        [ClassifierEstimator(), True],
        [EmptyTransformer(), False],
        [EmptyRegressor(), True],
        [BaseEstimator(), False],
    ],
)
def test_requires_y(estimator, value):
    assert get_tags(estimator).target_tags.required == value


def test_no___sklearn_tags__with_more_tags():
    """Test that calling `get_tags` on a class that defines `_more_tags` but not
    `__sklearn_tags__` raises an error.
    """

    class MoreTagsEstimator(BaseEstimator):
        def _more_tags(self):
            return {"requires_y": True}  # pragma: no cover

    with pytest.raises(
        TypeError, match="has defined either `_more_tags` or `_get_tags`"
    ):
        check_estimator_tags_renamed("MoreTagsEstimator", MoreTagsEstimator())


def test_tag_test_passes_with_inheritance():
    @dataclass
    class MyTags(Tags):
        my_tag: bool = True  # type: ignore[annotation-unchecked]

    class MyEstimator(BaseEstimator):
        def __sklearn_tags__(self):
            tags_orig = super().__sklearn_tags__()
            as_dict = {
                field.name: getattr(tags_orig, field.name)
                for field in fields(tags_orig)
            }
            tags = MyTags(**as_dict)
            tags.my_tag = True
            return tags

    check_valid_tag_types("MyEstimator", MyEstimator())


# TODO(1.8): Update this test to check for errors
def test_tags_no_sklearn_tags_concrete_implementation():
    """Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/30479

    Either the estimator doesn't implement `__sklearn_tags` or there is no class
    implementing `__sklearn_tags__` without calling `super().__sklearn_tags__()` in
    its mro. Thus, we raise a warning and request to inherit from
    `BaseEstimator` that implements `__sklearn_tags__`.
    """

    X = np.array([[1, 2], [2, 3], [3, 4]])
    y = np.array([1, 0, 1])

    # 1st case, the estimator inherits from a class that only implements
    # `__sklearn_tags__` by calling `super().__sklearn_tags__()`.
    class MyEstimator(ClassifierMixin):
        def __init__(self, *, param=1):
            self.param = param

        def fit(self, X, y=None):
            self.is_fitted_ = True
            return self

        def predict(self, X):
            return np.full(shape=X.shape[0], fill_value=self.param)

    my_pipeline = Pipeline([("estimator", MyEstimator(param=1))])
    with pytest.warns(DeprecationWarning, match="The following error was raised"):
        my_pipeline.fit(X, y).predict(X)

    # 2nd case, the estimator doesn't implement `__sklearn_tags__` at all.
    class MyEstimator2:
        def __init__(self, *, param=1):
            self.param = param

        def fit(self, X, y=None):
            self.is_fitted_ = True
            return self

        def predict(self, X):
            return np.full(shape=X.shape[0], fill_value=self.param)

    my_pipeline = Pipeline([("estimator", MyEstimator2(param=1))])
    with pytest.warns(DeprecationWarning, match="The following error was raised"):
        my_pipeline.fit(X, y).predict(X)

    # check that we still raise an error if it is not a AttributeError or related to
    # __sklearn_tags__
    class MyEstimator3(MyEstimator, BaseEstimator):
        def __init__(self, *, param=1, error_type=AttributeError):
            self.param = param
            self.error_type = error_type

        def __sklearn_tags__(self):
            super().__sklearn_tags__()
            raise self.error_type("test")

    for error_type in (AttributeError, TypeError, ValueError):
        estimator = MyEstimator3(param=1, error_type=error_type)
        with pytest.raises(error_type):
            get_tags(estimator)
