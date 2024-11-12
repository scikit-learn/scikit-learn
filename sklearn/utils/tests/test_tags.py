from dataclasses import asdict, dataclass

import pytest

from sklearn.base import (
    BaseEstimator,
    RegressorMixin,
    TransformerMixin,
)
from sklearn.utils import Tags, get_tags
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


@pytest.mark.filterwarnings("ignore:.*no __sklearn_tags__ attribute.*:FutureWarning")
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
            return {"requires_y": True}

    with pytest.raises(
        TypeError, match="has defined either `_more_tags` or `_get_tags`"
    ):
        check_estimator_tags_renamed("MoreTagsEstimator", MoreTagsEstimator())


def test_tag_test_passes_with_inheritance():
    @dataclass
    class MyTags(Tags):
        my_tag: bool = True

    class MyEstimator(BaseEstimator):
        def __sklearn_tags__(self):
            tags_orig = super().__sklearn_tags__()
            tags = MyTags(**asdict(tags_orig))
            tags.input_tags = tags_orig.input_tags
            tags.target_tags = tags_orig.target_tags
            tags.classifier_tags = tags_orig.classifier_tags
            tags.regressor_tags = tags_orig.regressor_tags
            tags.transformer_tags = tags_orig.transformer_tags
            tags.my_tag = True
            return tags

    check_valid_tag_types("MyEstimator", MyEstimator())
