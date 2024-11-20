from dataclasses import dataclass, fields

import pytest

from sklearn.base import (
    BaseEstimator,
    RegressorMixin,
    TransformerMixin,
)
from sklearn.utils import Tags, get_tags
from sklearn.utils._tags import _safe_tags, _to_new_tags
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
            return {"requires_y": True}  # pragma: no cover

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
            as_dict = {
                field.name: getattr(tags_orig, field.name)
                for field in fields(tags_orig)
            }
            tags = MyTags(**as_dict)
            tags.my_tag = True
            return tags

    check_valid_tag_types("MyEstimator", MyEstimator())


########################################################################################
# Test for the deprecation
# TODO(1.7): Remove this
########################################################################################


def test_tags_deprecation():
    class ChildClass(RegressorMixin, BaseEstimator):
        """Child implementing the old tags API together with our new API."""

        def _more_tags(self):
            return {"allow_nan": True}

    main_warn_msg = "only use `_get_tags` and `_more_tags`"
    with pytest.warns(FutureWarning, match=main_warn_msg):
        tags = ChildClass().__sklearn_tags__()
    assert tags.input_tags.allow_nan

    with pytest.warns(FutureWarning) as warning_list:
        tags = _safe_tags(ChildClass())
    assert len(warning_list) == 2, len(warning_list)
    assert str(warning_list[0].message).startswith(
        "The `_safe_tags` utility function is deprecated"
    )
    assert main_warn_msg in str(warning_list[1].message)

    assert isinstance(tags, dict)
    assert _to_new_tags(tags).input_tags.allow_nan

    with pytest.warns(FutureWarning, match=main_warn_msg):
        tags = get_tags(ChildClass())
    assert tags.input_tags.allow_nan

    class ChildClass(RegressorMixin, BaseEstimator):
        """Child implementing the old and new tags API during the transition period."""

        def _more_tags(self):
            return {"allow_nan": True}

        def __sklearn_tags__(self):
            tags = super().__sklearn_tags__()
            tags.input_tags.allow_nan = True
            return tags

    tags = get_tags(ChildClass())
    assert tags.input_tags.allow_nan

    warn_msg = "`_get_tags` tag provider is deprecated"
    with pytest.warns(FutureWarning, match=warn_msg):
        tags = ChildClass()._get_tags()
    assert isinstance(tags, dict)
    assert _to_new_tags(tags).input_tags.allow_nan

    warn_msg = "`_safe_tags` utility function is deprecated"
    with pytest.warns(FutureWarning, match=warn_msg):
        tags = _safe_tags(ChildClass())
    assert isinstance(tags, dict)
    assert _to_new_tags(tags).input_tags.allow_nan

    class ChildClass(RegressorMixin, BaseEstimator):
        """Child not setting any tags."""

    tags = get_tags(ChildClass())
    assert tags.target_tags.required

    warn_msg = "`_get_tags` tag provider is deprecated"
    with pytest.warns(FutureWarning, match=warn_msg):
        tags = ChildClass()._get_tags()
    assert isinstance(tags, dict)
    assert _to_new_tags(tags).target_tags.required

    warn_msg = "`_safe_tags` utility function is deprecated"
    with pytest.warns(FutureWarning, match=warn_msg):
        tags = _safe_tags(ChildClass())
    assert isinstance(tags, dict)
    assert _to_new_tags(tags).target_tags.required

    class Mixin:
        def _more_tags(self):
            return {"allow_nan": True}

    class ChildClass(Mixin, BaseEstimator):
        """Child following the new API with mixin following the old API."""

        def __sklearn_tags__(self):
            tags = super().__sklearn_tags__()
            tags.target_tags.required = True
            return tags

    err_msg = (
        "Some classes from which ChildClass inherits only use `_get_tags` and "
        "`_more_tags`"
    )
    with pytest.raises(ValueError, match=err_msg):
        tags = get_tags(ChildClass())
    with pytest.raises(ValueError, match=err_msg):
        with pytest.warns(FutureWarning):
            tags = ChildClass()._get_tags()
    with pytest.raises(ValueError, match=err_msg):
        with pytest.warns(FutureWarning):
            tags = _safe_tags(ChildClass())

    class Mixin:
        def _more_tags(self):
            return {"allow_nan": True}

    class ChildClass(Mixin, BaseEstimator):
        """Child following the old API with mixin following the old API."""

        def _more_tags(self):
            return {"requires_y": True}

    with pytest.warns(FutureWarning, match=main_warn_msg):
        tags = ChildClass().__sklearn_tags__()
    assert tags.input_tags.allow_nan

    with pytest.warns(FutureWarning) as warning_list:
        tags = _safe_tags(ChildClass())
    assert len(warning_list) == 2, len(warning_list)
    assert str(warning_list[0].message).startswith(
        "The `_safe_tags` utility function is deprecated"
    )
    assert main_warn_msg in str(warning_list[1].message)

    assert isinstance(tags, dict)
    assert _to_new_tags(tags).input_tags.allow_nan
