from dataclasses import dataclass, fields

import pytest

from sklearn.base import (
    BaseEstimator,
    RegressorMixin,
    TransformerMixin,
)
from sklearn.utils import Tags, get_tags
from sklearn.utils._tags import _safe_tags, _to_new_tags, _to_old_tags, default_tags
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


class MixinAllowNanOldTags:
    def _more_tags(self):
        return {"allow_nan": True}


class MixinAllowNanNewTags:
    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.allow_nan = True
        return tags


class MixinAllowNanOldNewTags:
    def _more_tags(self):
        return {"allow_nan": True}  # pragma: no cover

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.allow_nan = True
        return tags


class MixinArrayApiSupportOldTags:
    def _more_tags(self):
        return {"array_api_support": True}


class MixinArrayApiSupportNewTags:
    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        return tags


class MixinArrayApiSupportOldNewTags:
    def _more_tags(self):
        return {"array_api_support": True}  # pragma: no cover

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        return tags


class PredictorOldTags(BaseEstimator):
    def _more_tags(self):
        return {"requires_fit": True}


class PredictorNewTags(BaseEstimator):
    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.requires_fit = True
        return tags


class PredictorOldNewTags(BaseEstimator):
    def _more_tags(self):
        return {"requires_fit": True}  # pragma: no cover

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.requires_fit = True
        return tags


def test_get_tags_backward_compatibility():
    warn_msg = "Please define the `__sklearn_tags__` method"

    ####################################################################################
    # only predictor inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    for predictor_cls in predictor_classes:
        if predictor_cls.__name__.endswith("OldTags"):
            with pytest.warns(FutureWarning, match=warn_msg):
                tags = get_tags(predictor_cls())
        else:
            tags = get_tags(predictor_cls())
        assert tags.requires_fit

    ####################################################################################
    # one mixin and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for allow_nan_cls in allow_nan_classes:
        for predictor_cls in predictor_classes:

            class ChildClass(allow_nan_cls, predictor_cls):
                pass

            if any(
                base_cls.__name__.endswith("OldTags")
                for base_cls in (predictor_cls, allow_nan_cls)
            ):
                with pytest.warns(FutureWarning, match=warn_msg):
                    tags = get_tags(ChildClass())
            else:
                tags = get_tags(ChildClass())

            assert tags.input_tags.allow_nan
            assert tags.requires_fit

    ####################################################################################
    # two mixins and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    array_api_classes = [
        MixinArrayApiSupportNewTags,
        MixinArrayApiSupportOldNewTags,
        MixinArrayApiSupportOldTags,
    ]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for predictor_cls in predictor_classes:
        for array_api_cls in array_api_classes:
            for allow_nan_cls in allow_nan_classes:

                class ChildClass(allow_nan_cls, array_api_cls, predictor_cls):
                    pass

                if any(
                    base_cls.__name__.endswith("OldTags")
                    for base_cls in (predictor_cls, array_api_cls, allow_nan_cls)
                ):
                    with pytest.warns(FutureWarning, match=warn_msg):
                        tags = get_tags(ChildClass())
                else:
                    tags = get_tags(ChildClass())

                assert tags.input_tags.allow_nan
                assert tags.array_api_support
                assert tags.requires_fit


@pytest.mark.filterwarnings(
    "ignore:.*Please define the `__sklearn_tags__` method.*:FutureWarning"
)
def test_safe_tags_backward_compatibility():
    warn_msg = "The `_safe_tags` function is deprecated in 1.6"

    ####################################################################################
    # only predictor inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    for predictor_cls in predictor_classes:
        with pytest.warns(FutureWarning, match=warn_msg):
            tags = _safe_tags(predictor_cls())
        assert tags["requires_fit"]

    ####################################################################################
    # one mixin and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for allow_nan_cls in allow_nan_classes:
        for predictor_cls in predictor_classes:

            class ChildClass(allow_nan_cls, predictor_cls):
                pass

            with pytest.warns(FutureWarning, match=warn_msg):
                tags = _safe_tags(ChildClass())

            assert tags["allow_nan"]
            assert tags["requires_fit"]

    ####################################################################################
    # two mixins and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    array_api_classes = [
        MixinArrayApiSupportNewTags,
        MixinArrayApiSupportOldNewTags,
        MixinArrayApiSupportOldTags,
    ]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for predictor_cls in predictor_classes:
        for array_api_cls in array_api_classes:
            for allow_nan_cls in allow_nan_classes:

                class ChildClass(allow_nan_cls, array_api_cls, predictor_cls):
                    pass

                with pytest.warns(FutureWarning, match=warn_msg):
                    tags = _safe_tags(ChildClass())

                assert tags["allow_nan"]
                assert tags["array_api_support"]
                assert tags["requires_fit"]


@pytest.mark.filterwarnings(
    "ignore:.*Please define the `__sklearn_tags__` method.*:FutureWarning"
)
def test__get_tags_backward_compatibility():
    warn_msg = "The `_get_tags` method is deprecated in 1.6"

    ####################################################################################
    # only predictor inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    for predictor_cls in predictor_classes:
        with pytest.warns(FutureWarning, match=warn_msg):
            tags = predictor_cls()._get_tags()
        assert tags["requires_fit"]

    ####################################################################################
    # one mixin and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for allow_nan_cls in allow_nan_classes:
        for predictor_cls in predictor_classes:

            class ChildClass(allow_nan_cls, predictor_cls):
                pass

            with pytest.warns(FutureWarning, match=warn_msg):
                tags = ChildClass()._get_tags()

            assert tags["allow_nan"]
            assert tags["requires_fit"]

    ####################################################################################
    # two mixins and one predictor all inheriting from BaseEstimator
    predictor_classes = [PredictorNewTags, PredictorOldNewTags, PredictorOldTags]
    array_api_classes = [
        MixinArrayApiSupportNewTags,
        MixinArrayApiSupportOldNewTags,
        MixinArrayApiSupportOldTags,
    ]
    allow_nan_classes = [
        MixinAllowNanNewTags,
        MixinAllowNanOldNewTags,
        MixinAllowNanOldTags,
    ]

    for predictor_cls in predictor_classes:
        for array_api_cls in array_api_classes:
            for allow_nan_cls in allow_nan_classes:

                class ChildClass(allow_nan_cls, array_api_cls, predictor_cls):
                    pass

                with pytest.warns(FutureWarning, match=warn_msg):
                    tags = ChildClass()._get_tags()

                assert tags["allow_nan"]
                assert tags["array_api_support"]
                assert tags["requires_fit"]


def test_roundtrip_tags():
    estimator = PredictorNewTags()
    tags = default_tags(estimator)
    assert _to_new_tags(_to_old_tags(tags), estimator=estimator) == tags


def test_base_estimator_more_tags():
    """Test that the `_more_tags` and `_get_tags` methods are equivalent for
    `BaseEstimator`.
    """
    estimator = BaseEstimator()
    with pytest.warns(FutureWarning, match="The `_more_tags` method is deprecated"):
        more_tags = BaseEstimator._more_tags(estimator)

    with pytest.warns(FutureWarning, match="The `_get_tags` method is deprecated"):
        get_tags = BaseEstimator._get_tags(estimator)

    assert more_tags == get_tags


def test_safe_tags():
    estimator = PredictorNewTags()
    with pytest.warns(FutureWarning, match="The `_safe_tags` function is deprecated"):
        tags = _safe_tags(estimator)

    with pytest.warns(FutureWarning, match="The `_safe_tags` function is deprecated"):
        tags_requires_fit = _safe_tags(estimator, key="requires_fit")

    assert tags_requires_fit == tags["requires_fit"]
