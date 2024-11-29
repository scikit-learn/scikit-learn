from dataclasses import dataclass, fields

import pytest

from sklearn.base import (
    BaseEstimator,
    RegressorMixin,
    TransformerMixin,
)
from sklearn.utils import (
    ClassifierTags,
    InputTags,
    RegressorTags,
    Tags,
    TargetTags,
    TransformerTags,
    get_tags,
)
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

    err_msg = "The key unknown_key is not defined"
    with pytest.raises(ValueError, match=err_msg):
        with pytest.warns(
            FutureWarning, match="The `_safe_tags` function is deprecated"
        ):
            _safe_tags(estimator, key="unknown_key")


def test_old_tags():
    """Set to non-default values and check that we get the expected old tags."""

    class MyClass:
        _estimator_type = "regressor"

        def __sklearn_tags__(self):
            input_tags = InputTags(
                one_d_array=True,
                two_d_array=False,
                three_d_array=True,
                sparse=True,
                categorical=True,
                string=True,
                dict=True,
                positive_only=True,
                allow_nan=True,
                pairwise=True,
            )
            target_tags = TargetTags(
                required=False,
                one_d_labels=True,
                two_d_labels=True,
                positive_only=True,
                multi_output=True,
                single_output=False,
            )
            transformer_tags = None
            classifier_tags = None
            regressor_tags = RegressorTags(
                poor_score=True,
            )
            return Tags(
                estimator_type=self._estimator_type,
                input_tags=input_tags,
                target_tags=target_tags,
                transformer_tags=transformer_tags,
                classifier_tags=classifier_tags,
                regressor_tags=regressor_tags,
            )

    estimator = MyClass()
    new_tags = get_tags(estimator)
    old_tags = _to_old_tags(new_tags)
    expected_tags = {
        "allow_nan": True,
        "array_api_support": False,
        "binary_only": False,
        "multilabel": False,
        "multioutput": True,
        "multioutput_only": True,
        "no_validation": False,
        "non_deterministic": False,
        "pairwise": True,
        "preserves_dtype": ["float64"],
        "poor_score": True,
        "requires_fit": True,
        "requires_positive_X": True,
        "requires_y": False,
        "requires_positive_y": True,
        "_skip_test": False,
        "stateless": True,
        "X_types": [
            "1darray",
            "3darray",
            "sparse",
            "categorical",
            "string",
            "dict",
            "1dlabels",
            "2dlabels",
        ],
    }
    assert old_tags == expected_tags
    assert _to_new_tags(_to_old_tags(new_tags), estimator=estimator) == new_tags

    class MyClass:
        _estimator_type = "classifier"

        def __sklearn_tags__(self):
            input_tags = InputTags(
                one_d_array=True,
                two_d_array=False,
                three_d_array=True,
                sparse=True,
                categorical=True,
                string=True,
                dict=True,
                positive_only=True,
                allow_nan=True,
                pairwise=True,
            )
            target_tags = TargetTags(
                required=False,
                one_d_labels=True,
                two_d_labels=False,
                positive_only=True,
                multi_output=True,
                single_output=False,
            )
            transformer_tags = None
            classifier_tags = ClassifierTags(
                poor_score=True,
                multi_class=False,
                multi_label=True,
            )
            regressor_tags = None
            return Tags(
                estimator_type=self._estimator_type,
                input_tags=input_tags,
                target_tags=target_tags,
                transformer_tags=transformer_tags,
                classifier_tags=classifier_tags,
                regressor_tags=regressor_tags,
            )

    estimator = MyClass()
    new_tags = get_tags(estimator)
    old_tags = _to_old_tags(new_tags)
    expected_tags = {
        "allow_nan": True,
        "array_api_support": False,
        "binary_only": True,
        "multilabel": True,
        "multioutput": True,
        "multioutput_only": True,
        "no_validation": False,
        "non_deterministic": False,
        "pairwise": True,
        "preserves_dtype": ["float64"],
        "poor_score": True,
        "requires_fit": True,
        "requires_positive_X": True,
        "requires_y": False,
        "requires_positive_y": True,
        "_skip_test": False,
        "stateless": True,
        "X_types": [
            "1darray",
            "3darray",
            "sparse",
            "categorical",
            "string",
            "dict",
            "1dlabels",
        ],
    }
    assert old_tags == expected_tags
    assert _to_new_tags(_to_old_tags(new_tags), estimator=estimator) == new_tags

    class MyClass:

        def fit(self, X, y=None):
            return self  # pragma: no cover

        def transform(self, X):
            return X  # pragma: no cover

        def __sklearn_tags__(self):
            input_tags = InputTags(
                one_d_array=True,
                two_d_array=False,
                three_d_array=True,
                sparse=True,
                categorical=True,
                string=True,
                dict=True,
                positive_only=True,
                allow_nan=True,
                pairwise=True,
            )
            target_tags = TargetTags(
                required=False,
                one_d_labels=True,
                two_d_labels=False,
                positive_only=True,
                multi_output=True,
                single_output=False,
            )
            transformer_tags = TransformerTags(
                preserves_dtype=["float64"],
            )
            classifier_tags = None
            regressor_tags = None
            return Tags(
                estimator_type=None,
                input_tags=input_tags,
                target_tags=target_tags,
                transformer_tags=transformer_tags,
                classifier_tags=classifier_tags,
                regressor_tags=regressor_tags,
            )

    estimator = MyClass()
    new_tags = get_tags(estimator)
    old_tags = _to_old_tags(new_tags)
    expected_tags = {
        "allow_nan": True,
        "array_api_support": False,
        "binary_only": False,
        "multilabel": False,
        "multioutput": True,
        "multioutput_only": True,
        "no_validation": False,
        "non_deterministic": False,
        "pairwise": True,
        "preserves_dtype": ["float64"],
        "poor_score": False,
        "requires_fit": True,
        "requires_positive_X": True,
        "requires_y": False,
        "requires_positive_y": True,
        "_skip_test": False,
        "stateless": True,
        "X_types": [
            "1darray",
            "3darray",
            "sparse",
            "categorical",
            "string",
            "dict",
            "1dlabels",
        ],
    }
    assert old_tags == expected_tags
    assert _to_new_tags(_to_old_tags(new_tags), estimator=estimator) == new_tags
