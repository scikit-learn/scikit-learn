import types
import warnings

from sklearn.utils import all_estimators


# TODO: remove in 0.24
def test_partial_dependence_no_shadowing():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/15842
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        from sklearn.inspection.partial_dependence import partial_dependence as _  # noqa

        # Calling all_estimators() also triggers a recursive import of all
        # submodules, including deprecated ones.
        all_estimators()

    from sklearn.inspection import partial_dependence
    assert isinstance(partial_dependence, types.FunctionType)


# TODO: remove in 0.24
def test_dict_learning_no_shadowing():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/15842
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        from sklearn.decomposition.dict_learning import dict_learning as _  # noqa

        # Calling all_estimators() also triggers a recursive import of all
        # submodules, including deprecated ones.
        all_estimators()

    from sklearn.decomposition import dict_learning
    assert isinstance(dict_learning, types.FunctionType)
