"""Checks for minimal scikit-learn estimator support."""
from inspect import signature

import numpy as np

# from . import IS_PYPY
from ..base import clone
from ._testing import _get_args, ignore_warnings, set_random_state

from .common_utils_checks import (
    # _enforce_estimator_tags_x,
    _enforce_estimator_tags_y,
    _pairwise_estimator_convert_X,
)


def _yield_api_estimator_checks(estimator):
    # name = estimator.__class__.__name__
    # tags = _safe_tags(estimator)
    # pairwise = _is_pairwise(estimator)

    yield check_no_attributes_set_in_init
    yield check_takes_at_least_optional_y


@ignore_warnings(category=FutureWarning)
def check_no_attributes_set_in_init(name, estimator_orig):
    """Check attribute setting at `__init__`."""
    try:
        # Clone fails if the estimator does not store
        # all parameters as an attribute during init
        estimator = clone(estimator_orig)
    except AttributeError:
        raise AttributeError(
            f"Estimator {name} should store all parameters as an attribute during init."
            " Cloning mechanism will not work otherwise."
        )

    init_params = _get_args(type(estimator).__init__)
    # TODO: check if we can get more generic and not have a special case for PyPy
    # if IS_PYPY:
    #     # __init__ signature has additional objects in PyPy
    #     for key in ["obj"]:
    #         if key in init_params:
    #             init_params.remove(key)
    parents_init_params = [
        param
        for params_parent in (_get_args(parent) for parent in type(estimator).__mro__)
        for param in params_parent
    ]

    # Test for no setting apart from parameters during init
    invalid_attr = set(vars(estimator)) - set(init_params) - set(parents_init_params)
    assert not invalid_attr, (
        f"Estimator {name} should not set any attribute apart"
        f" from parameters during init. Found attributes {sorted(invalid_attr)}."
    )


@ignore_warnings
def check_takes_at_least_optional_y(name, estimator_orig):
    """Check that estimator accepts an optional `y` to be compatible with
    `Pipeline`.

    Tha parameter `y` should be available in the following methods: `fit`,
    `score`, `partial_fit`, `fit_predict`, `fit_transform`.
    """
    estimator = clone(estimator_orig)
    set_random_state(estimator)

    rnd = np.random.RandomState(0)
    n_samples = 30
    X = rnd.uniform(size=(n_samples, 3))
    X = _pairwise_estimator_convert_X(X, estimator_orig)
    y = np.arange(n_samples) % 3
    y = _enforce_estimator_tags_y(estimator, y)

    supported_methods = ["fit", "score", "partial_fit", "fit_predict", "fit_transform"]
    for method_name in supported_methods:
        method = getattr(estimator, method_name, None)
        if method is not None:
            method(X, y)
            args = [p.name for p in signature(method).parameters.values()]
            if args[0] == "self":
                # `if_delegate_has_method` or `available_if` makes methods
                # into functions with an explicit "self", so need to shift
                # arguments
                args = args[1:]
            assert args[1] in ["y", "Y"], (
                "Expected `y` or `Y` as second argument for method "
                f"{method_name}of {name}. Got arguments: {repr(args)}."
            )
