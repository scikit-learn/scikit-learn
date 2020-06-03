import pprint
import itertools
from inspect import signature
from typing import Optional, List, Dict, Any

from sklearn.datasets import load_iris
from sklearn.base import is_classifier, is_regressor
from sklearn.tree._classes import BaseDecisionTree
from sklearn.ensemble._forest import BaseForest
from sklearn.ensemble._gb import BaseGradientBoosting
from sklearn.utils import all_estimators
from sklearn.utils.estimator_checks import (
    _enforce_estimator_tags_y,
    parametrize_with_checks,
)


categorical_params = [
    "solver",
    "algorithm",
    "loss",
    "strategy",
    "selection",
    "criterion",
    "multi_class",
    "kernel",
]
bool_params = [
    "fit_prior",
    "fit_intercept",
    "positive",
    "normalize",
    "dual",
    "average",
    "shuffle",
]


class FakeParam(str):
    """Fake string parameter

    This object always returns False when compared to other values,
    but remember values it was compared with. It is used to determine
    valid values for a categorical parameter.

    Examples
    --------
    >>> fp = FakeParam()
    >>> isinstance(fp, str)
    True
    >>> fp
    fake-param
    >>> fp == 'a'
    False
    >>> fp.values
    {'a'}
    """

    def __init__(self):
        self.values = set()

    def __repr__(self):
        return "fake-param"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        self.values.add(other)
        return False

    def __getattr__(self, key):
        if key == "requires_vector_input":
            # Workaround for GaussianProcessClassifier
            return False
        raise AttributeError


def detect_valid_categorical_params(
    Estimator, param: str = "solver"
) -> Optional[List]:
    """Detect valid parameters for an estimator


    Returns
    -------
    set or None: a set of valid parameters of None if
                 it they could not be determined (or the
                 estimator has no such parameter name)

    Example
    -------
    >>> from sklearn.linear_model import LogisticRegression
    >>> detect_valid_categorical_params(LogisticRegression, param="solver")
    ['lbfgs', 'liblinear', 'newton-cg', 'sag', 'saga']
    """
    if not (is_regressor(Estimator) or is_classifier(Estimator)):
        return None

    est_signature = signature(Estimator)
    if param not in est_signature.parameters:
        return None

    # Special cases
    if param == "criterion" and issubclass(
        Estimator, (BaseDecisionTree, BaseForest, BaseGradientBoosting)
    ):
        # hardcode this case as the FakeParam apporach then doesn't work with
        # `param in valid_list` checks.
        from sklearn.tree._classes import CRITERIA_CLF, CRITERIA_REG

        if is_regressor(Estimator):
            return list(CRITERIA_REG)
        else:
            return list(CRITERIA_CLF)
    elif param == "loss":
        # hardcode a few other cases that can't be auto-detected
        from sklearn.linear_model import (
            PassiveAggressiveClassifier,
            PassiveAggressiveRegressor,
            SGDClassifier,
            SGDRegressor,
        )

        if issubclass(Estimator, PassiveAggressiveClassifier):
            return ["hinge", "squared_hinge"]
        elif issubclass(Estimator, SGDClassifier):
            return [
                "hinge",
                "log",
                "modified_huber",
                "squared_hinge",
                "perceptron",
            ]
        elif issubclass(Estimator, PassiveAggressiveRegressor):
            return ["epsilon_insensitive", "squared_epsilon_insensitive"]
        elif issubclass(Estimator, SGDRegressor):
            return [
                "squared_loss",
                "huber",
                "epsilon_insensitive",
                "squared_epsilon_insensitive",
            ]
    elif param == "kernel":
        from sklearn.gaussian_process import (
            GaussianProcessClassifier,
            GaussianProcessRegressor,
        )

        if issubclass(
            Estimator, (GaussianProcessClassifier, GaussianProcessRegressor)
        ):
            # these kernels are not string but instances, skip
            return None
    fp = FakeParam()

    X, y = load_iris(return_X_y=True)

    # Auto-detect valid string parameters with FakeParam
    try:
        args = {param: fp}
        est = Estimator(**args)
        y = _enforce_estimator_tags_y(est, y)
        est.fit(X, y)
    except Exception:
        if not fp.values:
            raise
    return list(sorted(fp.values))


def detect_all_params(Estimator) -> Dict[str, List]:
    """Detect all valid parameters for an estimator

    Example
    -------
    >>> from sklearn.linear_model import LogisticRegression
    >>> detect_all_params(LogisticRegression)
    {'solver': ['lbfgs', 'liblinear', 'newton-cg', 'sag', 'saga'],
     'multi_class': ['auto', 'multinomial', 'ovr'],
     'fit_intercept': [False, True],
     'dual': [False, True]}
    """
    res = {}
    for param_name in categorical_params:
        values = detect_valid_categorical_params(Estimator, param=param_name)
        if values is not None:
            res[param_name] = values
    for param_name in bool_params:
        est_signature = signature(Estimator)
        if param_name in est_signature.parameters:
            res[param_name] = [False, True]
    return res


def _merge_dict_product(**kwargs) -> List[Dict]:
    """Merge the cathesian product of dictionaries with lists

    Example
    -------
    >>> _merge_dict_product(a=[1, 2], b=[True, False], c=['O'])
    [{'a': 1, 'b': True, 'c': 'O'},
     {'a': 1, 'b': False, 'c': 'O'},
     {'a': 2, 'b': True, 'c': 'O'},
     {'a': 2, 'b': False, 'c': 'O'}]
    """
    tmp = []
    for key, val in kwargs.items():
        tmp.append([{key: el} for el in val])

    res = []
    for val in itertools.product(*tmp):
        row: Dict[str, Any] = {}
        for el in val:
            row.update(el)
        res.append(row)

    return res


def _make_all_estimator_instances(verbose=False):
    for name, Estimator in all_estimators(
        type_filter=["classifier", "regressor"]
    ):
        valid_params = detect_all_params(Estimator)
        if valid_params:
            for params in _merge_dict_product(**valid_params):
                yield Estimator(**params)
            if verbose:
                print(f"{name}")
                pprint.pp(valid_params, sort_dicts=True)


@parametrize_with_checks(list(_make_all_estimator_instances()))
def test_common_non_default(estimator, check):
    check(estimator)


if __name__ == "__main__":
    # Print the list of tested estimator
    list(_make_all_estimator_instances(verbose=True))
