import pprint
import itertools
from inspect import signature
from typing import Optional, List, Dict, Any

from sklearn.base import is_regressor
from sklearn.datasets import load_iris
from sklearn.tree._classes import BaseDecisionTree
from sklearn.ensemble._forest import BaseForest
from sklearn.ensemble._gb import BaseGradientBoosting
from sklearn.utils import all_estimators
from sklearn.utils.estimator_checks import (
    _enforce_estimator_tags_y,
    parametrize_with_checks,
)
from sklearn.utils._testing import ignore_warnings


categorical_params = [
    "solver",
    "algorithm",
    "loss",
    "strategy",
    "selection",
    "criterion",
    "multi_class",
    "kernel",
    "affinity",
    "linkage",
    "metric",
    "init",
    "eigen_solver",
    "initial_strategy",
    "imputation_order",
    "encode",
    "learning_method",
    "method",
    "fit_algorithm",
    "norm",
    "svd_solver",
    "order",
    "mode",
    "output_distribution",
]
bool_params = [
    "fit_prior",
    "fit_intercept",
    "positive",
    "normalize",
    "dual",
    "average",
    "shuffle",
    "whiten",
    "path_method",
    "include_bias",
    "interaction_only",
    "standardize",
    "with_centering",
    "with_scaling",
    "with_mean",
    "with_std",
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
    name = Estimator.__name__
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
            SGDClassifier,
            SGDRegressor,
        )

        if name == "PassiveAggressiveClassifier":
            return ["hinge", "squared_hinge"]
        elif issubclass(Estimator, SGDClassifier):
            return [
                "hinge",
                "log",
                "modified_huber",
                "squared_hinge",
                "perceptron",
            ]
        elif name == "PassiveAggressiveRegressor":
            return ["epsilon_insensitive", "squared_epsilon_insensitive"]
        elif issubclass(Estimator, SGDRegressor):
            return [
                "squared_loss",
                "huber",
                "epsilon_insensitive",
                "squared_epsilon_insensitive",
            ]
    elif param == "kernel":
        if name in [
            "GaussianProcessClassifier",
            "GaussianProcessRegressor",
            "Nystroem",
        ]:
            # these kernels are not string but instances, skip
            return None
        elif name == "KernelPCA":
            return ["linear", "poly", "rbf", "sigmoid", "cosine"]
    elif param == "affinity":
        if name in ["AgglomerativeClustering", "FeatureAgglomeration"]:
            return ["euclidean", "l1", "l2", "manhattan", "cosine"]
    elif param == "linkage":
        if name in ["AgglomerativeClustering", "FeatureAgglomeration"]:
            return ["ward", "complete", "average", "single"]
    elif param == "norm":
        if name in ["HashingVectorizer", "TfidfTransformer"]:
            # Vectorizers are not suppored in common tests for now
            return None
        elif name == "Normalizer":
            return ["l1", "l2", "max"]
        elif name == "ComplementNB":
            return [True, False]
    elif param == "order":
        if name == "PolynomialFeatures":
            return ["F", "C"]
    elif param == "mode":
        if name == "GenericUnivariateSelect":
            return ["percentile", "k_best", "fpr", "fdr", "fwe"]
        elif name in ["KNeighborsTransformer", "RadiusNeighborsTransformer"]:
            return ["distance", "connectivity"]

    fp = FakeParam()

    X, y = load_iris(return_X_y=True)

    # Auto-detect valid string parameters with FakeParam
    try:
        args = {param: fp}
        est = Estimator(**args)
        y = _enforce_estimator_tags_y(est, y)
        with ignore_warnings():
            est.fit(X, y)
    except Exception:
        if not fp.values:
            raise
    if not fp.values:
        raise ValueError(
            f"{Estimator.__name__}: {param}={fp.values} "
            f"should contain at least one element"
        )

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
    res: Dict[str, Any] = {}
    name = Estimator.__name__
    if name in ["ClassifierChain", "RegressorChain"]:
        # skip meta-estimators for now
        return res

    est_signature = signature(Estimator)
    for param_name in est_signature.parameters:
        if param_name in categorical_params:
            values = detect_valid_categorical_params(
                Estimator, param=param_name
            )
            if values is not None:
                res[param_name] = values
        elif param_name in bool_params:
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

    X_raw, y_raw = load_iris(return_X_y=True)

    for name, Estimator in all_estimators(
        type_filter=["transformer", "cluster", "classifier", "regressor"]
    ):
        valid_params = detect_all_params(Estimator)
        if verbose:
            print(f"{name}")
        if valid_params:
            for params in _merge_dict_product(**valid_params):
                # Check that we can train Iris, otherwise parameters
                # are likely incompatible
                try:
                    est = Estimator(**params)
                    y = _enforce_estimator_tags_y(est, y_raw)
                    with ignore_warnings():
                        est.fit(X_raw.copy(), y)
                    # Parameters should be OK
                    yield Estimator(**params)
                except Exception:
                    # Likely wrong parameters, skipping
                    pass

            if verbose:
                pprint.pp(valid_params, sort_dicts=True)


@parametrize_with_checks(list(_make_all_estimator_instances()))
def test_common_non_default(estimator, check):
    check(estimator)


if __name__ == "__main__":
    # Print the list of tested estimator
    list(_make_all_estimator_instances(verbose=True))
