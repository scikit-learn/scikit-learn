# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Common tests for estimators that implement callback support.

These checks are meant to run against every estimator that inherits from
`CallbackSupportMixin`. They verify invariants that must hold regardless of
the estimator's internals: lifecycle hook counts, begin/end balance, and
the identity of the estimator passed to each hook.

Estimator-specific checks (e.g. exact hook counts, fitted_estimator quality)
must remain in each estimator's own test file.
"""

import pytest
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning

from sklearn.callback.tests._utils import (
    RecordingCallback,
    skip_callback_test_if_wasm,
)
from sklearn.datasets import make_classification
from sklearn.utils._test_common.instance_generator import (
    _get_check_estimator_ids,
    _tested_estimators,
)
from sklearn.utils._testing import ignore_warnings

# Exhaustive list of estimators that do NOT yet support callbacks.
# When an estimator gains callback support, remove it from this set; the three
# common checks below will then run against it automatically.
_NO_CALLBACK_SUPPORT = {
    "ARDRegression",
    "AdaBoostClassifier",
    "AdaBoostRegressor",
    "AdditiveChi2Sampler",
    "AffinityPropagation",
    "AgglomerativeClustering",
    "BaggingClassifier",
    "BaggingRegressor",
    "BayesianGaussianMixture",
    "BayesianRidge",
    "BernoulliNB",
    "BernoulliRBM",
    "Binarizer",
    "Birch",
    "BisectingKMeans",
    "CCA",
    "CalibratedClassifierCV",
    "CategoricalNB",
    "ClassicalMDS",
    "ClassifierChain",
    "ColumnTransformer",
    "ComplementNB",
    "CountVectorizer",
    "DBSCAN",
    "DecisionTreeClassifier",
    "DecisionTreeRegressor",
    "DictVectorizer",
    "DictionaryLearning",
    "DummyClassifier",
    "DummyRegressor",
    "ElasticNet",
    "ElasticNetCV",
    "EllipticEnvelope",
    "EmpiricalCovariance",
    "ExtraTreeClassifier",
    "ExtraTreeRegressor",
    "ExtraTreesClassifier",
    "ExtraTreesRegressor",
    "FactorAnalysis",
    "FastICA",
    "FeatureAgglomeration",
    "FeatureHasher",
    "FeatureUnion",
    "FixedThresholdClassifier",
    "FunctionTransformer",
    "GammaRegressor",
    "GaussianMixture",
    "GaussianNB",
    "GaussianProcessClassifier",
    "GaussianProcessRegressor",
    "GaussianRandomProjection",
    "GenericUnivariateSelect",
    "GradientBoostingClassifier",
    "GradientBoostingRegressor",
    "GraphicalLasso",
    "GraphicalLassoCV",
    "HDBSCAN",
    "HashingVectorizer",
    "HistGradientBoostingClassifier",
    "HistGradientBoostingRegressor",
    "HuberRegressor",
    "IncrementalPCA",
    "IsolationForest",
    "Isomap",
    "IsotonicRegression",
    "IterativeImputer",
    "KBinsDiscretizer",
    "KMeans",
    "KNNImputer",
    "KNeighborsClassifier",
    "KNeighborsRegressor",
    "KNeighborsTransformer",
    "KernelCenterer",
    "KernelDensity",
    "KernelPCA",
    "KernelRidge",
    "LabelBinarizer",
    "LabelEncoder",
    "LabelPropagation",
    "LabelSpreading",
    "Lars",
    "LarsCV",
    "Lasso",
    "LassoCV",
    "LassoLars",
    "LassoLarsCV",
    "LassoLarsIC",
    "LatentDirichletAllocation",
    "LedoitWolf",
    "LinearDiscriminantAnalysis",
    "LinearRegression",
    "LinearSVC",
    "LinearSVR",
    "LocalOutlierFactor",
    "LocallyLinearEmbedding",
    "LogisticRegressionCV",
    "MDS",
    "MLPClassifier",
    "MLPRegressor",
    "MaxAbsScaler",
    "MeanShift",
    "MinCovDet",
    "MinMaxScaler",
    "MiniBatchDictionaryLearning",
    "MiniBatchKMeans",
    "MiniBatchNMF",
    "MiniBatchSparsePCA",
    "MissingIndicator",
    "MultiLabelBinarizer",
    "MultiOutputClassifier",
    "MultiOutputRegressor",
    "MultiTaskElasticNet",
    "MultiTaskElasticNetCV",
    "MultiTaskLasso",
    "MultiTaskLassoCV",
    "MultinomialNB",
    "NMF",
    "NearestCentroid",
    "NearestNeighbors",
    "NeighborhoodComponentsAnalysis",
    "Normalizer",
    "NuSVC",
    "NuSVR",
    "Nystroem",
    "OAS",
    "OPTICS",
    "OneClassSVM",
    "OneHotEncoder",
    "OneVsOneClassifier",
    "OneVsRestClassifier",
    "OrdinalEncoder",
    "OrthogonalMatchingPursuit",
    "OrthogonalMatchingPursuitCV",
    "OutputCodeClassifier",
    "PCA",
    "PLSCanonical",
    "PLSRegression",
    "PLSSVD",
    "PassiveAggressiveClassifier",
    "PassiveAggressiveRegressor",
    "PatchExtractor",
    "Perceptron",
    "PoissonRegressor",
    "PolynomialCountSketch",
    "PolynomialFeatures",
    "PowerTransformer",
    "QuadraticDiscriminantAnalysis",
    "QuantileRegressor",
    "QuantileTransformer",
    "RANSACRegressor",
    "RBFSampler",
    "RFE",
    "RFECV",
    "RadiusNeighborsClassifier",
    "RadiusNeighborsRegressor",
    "RadiusNeighborsTransformer",
    "RandomForestClassifier",
    "RandomForestRegressor",
    "RandomTreesEmbedding",
    "RegressorChain",
    "Ridge",
    "RidgeCV",
    "RidgeClassifier",
    "RidgeClassifierCV",
    "RobustScaler",
    "SGDClassifier",
    "SGDOneClassSVM",
    "SGDRegressor",
    "SVC",
    "SVR",
    "SelectFdr",
    "SelectFpr",
    "SelectFromModel",
    "SelectFwe",
    "SelectKBest",
    "SelectPercentile",
    "SelfTrainingClassifier",
    "SequentialFeatureSelector",
    "ShrunkCovariance",
    "SimpleImputer",
    "SkewedChi2Sampler",
    "SparseCoder",
    "SparsePCA",
    "SparseRandomProjection",
    "SpectralBiclustering",
    "SpectralClustering",
    "SpectralCoclustering",
    "SpectralEmbedding",
    "SplineTransformer",
    "StackingClassifier",
    "StackingRegressor",
    "TSNE",
    "TargetEncoder",
    "TfidfTransformer",
    "TfidfVectorizer",
    "TheilSenRegressor",
    "TransformedTargetRegressor",
    "TruncatedSVD",
    "TunedThresholdClassifierCV",
    "TweedieRegressor",
    "VarianceThreshold",
    "VotingClassifier",
    "VotingRegressor",
}


def _has_callback_support(estimator):
    """Return True if this estimator instance has working callback support.

    Class-level exclusions are declared in `_NO_CALLBACK_SUPPORT`. Add
    instance-level conditions below for estimators whose callback support
    depends on their parameters.
    """
    if type(estimator).__name__ in _NO_CALLBACK_SUPPORT:
        return False
    # LogisticRegression only supports callbacks with the lbfgs solver.
    if (
        type(estimator).__name__ == "LogisticRegression"
        and estimator.solver != "lbfgs"
    ):
        return False
    return True


_CALLBACK_ESTIMATORS = [
    est for est in _tested_estimators() if _has_callback_support(est)
]


def _make_data(estimator):
    """`make_classification` is used instead of `make_blobs` because it
    guarantees class balance by construction, which is necessary for
    HalvingSearch whose first rounds use very few samples.
    """
    return make_classification(
        n_samples=200, n_features=4, n_informative=4, n_redundant=0, random_state=0
    )


@ignore_warnings(category=(ConvergenceWarning, UserWarning))
def check_callback_setup_teardown_called_once(estimator, X, y):
    """setup and teardown are each called exactly once per fit, in that order.

    This verifies that the estimator correctly wraps its fit method with
    `@with_callbacks` or `callback_management_context`, which guarantees the
    lifecycle hooks are called exactly once regardless of what happens inside
    fit.
    """
    name = type(estimator).__name__
    cb = RecordingCallback()
    clone(estimator).set_callbacks(cb).fit(X, y)

    msg = f"{name}: expected setup to be called once, got {cb.count_hooks('setup')}"
    assert cb.count_hooks("setup") == 1, msg

    msg = f"{name}: expected teardown to be called once, got {cb.count_hooks('teardown')}"
    assert cb.count_hooks("teardown") == 1, msg

    hook_names = [entry["name"] for entry in cb.record]
    msg = f"{name}: teardown was recorded before setup"
    assert hook_names.index("setup") < hook_names.index("teardown"), msg


@ignore_warnings(category=(ConvergenceWarning, UserWarning))
def check_callback_begin_end_balanced(estimator, X, y):
    """on_fit_task_begin / on_fit_task_end calls form a single-rooted well-parenthesized sequence.

    Because a single `fit` call corresponds to exactly one root task, the
    sequence of task events must be a *primitive* Dyck word: the running
    balance (n_begin - n_end so far) must stay ≥ 1 for every event except the
    last, and be 0 only at the very end. This rules out both unmatched ends
    (balance goes negative) and multiple sibling root tasks, e.g. `()()`.
    """
    name = type(estimator).__name__
    cb = RecordingCallback()
    clone(estimator).set_callbacks(cb).fit(X, y)

    task_events = [
        entry
        for entry in cb.record
        if entry["name"] in ("on_fit_task_begin", "on_fit_task_end")
    ]

    balance = 0
    for i, entry in enumerate(task_events):
        if entry["name"] == "on_fit_task_begin":
            balance += 1
        else:
            balance -= 1
            msg = f"{name}: on_fit_task_end called without a matching on_fit_task_begin"
            assert balance >= 0, msg
            if balance == 0:
                msg = (
                    f"{name}: balance returned to 0 before the last task event; "
                    f"multiple top-level tasks detected"
                )
                assert i == len(task_events) - 1, msg
    msg = f"{name}: {balance} on_fit_task_begin call(s) have no matching on_fit_task_end"
    assert balance == 0, msg


@ignore_warnings(category=(ConvergenceWarning, UserWarning))
def check_callback_estimator_is_self(estimator, X, y):
    """Every hook receives the estimator instance that fit was called on.

    Each call to setup, on_fit_task_begin, on_fit_task_end, and teardown must
    pass the estimator that owns the fit (i.e. `self`), not a clone, not a
    sub-estimator. This is currently only tested for SearchCV estimators; this
    check extends the same invariant to all callback-supporting estimators.
    """
    cb = RecordingCallback()
    est = clone(estimator).set_callbacks(cb)
    est.fit(X, y)

    for entry in cb.record:
        msg = (
            f"{type(estimator).__name__}: hook '{entry['name']}' received "
            f"{entry['estimator']!r} as estimator, expected the fitted instance {est!r}"
        )
        assert entry["estimator"] is est, msg


_CHECKS = [
    check_callback_setup_teardown_called_once,
    check_callback_begin_end_balanced,
    check_callback_estimator_is_self,
]




@pytest.mark.parametrize("check", _CHECKS, ids=lambda fn: fn.__name__)
@pytest.mark.parametrize("estimator", _CALLBACK_ESTIMATORS, ids=_get_check_estimator_ids)
@skip_callback_test_if_wasm
def test_callback(estimator, check):
    X, y = _make_data(estimator)
    check(estimator, X, y)
