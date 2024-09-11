# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import re
import warnings
from functools import partial
from inspect import isfunction

from sklearn import clone, config_context
from sklearn.calibration import CalibratedClassifierCV
from sklearn.cluster import (
    HDBSCAN,
    AffinityPropagation,
    AgglomerativeClustering,
    Birch,
    BisectingKMeans,
    FeatureAgglomeration,
    KMeans,
    MeanShift,
    MiniBatchKMeans,
    SpectralBiclustering,
    SpectralClustering,
    SpectralCoclustering,
)
from sklearn.compose import ColumnTransformer
from sklearn.covariance import GraphicalLasso, GraphicalLassoCV
from sklearn.cross_decomposition import CCA, PLSSVD, PLSCanonical, PLSRegression
from sklearn.decomposition import (
    NMF,
    PCA,
    DictionaryLearning,
    FactorAnalysis,
    FastICA,
    IncrementalPCA,
    KernelPCA,
    LatentDirichletAllocation,
    MiniBatchDictionaryLearning,
    MiniBatchNMF,
    MiniBatchSparsePCA,
    SparseCoder,
    SparsePCA,
    TruncatedSVD,
)
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import (
    AdaBoostClassifier,
    AdaBoostRegressor,
    BaggingClassifier,
    BaggingRegressor,
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
    IsolationForest,
    RandomForestClassifier,
    RandomForestRegressor,
    RandomTreesEmbedding,
    StackingClassifier,
    StackingRegressor,
    VotingClassifier,
    VotingRegressor,
)
from sklearn.exceptions import SkipTestWarning
from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.feature_selection import (
    RFE,
    RFECV,
    SelectFdr,
    SelectFromModel,
    SelectKBest,
    SequentialFeatureSelector,
)
from sklearn.kernel_approximation import (
    Nystroem,
    PolynomialCountSketch,
    RBFSampler,
    SkewedChi2Sampler,
)
from sklearn.linear_model import (
    ARDRegression,
    BayesianRidge,
    ElasticNet,
    ElasticNetCV,
    GammaRegressor,
    HuberRegressor,
    LarsCV,
    Lasso,
    LassoCV,
    LassoLars,
    LassoLarsCV,
    LassoLarsIC,
    LinearRegression,
    LogisticRegression,
    LogisticRegressionCV,
    MultiTaskElasticNet,
    MultiTaskElasticNetCV,
    MultiTaskLasso,
    MultiTaskLassoCV,
    OrthogonalMatchingPursuitCV,
    PassiveAggressiveClassifier,
    PassiveAggressiveRegressor,
    Perceptron,
    PoissonRegressor,
    RANSACRegressor,
    Ridge,
    SGDClassifier,
    SGDOneClassSVM,
    SGDRegressor,
    TheilSenRegressor,
    TweedieRegressor,
)
from sklearn.manifold import (
    MDS,
    TSNE,
    Isomap,
    LocallyLinearEmbedding,
    SpectralEmbedding,
)
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture
from sklearn.model_selection import (
    FixedThresholdClassifier,
    GridSearchCV,
    HalvingGridSearchCV,
    HalvingRandomSearchCV,
    RandomizedSearchCV,
    TunedThresholdClassifierCV,
)
from sklearn.multiclass import (
    OneVsOneClassifier,
    OneVsRestClassifier,
    OutputCodeClassifier,
)
from sklearn.multioutput import (
    ClassifierChain,
    MultiOutputClassifier,
    MultiOutputRegressor,
    RegressorChain,
)
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.neural_network import BernoulliRBM, MLPClassifier, MLPRegressor
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler, TargetEncoder
from sklearn.random_projection import (
    GaussianRandomProjection,
    SparseRandomProjection,
)
from sklearn.semi_supervised import (
    LabelPropagation,
    LabelSpreading,
    SelfTrainingClassifier,
)
from sklearn.svm import SVC, SVR, LinearSVC, LinearSVR, NuSVC, NuSVR, OneClassSVM
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils import all_estimators
from sklearn.utils._testing import SkipTest

CROSS_DECOMPOSITION = ["PLSCanonical", "PLSRegression", "CCA", "PLSSVD"]

# The following dictionary is to indicate constructor arguments suitable for the test
# suite, which uses very small datasets, and is intended to run rather quickly.
INIT_PARAMS = {
    AdaBoostClassifier: dict(n_estimators=5),
    AdaBoostRegressor: dict(n_estimators=5),
    AffinityPropagation: dict(max_iter=5),
    AgglomerativeClustering: dict(n_clusters=2),
    ARDRegression: dict(max_iter=5),
    BaggingClassifier: dict(n_estimators=5),
    BaggingRegressor: dict(n_estimators=5),
    BayesianGaussianMixture: dict(n_init=2, max_iter=5),
    BayesianRidge: dict(max_iter=5),
    BernoulliRBM: dict(n_iter=5, batch_size=10),
    Birch: dict(n_clusters=2),
    BisectingKMeans: dict(n_init=2, n_clusters=2, max_iter=5),
    CalibratedClassifierCV: dict(estimator=LogisticRegression(C=1), cv=3),
    CCA: dict(n_components=1, max_iter=5),
    ClassifierChain: dict(base_estimator=LogisticRegression(C=1), cv=3),
    ColumnTransformer: dict(transformers=[("trans1", StandardScaler(), [0, 1])]),
    DictionaryLearning: dict(max_iter=20, transform_algorithm="lasso_lars"),
    # the default strategy prior would output constant predictions and fail
    # for check_classifiers_predictions
    DummyClassifier: dict(strategy="stratified"),
    ElasticNetCV: dict(max_iter=5, cv=3),
    ElasticNet: dict(max_iter=5),
    ExtraTreesClassifier: dict(n_estimators=5),
    ExtraTreesRegressor: dict(n_estimators=5),
    FactorAnalysis: dict(max_iter=5),
    FastICA: dict(max_iter=5),
    FeatureAgglomeration: dict(n_clusters=2),
    FeatureUnion: dict(transformer_list=[("trans1", StandardScaler())]),
    FixedThresholdClassifier: dict(estimator=LogisticRegression(C=1)),
    GammaRegressor: dict(max_iter=5),
    GaussianMixture: dict(n_init=2, max_iter=5),
    # Due to the jl lemma and often very few samples, the number
    # of components of the random matrix projection will be probably
    # greater than the number of features.
    # So we impose a smaller number (avoid "auto" mode)
    GaussianRandomProjection: dict(n_components=2),
    GradientBoostingClassifier: dict(n_estimators=5),
    GradientBoostingRegressor: dict(n_estimators=5),
    GraphicalLassoCV: dict(max_iter=5, cv=3),
    GraphicalLasso: dict(max_iter=5),
    GridSearchCV: [
        dict(
            cv=2,
            error_score="raise",
            estimator=Ridge(),
            param_grid={"alpha": [0.1, 1.0]},
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=LogisticRegression(),
            param_grid={"C": [0.1, 1.0]},
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(steps=[("pca", PCA()), ("ridge", Ridge())]),
            param_grid={"ridge__alpha": [0.1, 1.0]},
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(
                steps=[("pca", PCA()), ("logisticregression", LogisticRegression())]
            ),
            param_grid={"logisticregression__C": [0.1, 1.0]},
        ),
    ],
    HalvingGridSearchCV: [
        dict(
            cv=2,
            error_score="raise",
            estimator=Ridge(),
            min_resources="smallest",
            param_grid={"alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=LogisticRegression(),
            min_resources="smallest",
            param_grid={"C": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(steps=[("pca", PCA()), ("ridge", Ridge())]),
            min_resources="smallest",
            param_grid={"ridge__alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(
                steps=[("pca", PCA()), ("logisticregression", LogisticRegression())]
            ),
            min_resources="smallest",
            param_grid={"logisticregression__C": [0.1, 1.0]},
            random_state=0,
        ),
    ],
    HalvingRandomSearchCV: [
        dict(
            cv=2,
            error_score="raise",
            estimator=Ridge(),
            param_distributions={"alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=LogisticRegression(),
            param_distributions={"C": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(steps=[("pca", PCA()), ("ridge", Ridge())]),
            param_distributions={"ridge__alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(
                steps=[("pca", PCA()), ("logisticregression", LogisticRegression())]
            ),
            param_distributions={"logisticregression__C": [0.1, 1.0]},
            random_state=0,
        ),
    ],
    HDBSCAN: dict(min_samples=1),
    # The default min_samples_leaf (20) isn't appropriate for small
    # datasets (only very shallow trees are built) that the checks use.
    HistGradientBoostingClassifier: dict(max_iter=5, min_samples_leaf=5),
    HistGradientBoostingRegressor: dict(max_iter=5, min_samples_leaf=5),
    HuberRegressor: dict(max_iter=5),
    IncrementalPCA: dict(batch_size=10),
    IsolationForest: dict(n_estimators=5),
    KMeans: dict(n_init=2, n_clusters=2, max_iter=5),
    LabelPropagation: dict(max_iter=5),
    LabelSpreading: dict(max_iter=5),
    LarsCV: dict(max_iter=5, cv=3),
    LassoCV: dict(max_iter=5, cv=3),
    Lasso: dict(max_iter=5),
    LassoLarsCV: dict(max_iter=5, cv=3),
    LassoLars: dict(max_iter=5),
    # Noise variance estimation does not work when `n_samples < n_features`.
    # We need to provide the noise variance explicitly.
    LassoLarsIC: dict(max_iter=5, noise_variance=1.0),
    LatentDirichletAllocation: dict(max_iter=5, batch_size=10),
    LinearSVC: dict(max_iter=20),
    LinearSVR: dict(max_iter=20),
    LocallyLinearEmbedding: dict(max_iter=5),
    LogisticRegressionCV: dict(max_iter=5, cv=3),
    LogisticRegression: dict(max_iter=5),
    MDS: dict(n_init=2, max_iter=5),
    # In the case of check_fit2d_1sample, bandwidth is set to None and
    # is thus estimated. De facto it is 0.0 as a single sample is provided
    # and this makes the test fails. Hence we give it a placeholder value.
    MeanShift: dict(max_iter=5, bandwidth=1.0),
    MiniBatchDictionaryLearning: dict(batch_size=10, max_iter=5),
    MiniBatchKMeans: dict(n_init=2, n_clusters=2, max_iter=5, batch_size=10),
    MiniBatchNMF: dict(batch_size=10, max_iter=20, fresh_restarts=True),
    MiniBatchSparsePCA: dict(max_iter=5, batch_size=10),
    MLPClassifier: dict(max_iter=100),
    MLPRegressor: dict(max_iter=100),
    MultiOutputClassifier: dict(estimator=LogisticRegression(C=1)),
    MultiOutputRegressor: dict(estimator=Ridge()),
    MultiTaskElasticNetCV: dict(max_iter=5, cv=3),
    MultiTaskElasticNet: dict(max_iter=5),
    MultiTaskLassoCV: dict(max_iter=5, cv=3),
    MultiTaskLasso: dict(max_iter=5),
    NeighborhoodComponentsAnalysis: dict(max_iter=5),
    NMF: dict(max_iter=500),
    NuSVC: dict(max_iter=-1),
    NuSVR: dict(max_iter=-1),
    OneClassSVM: dict(max_iter=-1),
    OneHotEncoder: dict(handle_unknown="ignore"),
    OneVsOneClassifier: dict(estimator=LogisticRegression(C=1)),
    OneVsRestClassifier: dict(estimator=LogisticRegression(C=1)),
    OrthogonalMatchingPursuitCV: dict(cv=3),
    OutputCodeClassifier: dict(estimator=LogisticRegression(C=1)),
    PassiveAggressiveClassifier: dict(max_iter=5),
    PassiveAggressiveRegressor: dict(max_iter=5),
    Perceptron: dict(max_iter=5),
    Pipeline: [
        {"steps": [("scaler", StandardScaler()), ("final_estimator", Ridge())]},
        {
            "steps": [
                ("scaler", StandardScaler()),
                ("final_estimator", LogisticRegression()),
            ]
        },
    ],
    PLSCanonical: dict(n_components=1, max_iter=5),
    PLSRegression: dict(n_components=1, max_iter=5),
    PLSSVD: dict(n_components=1),
    PoissonRegressor: dict(max_iter=5),
    RandomForestClassifier: dict(n_estimators=5),
    RandomForestRegressor: dict(n_estimators=5),
    RandomizedSearchCV: [
        dict(
            cv=2,
            error_score="raise",
            estimator=Ridge(),
            param_distributions={"alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=LogisticRegression(),
            param_distributions={"C": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(steps=[("pca", PCA()), ("ridge", Ridge())]),
            param_distributions={"ridge__alpha": [0.1, 1.0]},
            random_state=0,
        ),
        dict(
            cv=2,
            error_score="raise",
            estimator=Pipeline(
                steps=[("pca", PCA()), ("logisticregression", LogisticRegression())]
            ),
            param_distributions={"logisticregression__C": [0.1, 1.0]},
            random_state=0,
        ),
    ],
    RandomTreesEmbedding: dict(n_estimators=5),
    # `RANSACRegressor` will raise an error with any model other
    # than `LinearRegression` if we don't fix the `min_samples` parameter.
    # For common tests, we can enforce using `LinearRegression` that
    # is the default estimator in `RANSACRegressor` instead of `Ridge`.
    RANSACRegressor: dict(estimator=LinearRegression(), max_trials=10),
    RegressorChain: dict(base_estimator=Ridge(), cv=3),
    RFECV: dict(estimator=LogisticRegression(C=1), cv=3),
    RFE: dict(estimator=LogisticRegression(C=1)),
    # be tolerant of noisy datasets (not actually speed)
    SelectFdr: dict(alpha=0.5),
    # Increases coverage because SGDRegressor has partial_fit
    SelectFromModel: dict(estimator=SGDRegressor(random_state=0)),
    # SelectKBest has a default of k=10
    # which is more feature than we have in most case.
    SelectKBest: dict(k=1),
    SelfTrainingClassifier: dict(estimator=LogisticRegression(C=1), max_iter=5),
    SequentialFeatureSelector: dict(estimator=LogisticRegression(C=1), cv=3),
    SGDClassifier: dict(max_iter=5),
    SGDOneClassSVM: dict(max_iter=5),
    SGDRegressor: dict(max_iter=5),
    SparsePCA: dict(max_iter=5),
    # Due to the jl lemma and often very few samples, the number
    # of components of the random matrix projection will be probably
    # greater than the number of features.
    # So we impose a smaller number (avoid "auto" mode)
    SparseRandomProjection: dict(n_components=2),
    SpectralBiclustering: dict(n_init=2, n_best=1, n_clusters=2),
    SpectralClustering: dict(n_init=2, n_clusters=2),
    SpectralCoclustering: dict(n_init=2, n_clusters=2),
    # Default "auto" parameter can lead to different ordering of eigenvalues on
    # windows: #24105
    SpectralEmbedding: dict(eigen_tol=1e-05),
    StackingClassifier: dict(
        estimators=[
            ("est1", DecisionTreeClassifier(max_depth=3, random_state=0)),
            ("est2", DecisionTreeClassifier(max_depth=3, random_state=1)),
        ],
        cv=3,
    ),
    StackingRegressor: dict(
        estimators=[
            ("est1", DecisionTreeRegressor(max_depth=3, random_state=0)),
            ("est2", DecisionTreeRegressor(max_depth=3, random_state=1)),
        ],
        cv=3,
    ),
    SVC: dict(max_iter=-1),
    SVR: dict(max_iter=-1),
    TargetEncoder: dict(cv=3),
    TheilSenRegressor: dict(max_iter=5, max_subpopulation=100),
    # TruncatedSVD doesn't run with n_components = n_features
    TruncatedSVD: dict(n_iter=5, n_components=1),
    TSNE: dict(perplexity=2),
    TunedThresholdClassifierCV: dict(estimator=LogisticRegression(C=1), cv=3),
    TweedieRegressor: dict(max_iter=5),
    VotingClassifier: dict(
        estimators=[
            ("est1", DecisionTreeClassifier(max_depth=3, random_state=0)),
            ("est2", DecisionTreeClassifier(max_depth=3, random_state=1)),
        ]
    ),
    VotingRegressor: dict(
        estimators=[
            ("est1", DecisionTreeRegressor(max_depth=3, random_state=0)),
            ("est2", DecisionTreeRegressor(max_depth=3, random_state=1)),
        ]
    ),
}

# This dictionary stores parameters for specific checks. It also enables running the
# same check with multiple instances of the same estimator with different parameters.
# The special key "*" allows to apply the parameters to all checks.
# TODO(devtools): allow third-party developers to pass test specific params to checks
PER_ESTIMATOR_CHECK_PARAMS: dict = {
    # TODO(devtools): check that function names here exist in checks for the estimator
    # TODO(devtools): write a test for the same thing with tags._xfail_checks
    AgglomerativeClustering: {"check_dict_unchanged": dict(n_clusters=1)},
    BayesianGaussianMixture: {"check_dict_unchanged": dict(max_iter=5, n_init=2)},
    BernoulliRBM: {"check_dict_unchanged": dict(n_components=1, n_iter=5)},
    Birch: {"check_dict_unchanged": dict(n_clusters=1)},
    BisectingKMeans: {"check_dict_unchanged": dict(max_iter=5, n_clusters=1, n_init=2)},
    CCA: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    DictionaryLearning: {
        "check_dict_unchanged": dict(
            max_iter=20, n_components=1, transform_algorithm="lasso_lars"
        )
    },
    FactorAnalysis: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    FastICA: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    FeatureAgglomeration: {"check_dict_unchanged": dict(n_clusters=1)},
    GaussianMixture: {"check_dict_unchanged": dict(max_iter=5, n_init=2)},
    GaussianRandomProjection: {"check_dict_unchanged": dict(n_components=1)},
    IncrementalPCA: {"check_dict_unchanged": dict(batch_size=10, n_components=1)},
    Isomap: {"check_dict_unchanged": dict(n_components=1)},
    KMeans: {"check_dict_unchanged": dict(max_iter=5, n_clusters=1, n_init=2)},
    KernelPCA: {"check_dict_unchanged": dict(n_components=1)},
    LatentDirichletAllocation: {
        "check_dict_unchanged": dict(batch_size=10, max_iter=5, n_components=1)
    },
    LinearDiscriminantAnalysis: {"check_dict_unchanged": dict(n_components=1)},
    LocallyLinearEmbedding: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    MDS: {"check_dict_unchanged": dict(max_iter=5, n_components=1, n_init=2)},
    MiniBatchDictionaryLearning: {
        "check_dict_unchanged": dict(batch_size=10, max_iter=5, n_components=1)
    },
    MiniBatchKMeans: {
        "check_dict_unchanged": dict(batch_size=10, max_iter=5, n_clusters=1, n_init=2)
    },
    MiniBatchNMF: {
        "check_dict_unchanged": dict(
            batch_size=10, fresh_restarts=True, max_iter=20, n_components=1
        )
    },
    MiniBatchSparsePCA: {
        "check_dict_unchanged": dict(batch_size=10, max_iter=5, n_components=1)
    },
    NMF: {"check_dict_unchanged": dict(max_iter=500, n_components=1)},
    NeighborhoodComponentsAnalysis: {
        "check_dict_unchanged": dict(max_iter=5, n_components=1)
    },
    Nystroem: {"check_dict_unchanged": dict(n_components=1)},
    PCA: {"check_dict_unchanged": dict(n_components=1)},
    PLSCanonical: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    PLSRegression: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    PLSSVD: {"check_dict_unchanged": dict(n_components=1)},
    PolynomialCountSketch: {"check_dict_unchanged": dict(n_components=1)},
    RBFSampler: {"check_dict_unchanged": dict(n_components=1)},
    SkewedChi2Sampler: {"check_dict_unchanged": dict(n_components=1)},
    SparsePCA: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    SparseRandomProjection: {"check_dict_unchanged": dict(n_components=1)},
    SpectralBiclustering: {
        "check_dict_unchanged": dict(n_best=1, n_clusters=1, n_components=1, n_init=2)
    },
    SpectralClustering: {
        "check_dict_unchanged": dict(n_clusters=1, n_components=1, n_init=2)
    },
    SpectralCoclustering: {"check_dict_unchanged": dict(n_clusters=1, n_init=2)},
    SpectralEmbedding: {"check_dict_unchanged": dict(eigen_tol=1e-05, n_components=1)},
    TSNE: {"check_dict_unchanged": dict(n_components=1, perplexity=2)},
    TruncatedSVD: {"check_dict_unchanged": dict(n_components=1)},
}


def _tested_estimators(type_filter=None):
    for name, Estimator in all_estimators(type_filter=type_filter):
        try:
            for estimator in _construct_instances(Estimator):
                yield estimator
        except SkipTest:
            continue


SKIPPED_ESTIMATORS = [SparseCoder]


def _construct_instances(Estimator):
    """Construct Estimator instances if possible.

    If parameter sets in INIT_PARAMS are provided, use them. If there are a list
    of parameter sets, return one instance for each set.
    """
    if Estimator in SKIPPED_ESTIMATORS:
        msg = f"Can't instantiate estimator {Estimator.__name__}"
        # raise additional warning to be shown by pytest
        warnings.warn(msg, SkipTestWarning)
        raise SkipTest(msg)

    if Estimator in INIT_PARAMS:
        param_sets = INIT_PARAMS[Estimator]
        if not isinstance(param_sets, list):
            param_sets = [param_sets]
        for params in param_sets:
            est = Estimator(**params)
            yield est
    else:
        yield Estimator()


def _get_check_estimator_ids(obj):
    """Create pytest ids for checks.

    When `obj` is an estimator, this returns the pprint version of the
    estimator (with `print_changed_only=True`). When `obj` is a function, the
    name of the function is returned with its keyword arguments.

    `_get_check_estimator_ids` is designed to be used as the `id` in
    `pytest.mark.parametrize` where `check_estimator(..., generate_only=True)`
    is yielding estimators and checks.

    Parameters
    ----------
    obj : estimator or function
        Items generated by `check_estimator`.

    Returns
    -------
    id : str or None

    See Also
    --------
    check_estimator
    """
    if isfunction(obj):
        return obj.__name__
    if isinstance(obj, partial):
        if not obj.keywords:
            return obj.func.__name__
        kwstring = ",".join(["{}={}".format(k, v) for k, v in obj.keywords.items()])
        return "{}({})".format(obj.func.__name__, kwstring)
    if hasattr(obj, "get_params"):
        with config_context(print_changed_only=True):
            return re.sub(r"\s", "", str(obj))


def _yield_instances_for_check(check, estimator_orig):
    """Yield instances for a check.

    For most estimators, this is a no-op.

    For estimators which have an entry in PER_ESTIMATOR_CHECK_PARAMS, this will yield
    an estimator for each parameter set in PER_ESTIMATOR_CHECK_PARAMS[estimator].
    """
    # TODO(devtools): enable this behavior for third party estimators as well
    if type(estimator_orig) not in PER_ESTIMATOR_CHECK_PARAMS:
        yield estimator_orig
        return

    check_params = PER_ESTIMATOR_CHECK_PARAMS[type(estimator_orig)]

    try:
        check_name = check.__name__
    except AttributeError:
        # partial tests
        check_name = check.func.__name__

    if check_name not in check_params:
        yield estimator_orig
        return

    param_set = check_params[check_name]
    if isinstance(param_set, dict):
        param_set = [param_set]

    for params in param_set:
        estimator = clone(estimator_orig)
        estimator.set_params(**params)
        yield estimator
