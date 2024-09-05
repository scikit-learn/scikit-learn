# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import re
import warnings
from functools import partial
from inspect import isfunction, signature
from itertools import product

from sklearn import config_context
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
    LatentDirichletAllocation,
    MiniBatchDictionaryLearning,
    MiniBatchNMF,
    MiniBatchSparsePCA,
    SparseCoder,
    SparsePCA,
    TruncatedSVD,
)
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
from sklearn.manifold import MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding
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
from sklearn.pipeline import FeatureUnion, Pipeline, make_pipeline
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
from sklearn.utils._testing import SkipTest, set_random_state

CROSS_DECOMPOSITION = ["PLSCanonical", "PLSRegression", "CCA", "PLSSVD"]

# The following dictionary is to indicate constructor arguments suitable for the test
# suite, which uses very small datasets, and is intended to run rather quickly.
INIT_PARAMS = {
    AdaBoostClassifier: {"n_estimators": 5},
    AdaBoostRegressor: {"n_estimators": 5},
    AffinityPropagation: {"max_iter": 5},
    AgglomerativeClustering: {"n_clusters": 2},
    ARDRegression: {"max_iter": 5},
    BaggingClassifier: {"n_estimators": 5},
    BaggingRegressor: {"n_estimators": 5},
    BayesianGaussianMixture: {"n_init": 2, "max_iter": 5},
    BayesianRidge: {"max_iter": 5},
    BernoulliRBM: {"n_iter": 5, "batch_size": 10},
    Birch: {"n_clusters": 2},
    BisectingKMeans: {"n_init": 2, "n_clusters": 2, "max_iter": 5},
    CalibratedClassifierCV: {"estimator": LogisticRegression(C=1), "cv": 3},
    CCA: {"n_components": 1, "max_iter": 5},
    ClassifierChain: {"base_estimator": LogisticRegression(C=1), "cv": 3},
    ColumnTransformer: {"transformers": [("trans1", StandardScaler(), [0, 1])]},
    DictionaryLearning: {"max_iter": 20, "transform_algorithm": "lasso_lars"},
    # the default strategy prior would output constant predictions and fail
    # for check_classifiers_predictions
    DummyClassifier: {"strategy": "stratified"},
    ElasticNetCV: {"max_iter": 5, "cv": 3},
    ElasticNet: {"max_iter": 5},
    ExtraTreesClassifier: {"n_estimators": 5},
    ExtraTreesRegressor: {"n_estimators": 5},
    FactorAnalysis: {"max_iter": 5},
    FastICA: {"max_iter": 5},
    FeatureAgglomeration: {"n_clusters": 2},
    FeatureUnion: {"transformer_list": [("trans1", StandardScaler())]},
    FixedThresholdClassifier: {"estimator": LogisticRegression(C=1)},
    GammaRegressor: {"max_iter": 5},
    GaussianMixture: {"n_init": 2, "max_iter": 5},
    # Due to the jl lemma and often very few samples, the number
    # of components of the random matrix projection will be probably
    # greater than the number of features.
    # So we impose a smaller number (avoid "auto" mode)
    GaussianRandomProjection: {"n_components": 2},
    GradientBoostingClassifier: {"n_estimators": 5},
    GradientBoostingRegressor: {"n_estimators": 5},
    GraphicalLassoCV: {"max_iter": 5, "cv": 3},
    GraphicalLasso: {"max_iter": 5},
    GridSearchCV: {
        "estimator": LogisticRegression(C=1),
        "param_grid": {"C": [1.0]},
        "cv": 3,
    },
    HalvingGridSearchCV: {
        "estimator": Ridge(),
        "min_resources": "smallest",
        "param_grid": {"alpha": [0.1, 1.0]},
        "random_state": 0,
        "cv": 2,
        "error_score": "raise",
    },
    HalvingRandomSearchCV: {
        "estimator": Ridge(),
        "param_distributions": {"alpha": [0.1, 1.0]},
        "min_resources": "smallest",
        "cv": 2,
        "error_score": "raise",
        "random_state": 0,
    },
    HDBSCAN: {"min_samples": 1},
    # The default min_samples_leaf (20) isn't appropriate for small
    # datasets (only very shallow trees are built) that the checks use.
    HistGradientBoostingClassifier: {"max_iter": 5, "min_samples_leaf": 5},
    HistGradientBoostingRegressor: {"max_iter": 5, "min_samples_leaf": 5},
    HuberRegressor: {"max_iter": 5},
    IncrementalPCA: {"batch_size": 10},
    IsolationForest: {"n_estimators": 5},
    KMeans: {"n_init": 2, "n_clusters": 2, "max_iter": 5},
    LabelPropagation: {"max_iter": 5},
    LabelSpreading: {"max_iter": 5},
    LarsCV: {"max_iter": 5, "cv": 3},
    LassoCV: {"max_iter": 5, "cv": 3},
    LassoLarsCV: {"max_iter": 5, "cv": 3},
    # Noise variance estimation does not work when `n_samples < n_features`.
    # We need to provide the noise variance explicitly.
    LassoLarsIC: {"max_iter": 5, "noise_variance": 1.0},
    LassoLars: {"max_iter": 5},
    Lasso: {"max_iter": 5},
    LatentDirichletAllocation: {"max_iter": 5, "batch_size": 10},
    LinearSVC: {"max_iter": 20},
    LinearSVR: {"max_iter": 20},
    LocallyLinearEmbedding: {"max_iter": 5},
    LogisticRegressionCV: {"max_iter": 5, "cv": 3},
    LogisticRegression: {"max_iter": 5},
    MDS: {"n_init": 2, "max_iter": 5},
    # In the case of check_fit2d_1sample, bandwidth is set to None and
    # is thus estimated. De facto it is 0.0 as a single sample is provided
    # and this makes the test fails. Hence we give it a placeholder value.
    MeanShift: {"max_iter": 5, "bandwidth": 1.0},
    MiniBatchDictionaryLearning: {"batch_size": 10, "max_iter": 5},
    MiniBatchKMeans: {"n_init": 2, "n_clusters": 2, "max_iter": 5, "batch_size": 10},
    MiniBatchNMF: {"batch_size": 10, "max_iter": 20, "fresh_restarts": True},
    MiniBatchSparsePCA: {"max_iter": 5, "batch_size": 10},
    MLPClassifier: {"max_iter": 100},
    MLPRegressor: {"max_iter": 100},
    MultiOutputClassifier: {"estimator": LogisticRegression(C=1)},
    MultiOutputRegressor: {"estimator": Ridge()},
    MultiTaskElasticNetCV: {"max_iter": 5, "cv": 3},
    MultiTaskElasticNet: {"max_iter": 5},
    MultiTaskLassoCV: {"max_iter": 5, "cv": 3},
    MultiTaskLasso: {"max_iter": 5},
    NeighborhoodComponentsAnalysis: {"max_iter": 5},
    NMF: {"max_iter": 500},
    NuSVC: {"max_iter": -1},
    NuSVR: {"max_iter": -1},
    OneClassSVM: {"max_iter": -1},
    OneHotEncoder: {"handle_unknown": "ignore"},
    OneVsOneClassifier: {"estimator": LogisticRegression(C=1)},
    OneVsRestClassifier: {"estimator": LogisticRegression(C=1)},
    OrthogonalMatchingPursuitCV: {"cv": 3},
    OutputCodeClassifier: {"estimator": LogisticRegression(C=1)},
    PassiveAggressiveClassifier: {"max_iter": 5},
    PassiveAggressiveRegressor: {"max_iter": 5},
    Perceptron: {"max_iter": 5},
    Pipeline: {"steps": [("scaler", StandardScaler()), ("est", Ridge())]},
    PLSCanonical: {"n_components": 1, "max_iter": 5},
    PLSRegression: {"n_components": 1, "max_iter": 5},
    PLSSVD: {"n_components": 1},
    PoissonRegressor: {"max_iter": 5},
    RandomForestClassifier: {"n_estimators": 5},
    RandomForestRegressor: {"n_estimators": 5},
    RandomizedSearchCV: {
        "estimator": LogisticRegression(C=1),
        "param_distributions": {"C": [1.0]},
        "n_iter": 5,
        "cv": 3,
    },
    RandomTreesEmbedding: {"n_estimators": 5},
    # `RANSACRegressor` will raise an error with any model other
    # than `LinearRegression` if we don't fix `min_samples` parameter.
    # For common test, we can enforce using `LinearRegression` that
    # is the default estimator in `RANSACRegressor` instead of `Ridge`.
    RANSACRegressor: {"estimator": LinearRegression(), "max_trials": 10},
    RegressorChain: {"base_estimator": Ridge(), "cv": 3},
    RFECV: {"estimator": LogisticRegression(C=1), "cv": 3},
    RFE: {"estimator": LogisticRegression(C=1)},
    # be tolerant of noisy datasets (not actually speed)
    SelectFdr: {"alpha": 0.5},
    # Increases coverage because SGDRegressor has partial_fit
    SelectFromModel: {"estimator": SGDRegressor(random_state=0)},
    # SelectKBest has a default of k=10
    # which is more feature than we have in most case.
    SelectKBest: {"k": 1},
    SelfTrainingClassifier: {"estimator": LogisticRegression(C=1), "max_iter": 5},
    SequentialFeatureSelector: {"estimator": LogisticRegression(C=1), "cv": 3},
    SGDClassifier: {"max_iter": 5},
    SGDOneClassSVM: {"max_iter": 5},
    SGDRegressor: {"max_iter": 5},
    SparsePCA: {"max_iter": 5},
    # Due to the jl lemma and often very few samples, the number
    # of components of the random matrix projection will be probably
    # greater than the number of features.
    # So we impose a smaller number (avoid "auto" mode)
    SparseRandomProjection: {"n_components": 2},
    SpectralBiclustering: {"n_init": 2, "n_best": 1, "n_clusters": 2},
    SpectralClustering: {"n_init": 2, "n_clusters": 2},
    SpectralCoclustering: {"n_init": 2, "n_clusters": 2},
    # Default "auto" parameter can lead to different ordering of eigenvalues on
    # windows: #24105
    SpectralEmbedding: {"eigen_tol": 1e-05},
    StackingClassifier: {
        "estimators": [
            ("est1", DecisionTreeClassifier(max_depth=3, random_state=0)),
            ("est2", DecisionTreeClassifier(max_depth=3, random_state=1)),
        ],
        "cv": 3,
    },
    StackingRegressor: {
        "estimators": [
            ("est1", DecisionTreeRegressor(max_depth=3, random_state=0)),
            ("est2", DecisionTreeRegressor(max_depth=3, random_state=1)),
        ],
        "cv": 3,
    },
    SVC: {"max_iter": -1},
    SVR: {"max_iter": -1},
    TargetEncoder: {"cv": 3},
    TheilSenRegressor: {"max_iter": 5, "max_subpopulation": 100},
    # TruncatedSVD doesn't run with n_components = n_features
    TruncatedSVD: {"n_iter": 5, "n_components": 1},
    TSNE: {"perplexity": 2},
    TunedThresholdClassifierCV: {"estimator": LogisticRegression(C=1), "cv": 3},
    TweedieRegressor: {"max_iter": 5},
    VotingClassifier: {
        "estimators": [
            ("est1", DecisionTreeClassifier(max_depth=3, random_state=0)),
            ("est2", DecisionTreeClassifier(max_depth=3, random_state=1)),
        ]
    },
    VotingRegressor: {
        "estimators": [
            ("est1", DecisionTreeRegressor(max_depth=3, random_state=0)),
            ("est2", DecisionTreeRegressor(max_depth=3, random_state=1)),
        ]
    },
}


def _tested_estimators(type_filter=None):
    for name, Estimator in all_estimators(type_filter=type_filter):
        try:
            estimator = _construct_instance(Estimator)
        except SkipTest:
            continue

        yield estimator


def _generate_pipeline():
    """Generator of simple pipeline to check compliance of the
    :class:`~sklearn.pipeline.Pipeline` class.
    """
    for final_estimator in [Ridge(), LogisticRegression()]:
        yield Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("final_estimator", final_estimator),
            ]
        )


SKIPPED_ESTIMATORS = [SparseCoder]


def _construct_instance(Estimator):
    """Construct Estimator instance if possible."""
    if Estimator in SKIPPED_ESTIMATORS:
        msg = f"Can't instantiate estimator {Estimator.__name__}"
        # raise additional warning to be shown by pytest
        warnings.warn(msg, SkipTestWarning)
        raise SkipTest(msg)

    if Estimator in INIT_PARAMS:
        estimator = Estimator(**INIT_PARAMS[Estimator])
    else:
        estimator = Estimator()
    return estimator


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


def _generate_search_cv_instances():
    """Generator of `SearchCV` instances to check their compliance with scikit-learn."""
    for SearchCV, (Estimator, param_grid) in product(
        [
            GridSearchCV,
            HalvingGridSearchCV,
            RandomizedSearchCV,
            HalvingGridSearchCV,
        ],
        [
            (Ridge, {"alpha": [0.1, 1.0]}),
            (LogisticRegression, {"C": [0.1, 1.0]}),
        ],
    ):
        init_params = signature(SearchCV).parameters
        extra_params = (
            {"min_resources": "smallest"} if "min_resources" in init_params else {}
        )
        search_cv = SearchCV(
            Estimator(), param_grid, cv=2, error_score="raise", **extra_params
        )
        set_random_state(search_cv)
        yield search_cv

    for SearchCV, (Estimator, param_grid) in product(
        [
            GridSearchCV,
            HalvingGridSearchCV,
            RandomizedSearchCV,
            HalvingRandomSearchCV,
        ],
        [
            (Ridge, {"ridge__alpha": [0.1, 1.0]}),
            (LogisticRegression, {"logisticregression__C": [0.1, 1.0]}),
        ],
    ):
        init_params = signature(SearchCV).parameters
        extra_params = (
            {"min_resources": "smallest"} if "min_resources" in init_params else {}
        )
        search_cv = SearchCV(
            make_pipeline(PCA(), Estimator()), param_grid, cv=2, **extra_params
        ).set_params(error_score="raise")
        set_random_state(search_cv)
        yield search_cv
