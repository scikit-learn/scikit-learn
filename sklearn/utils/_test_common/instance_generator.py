# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import re
import warnings
from functools import partial
from inspect import isfunction, signature
from itertools import product

from sklearn import config_context
from sklearn.base import RegressorMixin
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
)
from sklearn.exceptions import SkipTestWarning
from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.feature_selection import (
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
    GridSearchCV,
    HalvingGridSearchCV,
    HalvingRandomSearchCV,
    RandomizedSearchCV,
    TunedThresholdClassifierCV,
)
from sklearn.multioutput import ClassifierChain, RegressorChain
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.neural_network import BernoulliRBM, MLPClassifier, MLPRegressor
from sklearn.pipeline import Pipeline, make_pipeline
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
TEST_PARAMS = {
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
    CalibratedClassifierCV: dict(cv=3),
    CCA: dict(n_components=1, max_iter=5),
    ClassifierChain: dict(cv=3),
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
    GridSearchCV: dict(cv=3),
    HalvingGridSearchCV: dict(cv=3),
    HalvingRandomSearchCV: dict(cv=3),
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
    LinearSVR: dict(max_iter=20),
    LinearSVC: dict(max_iter=20),
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
    OrthogonalMatchingPursuitCV: dict(cv=3),
    PassiveAggressiveClassifier: dict(max_iter=5),
    PassiveAggressiveRegressor: dict(max_iter=5),
    Perceptron: dict(max_iter=5),
    PLSCanonical: dict(n_components=1, max_iter=5),
    PLSRegression: dict(n_components=1, max_iter=5),
    PLSSVD: dict(n_components=1),
    PoissonRegressor: dict(max_iter=5),
    RandomForestClassifier: dict(n_estimators=5),
    RandomForestRegressor: dict(n_estimators=5),
    RandomizedSearchCV: dict(n_iter=5, cv=3),
    RandomTreesEmbedding: dict(n_estimators=5),
    RANSACRegressor: dict(max_trials=10),
    RegressorChain: dict(cv=3),
    RFECV: dict(cv=3),
    # be tolerant of noisy datasets (not actually speed)
    SelectFdr: dict(alpha=0.5),
    # SelectKBest has a default of k=10
    # which is more feature than we have in most case.
    SelectKBest: dict(k=1),
    SelfTrainingClassifier: dict(max_iter=5),
    SequentialFeatureSelector: dict(cv=3),
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
    SpectralEmbedding: dict(eigen_tol=1e-5),
    StackingClassifier: dict(cv=3),
    StackingRegressor: dict(cv=3),
    SVC: dict(max_iter=-1),
    SVR: dict(max_iter=-1),
    TargetEncoder: dict(cv=3),
    TheilSenRegressor: dict(max_iter=5, max_subpopulation=100),
    # TruncatedSVD doesn't run with n_components = n_features
    TruncatedSVD: dict(n_iter=5, n_components=1),
    TSNE: dict(perplexity=2),
    TunedThresholdClassifierCV: dict(cv=3),
    TweedieRegressor: dict(max_iter=5),
}


def _set_checking_parameters(estimator):
    """Set the parameters of an estimator instance to speed-up tests and avoid
    deprecation warnings in common test."""
    if type(estimator) in TEST_PARAMS:
        test_params = TEST_PARAMS[type(estimator)]
        estimator.set_params(**test_params)


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


def _construct_instance(Estimator):
    """Construct Estimator instance if possible."""
    required_parameters = getattr(Estimator, "_required_parameters", [])
    if len(required_parameters):
        if required_parameters in (["estimator"], ["base_estimator"]):
            # `RANSACRegressor` will raise an error with any model other
            # than `LinearRegression` if we don't fix `min_samples` parameter.
            # For common test, we can enforce using `LinearRegression` that
            # is the default estimator in `RANSACRegressor` instead of `Ridge`.
            if issubclass(Estimator, RANSACRegressor):
                estimator = Estimator(LinearRegression())
            elif issubclass(Estimator, RegressorMixin):
                estimator = Estimator(Ridge())
            elif issubclass(Estimator, SelectFromModel):
                # Increases coverage because SGDRegressor has partial_fit
                estimator = Estimator(SGDRegressor(random_state=0))
            else:
                estimator = Estimator(LogisticRegression(C=1))
        elif required_parameters in (["estimators"],):
            # Heterogeneous ensemble classes (i.e. stacking, voting)
            if issubclass(Estimator, RegressorMixin):
                estimator = Estimator(
                    estimators=[
                        ("est1", DecisionTreeRegressor(max_depth=3, random_state=0)),
                        ("est2", DecisionTreeRegressor(max_depth=3, random_state=1)),
                    ]
                )
            else:
                estimator = Estimator(
                    estimators=[
                        ("est1", DecisionTreeClassifier(max_depth=3, random_state=0)),
                        ("est2", DecisionTreeClassifier(max_depth=3, random_state=1)),
                    ]
                )
        else:
            msg = (
                f"Can't instantiate estimator {Estimator.__name__} "
                f"parameters {required_parameters}"
            )
            # raise additional warning to be shown by pytest
            warnings.warn(msg, SkipTestWarning)
            raise SkipTest(msg)
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


def _generate_column_transformer_instances():
    """Generate a `ColumnTransformer` instance to check its compliance with
    scikit-learn."""
    yield ColumnTransformer(
        transformers=[
            ("trans1", StandardScaler(), [0, 1]),
        ]
    )


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
