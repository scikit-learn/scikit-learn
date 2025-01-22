# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import re
import warnings
from contextlib import suppress
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
from sklearn.frozen import FrozenEstimator
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
    QuantileRegressor,
    RANSACRegressor,
    Ridge,
    RidgeClassifier,
    RidgeCV,
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
from sklearn.naive_bayes import CategoricalNB
from sklearn.neighbors import (
    KernelDensity,
    KNeighborsClassifier,
    KNeighborsRegressor,
    KNeighborsTransformer,
    NeighborhoodComponentsAnalysis,
    RadiusNeighborsTransformer,
)
from sklearn.neural_network import BernoulliRBM, MLPClassifier, MLPRegressor
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.preprocessing import (
    KBinsDiscretizer,
    OneHotEncoder,
    SplineTransformer,
    StandardScaler,
    TargetEncoder,
)
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
from sklearn.utils._tags import get_tags
from sklearn.utils._testing import SkipTest
from sklearn.utils.fixes import parse_version, sp_base_version

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
    ClassifierChain: dict(estimator=LogisticRegression(C=1), cv=3),
    ColumnTransformer: dict(transformers=[("trans1", StandardScaler(), [0, 1])]),
    DictionaryLearning: dict(max_iter=20, transform_algorithm="lasso_lars"),
    # the default strategy prior would output constant predictions and fail
    # for check_classifiers_predictions
    DummyClassifier: [dict(strategy="stratified"), dict(strategy="most_frequent")],
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
    KNeighborsClassifier: [dict(n_neighbors=2), dict(metric="precomputed")],
    KNeighborsRegressor: [dict(n_neighbors=2), dict(metric="precomputed")],
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
    RegressorChain: dict(estimator=Ridge(), cv=3),
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
    SVC: [dict(max_iter=-1), dict(kernel="precomputed")],
    SVR: [dict(max_iter=-1), dict(kernel="precomputed")],
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
    AgglomerativeClustering: {"check_dict_unchanged": dict(n_clusters=1)},
    BayesianGaussianMixture: {"check_dict_unchanged": dict(max_iter=5, n_init=2)},
    BernoulliRBM: {"check_dict_unchanged": dict(n_components=1, n_iter=5)},
    Birch: {"check_dict_unchanged": dict(n_clusters=1)},
    BisectingKMeans: {"check_dict_unchanged": dict(max_iter=5, n_clusters=1, n_init=2)},
    CCA: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    DecisionTreeRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(criterion="squared_error"),
            dict(criterion="absolute_error"),
            dict(criterion="friedman_mse"),
            dict(criterion="poisson"),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(criterion="squared_error"),
            dict(criterion="absolute_error"),
            dict(criterion="friedman_mse"),
            dict(criterion="poisson"),
        ],
    },
    DecisionTreeClassifier: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(criterion="gini"),
            dict(criterion="log_loss"),
            dict(criterion="entropy"),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(criterion="gini"),
            dict(criterion="log_loss"),
            dict(criterion="entropy"),
        ],
    },
    DictionaryLearning: {
        "check_dict_unchanged": dict(
            max_iter=20, n_components=1, transform_algorithm="lasso_lars"
        )
    },
    FactorAnalysis: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    FastICA: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    FeatureAgglomeration: {"check_dict_unchanged": dict(n_clusters=1)},
    FeatureUnion: {
        "check_estimator_sparse_tag": [
            dict(transformer_list=[("trans1", StandardScaler())]),
            dict(
                transformer_list=[
                    ("trans1", StandardScaler(with_mean=False)),
                    ("trans2", "drop"),
                    ("trans3", "passthrough"),
                ]
            ),
        ]
    },
    GammaRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="newton-cholesky"),
            dict(solver="lbfgs"),
        ],
    },
    GaussianMixture: {"check_dict_unchanged": dict(max_iter=5, n_init=2)},
    GaussianRandomProjection: {"check_dict_unchanged": dict(n_components=1)},
    IncrementalPCA: {"check_dict_unchanged": dict(batch_size=10, n_components=1)},
    Isomap: {"check_dict_unchanged": dict(n_components=1)},
    KMeans: {"check_dict_unchanged": dict(max_iter=5, n_clusters=1, n_init=2)},
    KernelPCA: {"check_dict_unchanged": dict(n_components=1)},
    LassoLars: {"check_non_transformer_estimators_n_iter": dict(alpha=0.0)},
    LatentDirichletAllocation: {
        "check_dict_unchanged": dict(batch_size=10, max_iter=5, n_components=1)
    },
    LinearDiscriminantAnalysis: {"check_dict_unchanged": dict(n_components=1)},
    LinearRegression: {
        "check_estimator_sparse_tag": [dict(positive=False), dict(positive=True)],
        "check_sample_weight_equivalence_on_dense_data": [
            dict(positive=False),
            dict(positive=True),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [dict(tol=1e-12)],
    },
    LocallyLinearEmbedding: {"check_dict_unchanged": dict(max_iter=5, n_components=1)},
    LogisticRegression: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="lbfgs"),
            dict(solver="liblinear"),
            dict(solver="newton-cg"),
            dict(solver="newton-cholesky"),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(solver="liblinear"),
        ],
    },
    MDS: {"check_dict_unchanged": dict(max_iter=5, n_components=1, n_init=2)},
    MLPClassifier: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="lbfgs"),
        ]
    },
    MLPRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="sgd", tol=1e-2, random_state=42),
        ]
    },
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
    PoissonRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="newton-cholesky"),
            dict(solver="lbfgs"),
        ],
    },
    PolynomialCountSketch: {"check_dict_unchanged": dict(n_components=1)},
    QuantileRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(quantile=0.5),
            dict(quantile=0.75),
            dict(solver="highs-ds"),
            dict(solver="highs-ipm"),
        ],
    },
    RBFSampler: {"check_dict_unchanged": dict(n_components=1)},
    Ridge: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="svd"),
            dict(solver="cholesky"),
            dict(solver="sparse_cg"),
            dict(solver="lsqr"),
            dict(solver="lbfgs", positive=True),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(solver="sparse_cg"),
            dict(solver="lsqr"),
        ],
    },
    RidgeClassifier: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="svd"),
            dict(solver="cholesky"),
            dict(solver="sparse_cg"),
            dict(solver="lsqr"),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(solver="sparse_cg"),
            dict(solver="lsqr"),
        ],
    },
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
    StandardScaler: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(with_mean=True),
            dict(with_mean=False),
        ],
        "check_sample_weight_equivalence_on_sparse_data": [
            dict(with_mean=False),
        ],
    },
    TSNE: {"check_dict_unchanged": dict(n_components=1, perplexity=2)},
    TruncatedSVD: {"check_dict_unchanged": dict(n_components=1)},
    TweedieRegressor: {
        "check_sample_weight_equivalence_on_dense_data": [
            dict(solver="newton-cholesky"),
            dict(solver="lbfgs"),
        ],
    },
}


def _tested_estimators(type_filter=None):
    for _, Estimator in all_estimators(type_filter=type_filter):
        with suppress(SkipTest):
            for estimator in _construct_instances(Estimator):
                yield estimator


SKIPPED_ESTIMATORS = [SparseCoder, FrozenEstimator]


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


PER_ESTIMATOR_XFAIL_CHECKS = {
    AdaBoostClassifier: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    AdaBoostRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    BaggingClassifier: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    BaggingRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    BernoulliRBM: {
        "check_methods_subset_invariance": ("fails for the decision_function method"),
        "check_methods_sample_order_invariance": ("fails for the score_samples method"),
    },
    BisectingKMeans: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    CategoricalNB: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    ColumnTransformer: {
        "check_estimators_empty_data_messages": "FIXME",
        "check_estimators_nan_inf": "FIXME",
        "check_estimator_sparse_array": "FIXME",
        "check_estimator_sparse_matrix": "FIXME",
        "check_fit1d": "FIXME",
        "check_fit2d_predict1d": "FIXME",
        "check_complex_data": "FIXME",
        "check_fit2d_1feature": "FIXME",
    },
    DummyClassifier: {
        "check_methods_subset_invariance": "fails for the predict method",
        "check_methods_sample_order_invariance": "fails for the predict method",
    },
    FeatureUnion: {
        "check_estimators_overwrite_params": "FIXME",
        "check_estimators_nan_inf": "FIXME",
        "check_dont_overwrite_parameters": "FIXME",
    },
    FixedThresholdClassifier: {
        "check_classifiers_train": "Threshold at probability 0.5 does not hold",
        "check_sample_weight_equivalence_on_dense_data": (
            "Due to the cross-validation and sample ordering, removing a sample"
            " is not strictly equal to putting is weight to zero. Specific unit"
            " tests are added for TunedThresholdClassifierCV specifically."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    GradientBoostingClassifier: {
        # TODO: investigate failure see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    GradientBoostingRegressor: {
        # TODO: investigate failure see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    GridSearchCV: {
        "check_supervised_y_2d": "DataConversionWarning not caught",
        "check_requires_y_none": "Doesn't fail gracefully",
    },
    HalvingGridSearchCV: {
        "check_fit2d_1sample": (
            "Fail during parameter check since min/max resources requires"
            " more samples"
        ),
        "check_estimators_nan_inf": "FIXME",
        "check_classifiers_one_label_sample_weights": "FIXME",
        "check_fit2d_1feature": "FIXME",
        "check_supervised_y_2d": "DataConversionWarning not caught",
        "check_requires_y_none": "Doesn't fail gracefully",
    },
    HalvingRandomSearchCV: {
        "check_fit2d_1sample": (
            "Fail during parameter check since min/max resources requires"
            " more samples"
        ),
        "check_estimators_nan_inf": "FIXME",
        "check_classifiers_one_label_sample_weights": "FIXME",
        "check_fit2d_1feature": "FIXME",
        "check_supervised_y_2d": "DataConversionWarning not caught",
        "check_requires_y_none": "Doesn't fail gracefully",
    },
    HistGradientBoostingClassifier: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    HistGradientBoostingRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    IsolationForest: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    KBinsDiscretizer: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    KernelDensity: {
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight must have positive values"
        ),
    },
    KMeans: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    KNeighborsTransformer: {
        "check_methods_sample_order_invariance": "check is not applicable."
    },
    LinearSVC: {
        # TODO: replace by a statistical test when _dual=True, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_non_transformer_estimators_n_iter": (
            "n_iter_ cannot be easily accessed."
        ),
    },
    LinearSVR: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    LogisticRegression: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    MiniBatchKMeans: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    NuSVC: {
        "check_class_weight_classifiers": "class_weight is ignored.",
        # TODO: fix sample_weight handling of this estimator when probability=False
        # TODO: replace by a statistical test when probability=True
        # see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_classifiers_one_label_sample_weights": (
            "specified nu is infeasible for the fit."
        ),
    },
    NuSVR: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    Nystroem: {
        "check_transformer_preserves_dtypes": (
            "dtypes are preserved but not at a close enough precision"
        )
    },
    OneClassSVM: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    Perceptron: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    Pipeline: {
        "check_dont_overwrite_parameters": (
            "Pipeline changes the `steps` parameter, which it shouldn't."
            "Therefore this test is x-fail until we fix this."
        ),
        "check_estimators_overwrite_params": (
            "Pipeline changes the `steps` parameter, which it shouldn't."
            "Therefore this test is x-fail until we fix this."
        ),
    },
    RadiusNeighborsTransformer: {
        "check_methods_sample_order_invariance": "check is not applicable."
    },
    RandomForestClassifier: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    RandomForestRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    RandomizedSearchCV: {
        "check_supervised_y_2d": "DataConversionWarning not caught",
        "check_requires_y_none": "Doesn't fail gracefully",
    },
    RandomTreesEmbedding: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    RANSACRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    Ridge: {
        "check_non_transformer_estimators_n_iter": (
            "n_iter_ cannot be easily accessed."
        )
    },
    RidgeClassifier: {
        "check_non_transformer_estimators_n_iter": (
            "n_iter_ cannot be easily accessed."
        )
    },
    RidgeCV: {
        "check_sample_weight_equivalence_on_dense_data": (
            "GridSearchCV does not forward the weights to the scorer by default."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    SelfTrainingClassifier: {
        "check_non_transformer_estimators_n_iter": "n_iter_ can be 0."
    },
    SGDClassifier: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    SGDOneClassSVM: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    SGDRegressor: {
        # TODO: replace by a statistical test, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    SpectralCoclustering: {
        "check_estimators_dtypes": "raises nan error",
        "check_fit2d_1sample": "_scale_normalize fails",
        "check_fit2d_1feature": "raises apply_along_axis error",
        "check_estimator_sparse_matrix": "does not fail gracefully",
        "check_estimator_sparse_array": "does not fail gracefully",
        "check_methods_subset_invariance": "empty array passed inside",
        "check_dont_overwrite_parameters": "empty array passed inside",
        "check_fit2d_predict1d": "empty array passed inside",
        # ValueError: Found array with 0 feature(s) (shape=(23, 0))
        # while a minimum of 1 is required.
        "check_dict_unchanged": "FIXME",
    },
    SpectralBiclustering: {
        "check_estimators_dtypes": "raises nan error",
        "check_fit2d_1sample": "_scale_normalize fails",
        "check_fit2d_1feature": "raises apply_along_axis error",
        "check_estimator_sparse_matrix": "does not fail gracefully",
        "check_estimator_sparse_array": "does not fail gracefully",
        "check_methods_subset_invariance": "empty array passed inside",
        "check_dont_overwrite_parameters": "empty array passed inside",
        "check_fit2d_predict1d": "empty array passed inside",
    },
    SVC: {
        # TODO: fix sample_weight handling of this estimator when probability=False
        # TODO: replace by a statistical test when probability=True
        # see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    SVR: {
        # TODO: fix sample_weight handling of this estimator, see meta-issue #16298
        "check_sample_weight_equivalence_on_dense_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
        "check_sample_weight_equivalence_on_sparse_data": (
            "sample_weight is not equivalent to removing/repeating samples."
        ),
    },
    TunedThresholdClassifierCV: {
        "check_classifiers_train": "Threshold at probability 0.5 does not hold",
        "check_sample_weight_equivalence_on_dense_data": (
            "Due to the cross-validation and sample ordering, removing a sample"
            " is not strictly equal to putting is weight to zero. Specific unit"
            " tests are added for TunedThresholdClassifierCV specifically."
        ),
    },
}

# TODO: remove when scipy min version >= 1.11
if sp_base_version < parse_version("1.11"):
    PER_ESTIMATOR_XFAIL_CHECKS[SplineTransformer] = {
        "check_estimators_pickle": (
            "scipy < 1.11 implementation of _bsplines does not"
            "support const memory views."
        ),
    }


def _get_expected_failed_checks(estimator):
    """Get the expected failed checks for all estimators in scikit-learn."""
    failed_checks = PER_ESTIMATOR_XFAIL_CHECKS.get(type(estimator), {})

    tags = get_tags(estimator)

    # all xfail marks that depend on the instance, come here. As of now, we have only
    # these two cases.
    if type(estimator) in [KNeighborsClassifier, KNeighborsRegressor]:
        if tags.input_tags.pairwise:
            failed_checks.update(
                {
                    "check_n_features_in_after_fitting": "FIXME",
                    "check_dataframe_column_names_consistency": "FIXME",
                }
            )

    return failed_checks
