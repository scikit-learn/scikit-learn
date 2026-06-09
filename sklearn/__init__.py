"""Configure global settings and get information about the working environment."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# Machine learning module for Python
# ==================================
#
# sklearn is a Python module integrating classical machine
# learning algorithms in the tightly-knit world of scientific Python
# packages (numpy, scipy, matplotlib).
#
# It aims to provide simple and efficient solutions to learning problems
# that are accessible to everybody and reusable in various contexts:
# machine-learning as a versatile tool for science and engineering.
#
# See https://scikit-learn.org for complete documentation.

import importlib as _importlib
import logging
import os
import random

from sklearn._config import config_context, get_config, set_config

logger = logging.getLogger(__name__)


# PEP0440 compatible formatted version, see:
# https://www.python.org/dev/peps/pep-0440/
#
# Generic release markers:
#   X.Y.0   # For first release after an increment in Y
#   X.Y.Z   # For bugfix releases
#
# Admissible pre-release markers:
#   X.Y.ZaN   # Alpha release
#   X.Y.ZbN   # Beta release
#   X.Y.ZrcN  # Release Candidate
#   X.Y.Z     # Final release
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'
#
__version__ = "1.10.dev0"


# On OSX, we can get a runtime error due to multiple OpenMP libraries loaded
# simultaneously. This can happen for instance when calling BLAS inside a
# prange. Setting the following environment variable allows multiple OpenMP
# libraries to be loaded. It should not degrade performances since we manually
# take care of potential over-subscription performance issues, in sections of
# the code where nested OpenMP loops can happen, by dynamically reconfiguring
# the inner OpenMP runtime to temporarily disable it while under the scope of
# the outer OpenMP parallel section.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "True")

# Workaround issue discovered in intel-openmp 2019.5:
# https://github.com/ContinuumIO/anaconda-issues/issues/11294
os.environ.setdefault("KMP_INIT_AT_FORK", "FALSE")

# `_distributor_init` allows distributors to run custom init code.
# For instance, for the Windows wheel, this is used to pre-load the
# vcomp shared library runtime for OpenMP embedded in the sklearn/.libs
# sub-folder.
# It is necessary to do this prior to importing show_versions as the
# later is linked to the OpenMP runtime to make it possible to introspect
# it and importing it first would fail if the OpenMP dll cannot be found.
from sklearn import __check_build, _distributor_init  # noqa: E402 F401
from sklearn.base import clone  # noqa: E402
from sklearn.utils._show_versions import show_versions  # noqa: E402

_submodules = [
    "calibration",
    "callback",
    "cluster",
    "covariance",
    "cross_decomposition",
    "datasets",
    "decomposition",
    "dummy",
    "ensemble",
    "exceptions",
    "experimental",
    "externals",
    "feature_extraction",
    "feature_selection",
    "frozen",
    "gaussian_process",
    "inspection",
    "isotonic",
    "kernel_approximation",
    "kernel_ridge",
    "linear_model",
    "manifold",
    "metrics",
    "mixture",
    "model_selection",
    "multiclass",
    "multioutput",
    "naive_bayes",
    "neighbors",
    "neural_network",
    "pipeline",
    "preprocessing",
    "random_projection",
    "semi_supervised",
    "svm",
    "tree",
    "discriminant_analysis",
    "impute",
    "compose",
]

# Mapping from a public submodule to the estimators it exposes. This is used to
# make all estimators lazily accessible from the top-level `sklearn` namespace
# through the module-level `__getattr__` below, without eagerly importing every
# submodule when `sklearn` is first imported.
_estimators = {
    "sklearn.calibration": [
        "CalibratedClassifierCV",
    ],
    "sklearn.cluster": [
        "AffinityPropagation",
        "AgglomerativeClustering",
        "Birch",
        "BisectingKMeans",
        "DBSCAN",
        "FeatureAgglomeration",
        "HDBSCAN",
        "KMeans",
        "MeanShift",
        "MiniBatchKMeans",
        "OPTICS",
        "SpectralBiclustering",
        "SpectralClustering",
        "SpectralCoclustering",
    ],
    "sklearn.compose": [
        "ColumnTransformer",
        "TransformedTargetRegressor",
    ],
    "sklearn.covariance": [
        "EllipticEnvelope",
        "EmpiricalCovariance",
        "GraphicalLasso",
        "GraphicalLassoCV",
        "LedoitWolf",
        "MinCovDet",
        "OAS",
        "ShrunkCovariance",
    ],
    "sklearn.cross_decomposition": [
        "CCA",
        "PLSCanonical",
        "PLSRegression",
        "PLSSVD",
    ],
    "sklearn.decomposition": [
        "DictionaryLearning",
        "FactorAnalysis",
        "FastICA",
        "IncrementalPCA",
        "KernelPCA",
        "LatentDirichletAllocation",
        "MiniBatchDictionaryLearning",
        "MiniBatchNMF",
        "MiniBatchSparsePCA",
        "NMF",
        "PCA",
        "SparseCoder",
        "SparsePCA",
        "TruncatedSVD",
    ],
    "sklearn.discriminant_analysis": [
        "LinearDiscriminantAnalysis",
        "QuadraticDiscriminantAnalysis",
    ],
    "sklearn.dummy": [
        "DummyClassifier",
        "DummyRegressor",
    ],
    "sklearn.ensemble": [
        "AdaBoostClassifier",
        "AdaBoostRegressor",
        "BaggingClassifier",
        "BaggingRegressor",
        "ExtraTreesClassifier",
        "ExtraTreesRegressor",
        "GradientBoostingClassifier",
        "GradientBoostingRegressor",
        "HistGradientBoostingClassifier",
        "HistGradientBoostingRegressor",
        "IsolationForest",
        "RandomForestClassifier",
        "RandomForestRegressor",
        "RandomTreesEmbedding",
        "StackingClassifier",
        "StackingRegressor",
        "VotingClassifier",
        "VotingRegressor",
    ],
    "sklearn.feature_extraction": [
        "DictVectorizer",
        "FeatureHasher",
    ],
    "sklearn.feature_extraction.image": [
        "PatchExtractor",
    ],
    "sklearn.feature_extraction.text": [
        "CountVectorizer",
        "HashingVectorizer",
        "TfidfTransformer",
        "TfidfVectorizer",
    ],
    "sklearn.feature_selection": [
        "GenericUnivariateSelect",
        "RFE",
        "RFECV",
        "SelectFdr",
        "SelectFpr",
        "SelectFromModel",
        "SelectFwe",
        "SelectKBest",
        "SelectPercentile",
        "SequentialFeatureSelector",
        "VarianceThreshold",
    ],
    "sklearn.frozen": [
        "FrozenEstimator",
    ],
    "sklearn.gaussian_process": [
        "GaussianProcessClassifier",
        "GaussianProcessRegressor",
    ],
    "sklearn.impute": [
        "KNNImputer",
        "MissingIndicator",
        "SimpleImputer",
    ],
    "sklearn.isotonic": [
        "IsotonicRegression",
    ],
    "sklearn.kernel_approximation": [
        "AdditiveChi2Sampler",
        "Nystroem",
        "PolynomialCountSketch",
        "RBFSampler",
        "SkewedChi2Sampler",
    ],
    "sklearn.kernel_ridge": [
        "KernelRidge",
    ],
    "sklearn.linear_model": [
        "ARDRegression",
        "BayesianRidge",
        "ElasticNet",
        "ElasticNetCV",
        "GammaRegressor",
        "HuberRegressor",
        "Lars",
        "LarsCV",
        "Lasso",
        "LassoCV",
        "LassoLars",
        "LassoLarsCV",
        "LassoLarsIC",
        "LinearRegression",
        "LogisticRegression",
        "LogisticRegressionCV",
        "MultiTaskElasticNet",
        "MultiTaskElasticNetCV",
        "MultiTaskLasso",
        "MultiTaskLassoCV",
        "OrthogonalMatchingPursuit",
        "OrthogonalMatchingPursuitCV",
        "PassiveAggressiveClassifier",
        "PassiveAggressiveRegressor",
        "Perceptron",
        "PoissonRegressor",
        "QuantileRegressor",
        "RANSACRegressor",
        "Ridge",
        "RidgeCV",
        "RidgeClassifier",
        "RidgeClassifierCV",
        "SGDClassifier",
        "SGDOneClassSVM",
        "SGDRegressor",
        "TheilSenRegressor",
        "TweedieRegressor",
    ],
    "sklearn.manifold": [
        "ClassicalMDS",
        "Isomap",
        "LocallyLinearEmbedding",
        "MDS",
        "SpectralEmbedding",
        "TSNE",
    ],
    "sklearn.mixture": [
        "BayesianGaussianMixture",
        "GaussianMixture",
    ],
    "sklearn.model_selection": [
        "FixedThresholdClassifier",
        "GridSearchCV",
        "RandomizedSearchCV",
        "TunedThresholdClassifierCV",
    ],
    "sklearn.multiclass": [
        "OneVsOneClassifier",
        "OneVsRestClassifier",
        "OutputCodeClassifier",
    ],
    "sklearn.multioutput": [
        "ClassifierChain",
        "MultiOutputClassifier",
        "MultiOutputRegressor",
        "RegressorChain",
    ],
    "sklearn.naive_bayes": [
        "BernoulliNB",
        "CategoricalNB",
        "ComplementNB",
        "GaussianNB",
        "MultinomialNB",
    ],
    "sklearn.neighbors": [
        "KNeighborsClassifier",
        "KNeighborsRegressor",
        "KNeighborsTransformer",
        "KernelDensity",
        "LocalOutlierFactor",
        "NearestCentroid",
        "NearestNeighbors",
        "NeighborhoodComponentsAnalysis",
        "RadiusNeighborsClassifier",
        "RadiusNeighborsRegressor",
        "RadiusNeighborsTransformer",
    ],
    "sklearn.neural_network": [
        "BernoulliRBM",
        "MLPClassifier",
        "MLPRegressor",
    ],
    "sklearn.pipeline": [
        "FeatureUnion",
        "Pipeline",
    ],
    "sklearn.preprocessing": [
        "Binarizer",
        "FunctionTransformer",
        "KBinsDiscretizer",
        "KernelCenterer",
        "LabelBinarizer",
        "LabelEncoder",
        "MaxAbsScaler",
        "MinMaxScaler",
        "MultiLabelBinarizer",
        "Normalizer",
        "OneHotEncoder",
        "OrdinalEncoder",
        "PolynomialFeatures",
        "PowerTransformer",
        "QuantileTransformer",
        "RobustScaler",
        "SplineTransformer",
        "StandardScaler",
        "TargetEncoder",
    ],
    "sklearn.random_projection": [
        "GaussianRandomProjection",
        "SparseRandomProjection",
    ],
    "sklearn.semi_supervised": [
        "LabelPropagation",
        "LabelSpreading",
        "SelfTrainingClassifier",
    ],
    "sklearn.svm": [
        "LinearSVC",
        "LinearSVR",
        "NuSVC",
        "NuSVR",
        "OneClassSVM",
        "SVC",
        "SVR",
    ],
    "sklearn.tree": [
        "DecisionTreeClassifier",
        "DecisionTreeRegressor",
        "ExtraTreeClassifier",
        "ExtraTreeRegressor",
    ],
}

# Reverse mapping from an estimator name to the public submodule it lives in,
# used by `__getattr__` to lazily import the estimator on first access.
_estimator_to_module = {
    estimator: module
    for module, estimators in _estimators.items()
    for estimator in estimators
}

__all__ = (
    _submodules
    + sorted(_estimator_to_module)
    + [
        # Non-modules:
        "clone",
        "get_config",
        "set_config",
        "config_context",
        "show_versions",
    ]
)


def __dir__():
    return __all__


def __getattr__(name):
    if name in _submodules:
        return _importlib.import_module(f"sklearn.{name}")
    elif name in _estimator_to_module:
        module = _importlib.import_module(_estimator_to_module[name])
        estimator = getattr(module, name)
        # Cache on the module so subsequent accesses bypass `__getattr__`.
        globals()[name] = estimator
        return estimator
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(f"Module 'sklearn' has no attribute '{name}'")


def setup_module(module):
    """Fixture for the tests to assure globally controllable seeding of RNGs"""

    import numpy as np

    # Check if a random seed exists in the environment, if not create one.
    _random_seed = os.environ.get("SKLEARN_SEED", None)
    if _random_seed is None:
        _random_seed = np.random.uniform() * np.iinfo(np.int32).max
    _random_seed = int(_random_seed)
    print("I: Seeding RNGs with %r" % _random_seed)
    np.random.seed(_random_seed)
    random.seed(_random_seed)
