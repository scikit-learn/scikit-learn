"""Configuration for the API reference documentation.

CONFIGURING AN AUTOSUMMARY BLOCK
================================

An autosummary block is configured as a list of dictionaries, each containing "template"
that maps to a template name available under `doc/templates/` or `None` if no custom
template is needed, and "entries" that maps to a list of functions/classes/etc. that
uses the corresponding template. For instance, we have the following translation:

|---------------------------------------------| |--------------------------------------|
|  [                                          | |                                      |
|      {                                      | |  .. autosummary::                    |
|          "template": "class",               | |     :toctree: ../modules/generated/  |
|          "entries": [                       | |     :template: class.rst             |
|              "ColumnTransformer",           | |                                      |
|              "TransformedTargetRegressor",  | |     ColumnTransformer                |
|          ],                                 | |     TransformedTargetRegressor       |
|      },                                     | |                                      |
|      {                                      | |     :template: function.rst          |
|          "template": "function",            | |                                      |
|          "entries": [                       | |     make_column_selector             |
|              "make_column_transformer",     | |     make_column_transformer          |
|              "make_column_selector",        | |                                      |
|          ],                                 | |                                      |
|      },                                     | |                                      |
|  ],                                         | |                                      |
|---------------------------------------------| |--------------------------------------|

CONFIGURING API_REFERENCE
=========================

API_REFERENCE maps each module name to a dictionary that consists of the following
components:

short_summary (required)
    The text to be printed on the index page; it has nothing to do the API reference
    page of each module.
description (required, `None` if not needed)
    The additional description for the module to be placed under the module
    docstring, before the sections start.
sections (required)
    A list of sections, each of which consists of:
    - title (required, `None` if not needed): the section title, commonly it should
      not be `None` except for the first section of a module,
    - description (optional): the optional additional description for the section,
    - autosummary (required): an autosummary block, assuming current module is the
      current module name.

Essentially, the rendered page would look like the following:

|---------------------------------------------------------------------------------|
|     {{ module_name }}                                                           |
|     =================                                                           |
|     {{ module_docstring }}                                                      |
|     {{ description }}                                                           |
|                                                                                 |
|     {{ section_title_1 }}   <-------------- Optional if one wants the first     |
|     ---------------------                   section to directly follow          |
|     {{ section_description_1 }}             without a second-level heading.     |
|     {{ section_autosummary_1 }}                                                 |
|                                                                                 |
|     {{ section_title_2 }}                                                       |
|     ---------------------                                                       |
|     {{ section_description_2 }}                                                 |
|     {{ section_autosummary_2 }}                                                 |
|                                                                                 |
|     More sections...                                                            |
|---------------------------------------------------------------------------------|

CONFIGURING DEPRECATED_API_REFERENCE
====================================

DEPRECATED_API_REFERENCE maps each deprecation target version to a corresponding
autosummary block. It will be placed at the bottom of the API index page under the
"Recently deprecated" section. Essentially, the rendered section would look like the
following:

|------------------------------------------|
|     To be removed in {{ version_1 }}     |
|     --------------------------------     |
|     {{ autosummary_1 }}                  |
|                                          |
|     To be removed in {{ version_2 }}     |
|     --------------------------------     |
|     {{ autosummary_2 }}                  |
|                                          |
|     More versions...                     |
|------------------------------------------|

Note that the autosummary here assumes that the current module is `sklearn`, i.e., if
`sklearn.utils.Memory` is deprecated, one should put `utils.Memory` in the "entries"
slot of the autosummary block.
"""

from io import StringIO


def _get_guide(*refs, is_developer=False):
    if len(refs) == 1:
        ref_desc = f":ref:`{refs[0]}` section"
    elif len(refs) == 2:
        ref_desc = f":ref:`{refs[0]}` and :ref:`{refs[1]}` sections"
    else:
        ref_desc = ", ".join(f":ref:`{ref}`" for ref in refs[:-1])
        ref_desc += f", and :ref:`{refs[-1]}` sections"

    guide_name = "Developer" if is_developer else "User"
    return f"**{guide_name} guide.** See the {ref_desc} for further details."


API_REFERENCE = {
    "sklearn": {
        "short_summary": "Settings and information tools.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "config_context",
                            "get_config",
                            "set_config",
                            "show_versions",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.base": {
        "short_summary": "Base classes and utility functions.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "BaseEstimator",
                            "BiclusterMixin",
                            "ClassifierMixin",
                            "ClusterMixin",
                            "DensityMixin",
                            "RegressorMixin",
                            "TransformerMixin",
                            "MetaEstimatorMixin",
                            "OneToOneFeatureMixin",
                            "OutlierMixin",
                            "ClassNamePrefixFeaturesOutMixin",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "clone",
                            "is_classifier",
                            "is_regressor",
                        ],
                    },
                ],
            }
        ],
    },
    "sklearn.calibration": {
        "short_summary": "Probability calibration.",
        "description": _get_guide("calibration"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["CalibratedClassifierCV"],
                    },
                    {
                        "template": "function",
                        "entries": ["calibration_curve"],
                    },
                ],
            },
            {
                "title": "Visualization",
                "autosummary": [
                    {
                        "template": "display_all_class_methods",
                        "entries": ["CalibrationDisplay"],
                    },
                ],
            },
        ],
    },
    "sklearn.cluster": {
        "short_summary": "Clustering.",
        "description": _get_guide("clustering", "biclustering"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "AffinityPropagation",
                            "AgglomerativeClustering",
                            "Birch",
                            "DBSCAN",
                            "HDBSCAN",
                            "FeatureAgglomeration",
                            "KMeans",
                            "BisectingKMeans",
                            "MiniBatchKMeans",
                            "MeanShift",
                            "OPTICS",
                            "SpectralClustering",
                            "SpectralBiclustering",
                            "SpectralCoclustering",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "affinity_propagation",
                            "cluster_optics_dbscan",
                            "cluster_optics_xi",
                            "compute_optics_graph",
                            "dbscan",
                            "estimate_bandwidth",
                            "k_means",
                            "kmeans_plusplus",
                            "mean_shift",
                            "spectral_clustering",
                            "ward_tree",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.compose": {
        "short_summary": "Composite estimators.",
        "description": _get_guide("combining_estimators"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "ColumnTransformer",
                            "TransformedTargetRegressor",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "make_column_transformer",
                            "make_column_selector",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.covariance": {
        "short_summary": "Covariance estimation.",
        "description": _get_guide("covariance"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "EmpiricalCovariance",
                            "EllipticEnvelope",
                            "GraphicalLasso",
                            "GraphicalLassoCV",
                            "LedoitWolf",
                            "MinCovDet",
                            "OAS",
                            "ShrunkCovariance",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "empirical_covariance",
                            "graphical_lasso",
                            "ledoit_wolf",
                            "ledoit_wolf_shrinkage",
                            "oas",
                            "shrunk_covariance",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.cross_decomposition": {
        "short_summary": "Cross decomposition.",
        "description": _get_guide("cross_decomposition"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "CCA",
                            "PLSCanonical",
                            "PLSRegression",
                            "PLSSVD",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.datasets": {
        "short_summary": "Datasets.",
        "description": _get_guide("datasets"),
        "sections": [
            {
                "title": "Loaders",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "clear_data_home",
                            "dump_svmlight_file",
                            "fetch_20newsgroups",
                            "fetch_20newsgroups_vectorized",
                            "fetch_california_housing",
                            "fetch_covtype",
                            "fetch_kddcup99",
                            "fetch_lfw_pairs",
                            "fetch_lfw_people",
                            "fetch_olivetti_faces",
                            "fetch_openml",
                            "fetch_rcv1",
                            "fetch_species_distributions",
                            "get_data_home",
                            "load_breast_cancer",
                            "load_diabetes",
                            "load_digits",
                            "load_files",
                            "load_iris",
                            "load_linnerud",
                            "load_sample_image",
                            "load_sample_images",
                            "load_svmlight_file",
                            "load_svmlight_files",
                            "load_wine",
                        ],
                    },
                ],
            },
            {
                "title": "Sample generators",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "make_biclusters",
                            "make_blobs",
                            "make_checkerboard",
                            "make_circles",
                            "make_classification",
                            "make_friedman1",
                            "make_friedman2",
                            "make_friedman3",
                            "make_gaussian_quantiles",
                            "make_hastie_10_2",
                            "make_low_rank_matrix",
                            "make_moons",
                            "make_multilabel_classification",
                            "make_regression",
                            "make_s_curve",
                            "make_sparse_coded_signal",
                            "make_sparse_spd_matrix",
                            "make_sparse_uncorrelated",
                            "make_spd_matrix",
                            "make_swiss_roll",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.decomposition": {
        "short_summary": "Matrix decomposition.",
        "description": _get_guide("decompositions"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "DictionaryLearning",
                            "FactorAnalysis",
                            "FastICA",
                            "IncrementalPCA",
                            "KernelPCA",
                            "LatentDirichletAllocation",
                            "MiniBatchDictionaryLearning",
                            "MiniBatchSparsePCA",
                            "NMF",
                            "MiniBatchNMF",
                            "PCA",
                            "SparsePCA",
                            "SparseCoder",
                            "TruncatedSVD",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "dict_learning",
                            "dict_learning_online",
                            "fastica",
                            "non_negative_factorization",
                            "sparse_encode",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.discriminant_analysis": {
        "short_summary": "Discriminant analysis.",
        "description": _get_guide("lda_qda"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "LinearDiscriminantAnalysis",
                            "QuadraticDiscriminantAnalysis",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.dummy": {
        "short_summary": "Dummy estimators.",
        "description": _get_guide("model_evaluation"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "DummyClassifier",
                            "DummyRegressor",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.ensemble": {
        "short_summary": "Ensemble methods.",
        "description": _get_guide("ensemble"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "AdaBoostClassifier",
                            "AdaBoostRegressor",
                            "BaggingClassifier",
                            "BaggingRegressor",
                            "ExtraTreesClassifier",
                            "ExtraTreesRegressor",
                            "GradientBoostingClassifier",
                            "GradientBoostingRegressor",
                            "IsolationForest",
                            "RandomForestClassifier",
                            "RandomForestRegressor",
                            "RandomTreesEmbedding",
                            "StackingClassifier",
                            "StackingRegressor",
                            "VotingClassifier",
                            "VotingRegressor",
                            "HistGradientBoostingRegressor",
                            "HistGradientBoostingClassifier",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.exceptions": {
        "short_summary": "Exceptions and warnings.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "DataConversionWarning",
                            "ConvergenceWarning",
                            "DataDimensionalityWarning",
                            "EfficiencyWarning",
                            "FitFailedWarning",
                            "InconsistentVersionWarning",
                            "NotFittedError",
                            "UndefinedMetricWarning",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.experimental": {
        "short_summary": "Experimental tools.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": None,
                        "entries": [
                            "enable_iterative_imputer",
                            "enable_halving_search_cv",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.feature_extraction": {
        "short_summary": "Feature extraction.",
        "description": _get_guide("feature_extraction"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "DictVectorizer",
                            "FeatureHasher",
                        ],
                    },
                ],
            },
            {
                "title": "From images",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["image.PatchExtractor"],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "image.extract_patches_2d",
                            "image.grid_to_graph",
                            "image.img_to_graph",
                            "image.reconstruct_from_patches_2d",
                        ],
                    },
                ],
            },
            {
                "title": "From text",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "text.CountVectorizer",
                            "text.HashingVectorizer",
                            "text.TfidfTransformer",
                            "text.TfidfVectorizer",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.feature_selection": {
        "short_summary": "Feature selection.",
        "description": _get_guide("feature_selection"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "GenericUnivariateSelect",
                            "SelectPercentile",
                            "SelectKBest",
                            "SelectFpr",
                            "SelectFdr",
                            "SelectFromModel",
                            "SelectFwe",
                            "SequentialFeatureSelector",
                            "RFE",
                            "RFECV",
                            "VarianceThreshold",
                            "SelectorMixin",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "chi2",
                            "f_classif",
                            "f_regression",
                            "r_regression",
                            "mutual_info_classif",
                            "mutual_info_regression",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.gaussian_process": {
        "short_summary": "Gaussian processes.",
        "description": _get_guide("gaussian_process"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "GaussianProcessClassifier",
                            "GaussianProcessRegressor",
                        ],
                    },
                ],
            },
            {
                "title": "Kernels",
                "autosummary": [
                    {
                        "template": "class_with_call",
                        "entries": [
                            "kernels.CompoundKernel",
                            "kernels.ConstantKernel",
                            "kernels.DotProduct",
                            "kernels.ExpSineSquared",
                            "kernels.Exponentiation",
                            "kernels.Hyperparameter",
                            "kernels.Kernel",
                            "kernels.Matern",
                            "kernels.PairwiseKernel",
                            "kernels.Product",
                            "kernels.RBF",
                            "kernels.RationalQuadratic",
                            "kernels.Sum",
                            "kernels.WhiteKernel",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.impute": {
        "short_summary": "Imputation.",
        "description": _get_guide("impute"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "SimpleImputer",
                            "IterativeImputer",
                            "MissingIndicator",
                            "KNNImputer",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.inspection": {
        "short_summary": "Inspection.",
        "description": _get_guide("inspection"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "partial_dependence",
                            "permutation_importance",
                        ],
                    },
                ],
            },
            {
                "title": "Plotting",
                "autosummary": [
                    {
                        "template": "display_only_from_estimator",
                        "entries": [
                            "DecisionBoundaryDisplay",
                            "PartialDependenceDisplay",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.isotonic": {
        "short_summary": "Isotonic regression.",
        "description": _get_guide("isotonic"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["IsotonicRegression"],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "check_increasing",
                            "isotonic_regression",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.kernel_approximation": {
        "short_summary": "Isotonic regression.",
        "description": _get_guide("kernel_approximation"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "AdditiveChi2Sampler",
                            "Nystroem",
                            "PolynomialCountSketch",
                            "RBFSampler",
                            "SkewedChi2Sampler",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.kernel_ridge": {
        "short_summary": "Kernel ridge regression.",
        "description": _get_guide("kernel_ridge"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["KernelRidge"],
                    },
                ],
            },
        ],
    },
    "sklearn.linear_model": {
        "short_summary": "Generalized linear models.",
        "description": (
            _get_guide("linear_model")
            + "\n\nThe following subsections are only rough guidelines: the same "
            "estimator can fall into multiple categories, depending on its parameters."
        ),
        "sections": [
            {
                "title": "Linear classifiers",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "LogisticRegression",
                            "LogisticRegressionCV",
                            "PassiveAggressiveClassifier",
                            "Perceptron",
                            "RidgeClassifier",
                            "RidgeClassifierCV",
                            "SGDClassifier",
                            "SGDOneClassSVM",
                        ],
                    },
                ],
            },
            {
                "title": "Classical linear regressors",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "LinearRegression",
                            "Ridge",
                            "RidgeCV",
                            "SGDRegressor",
                        ],
                    },
                ],
            },
            {
                "title": "Regressors with variable selection",
                "description": (
                    "The following estimators have built-in variable selection fitting "
                    "procedures, but any estimator using a L1 or elastic-net penalty "
                    "also performs variable selection: typically "
                    ":class:`~linear_model.SGDRegressor` or "
                    ":class:`~sklearn.linear_model.SGDClassifier` with an appropriate "
                    "penalty."
                ),
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "ElasticNet",
                            "ElasticNetCV",
                            "Lars",
                            "LarsCV",
                            "Lasso",
                            "LassoCV",
                            "LassoLars",
                            "LassoLarsCV",
                            "LassoLarsIC",
                            "OrthogonalMatchingPursuit",
                            "OrthogonalMatchingPursuitCV",
                        ],
                    },
                ],
            },
            {
                "title": "Bayesian regressors",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "ARDRegression",
                            "BayesianRidge",
                        ],
                    },
                ],
            },
            {
                "title": "Multi-task linear regressors with variable selection",
                "description": (
                    "These estimators fit multiple regression problems (or tasks)"
                    " jointly, while inducing sparse coefficients. While the inferred"
                    " coefficients may differ between the tasks, they are constrained"
                    " to agree on the features that are selected (non-zero"
                    " coefficients)."
                ),
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "MultiTaskElasticNet",
                            "MultiTaskElasticNetCV",
                            "MultiTaskLasso",
                            "MultiTaskLassoCV",
                        ],
                    },
                ],
            },
            {
                "title": "Outlier-robust regressors",
                "description": (
                    "Any estimator using the Huber loss would also be robust to "
                    "outliers, e.g., :class:`~linear_model.SGDRegressor` with "
                    "``loss='huber'``."
                ),
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "HuberRegressor",
                            "QuantileRegressor",
                            "RANSACRegressor",
                            "TheilSenRegressor",
                        ],
                    },
                ],
            },
            {
                "title": "Generalized linear models (GLM) for regression",
                "description": (
                    "These models allow for response variables to have error "
                    "distributions other than a normal distribution."
                ),
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "PoissonRegressor",
                            "TweedieRegressor",
                            "GammaRegressor",
                        ],
                    },
                ],
            },
            {
                "title": "Miscellaneous",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["PassiveAggressiveRegressor"],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "enet_path",
                            "lars_path",
                            "lars_path_gram",
                            "lasso_path",
                            "orthogonal_mp",
                            "orthogonal_mp_gram",
                            "ridge_regression",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.manifold": {
        "short_summary": "Manifold learning.",
        "description": _get_guide("manifold"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "Isomap",
                            "LocallyLinearEmbedding",
                            "MDS",
                            "SpectralEmbedding",
                            "TSNE",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "locally_linear_embedding",
                            "smacof",
                            "spectral_embedding",
                            "trustworthiness",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.metrics": {
        "short_summary": "Metrics.",
        "description": _get_guide("model_evaluation", "metrics"),
        "sections": [
            {
                "title": "Model selection interface",
                "description": _get_guide("scoring_parameter"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "check_scoring",
                            "get_scorer",
                            "get_scorer_names",
                            "make_scorer",
                        ],
                    },
                ],
            },
            {
                "title": "Classification metrics",
                "description": _get_guide("classification_metrics"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "accuracy_score",
                            "auc",
                            "average_precision_score",
                            "balanced_accuracy_score",
                            "brier_score_loss",
                            "class_likelihood_ratios",
                            "classification_report",
                            "cohen_kappa_score",
                            "confusion_matrix",
                            "dcg_score",
                            "det_curve",
                            "f1_score",
                            "fbeta_score",
                            "hamming_loss",
                            "hinge_loss",
                            "jaccard_score",
                            "log_loss",
                            "matthews_corrcoef",
                            "multilabel_confusion_matrix",
                            "ndcg_score",
                            "precision_recall_curve",
                            "precision_recall_fscore_support",
                            "precision_score",
                            "recall_score",
                            "roc_auc_score",
                            "roc_curve",
                            "top_k_accuracy_score",
                            "zero_one_loss",
                        ],
                    },
                ],
            },
            {
                "title": "Regression metrics",
                "description": _get_guide("regression_metrics"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "explained_variance_score",
                            "max_error",
                            "mean_absolute_error",
                            "mean_squared_error",
                            "mean_squared_log_error",
                            "median_absolute_error",
                            "mean_absolute_percentage_error",
                            "r2_score",
                            "root_mean_squared_log_error",
                            "root_mean_squared_error",
                            "mean_poisson_deviance",
                            "mean_gamma_deviance",
                            "mean_tweedie_deviance",
                            "d2_tweedie_score",
                            "mean_pinball_loss",
                            "d2_pinball_score",
                            "d2_absolute_error_score",
                        ],
                    },
                ],
            },
            {
                "title": "Multilabel ranking metrics",
                "description": _get_guide("multilabel_ranking_metrics"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "coverage_error",
                            "label_ranking_average_precision_score",
                            "label_ranking_loss",
                        ],
                    },
                ],
            },
            {
                "title": "Clustering metrics",
                "description": (
                    "There are two forms of evaluation for cluster analysis results:\n"
                    "\n- supervised, which uses a ground truth class values for each "
                    "sample,\n- unsupervised, which does not use a ground truth and "
                    "measures the 'quality' of the model itself.\n\n"
                    + _get_guide("clustering_evaluation")
                ),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "adjusted_mutual_info_score",
                            "adjusted_rand_score",
                            "calinski_harabasz_score",
                            "davies_bouldin_score",
                            "completeness_score",
                            "cluster.contingency_matrix",
                            "cluster.pair_confusion_matrix",
                            "fowlkes_mallows_score",
                            "homogeneity_completeness_v_measure",
                            "homogeneity_score",
                            "mutual_info_score",
                            "normalized_mutual_info_score",
                            "rand_score",
                            "silhouette_score",
                            "silhouette_samples",
                            "v_measure_score",
                        ],
                    },
                ],
            },
            {
                "title": "Biclustering metrics",
                "description": _get_guide("biclustering_evaluation"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["consensus_score"],
                    },
                ],
            },
            {
                "title": "Distance metrics",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["DistanceMetric"],
                    },
                ],
            },
            {
                "title": "Pairwise metrics",
                "description": _get_guide("metrics"),
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "pairwise.additive_chi2_kernel",
                            "pairwise.chi2_kernel",
                            "pairwise.cosine_similarity",
                            "pairwise.cosine_distances",
                            "pairwise.distance_metrics",
                            "pairwise.euclidean_distances",
                            "pairwise.haversine_distances",
                            "pairwise.kernel_metrics",
                            "pairwise.laplacian_kernel",
                            "pairwise.linear_kernel",
                            "pairwise.manhattan_distances",
                            "pairwise.nan_euclidean_distances",
                            "pairwise.pairwise_kernels",
                            "pairwise.polynomial_kernel",
                            "pairwise.rbf_kernel",
                            "pairwise.sigmoid_kernel",
                            "pairwise.paired_euclidean_distances",
                            "pairwise.paired_manhattan_distances",
                            "pairwise.paired_cosine_distances",
                            "pairwise.paired_distances",
                            "pairwise_distances",
                            "pairwise_distances_argmin",
                            "pairwise_distances_argmin_min",
                            "pairwise_distances_chunked",
                        ],
                    },
                ],
            },
            {
                "title": "Plotting",
                "description": _get_guide("visualizations"),
                "autosummary": [
                    {
                        "template": "display_all_class_methods",
                        "entries": [
                            "ConfusionMatrixDisplay",
                            "DetCurveDisplay",
                            "PrecisionRecallDisplay",
                            "PredictionErrorDisplay",
                            "RocCurveDisplay",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.mixture": {
        "short_summary": "Gaussian mixture models.",
        "description": _get_guide("mixture"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "BayesianGaussianMixture",
                            "GaussianMixture",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.model_selection": {
        "short_summary": "Model selection.",
        "description": _get_guide("cross_validation", "grid_search", "learning_curve"),
        "sections": [
            {
                "title": "Splitters",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "GroupKFold",
                            "GroupShuffleSplit",
                            "KFold",
                            "LeaveOneGroupOut",
                            "LeavePGroupsOut",
                            "LeaveOneOut",
                            "LeavePOut",
                            "PredefinedSplit",
                            "RepeatedKFold",
                            "RepeatedStratifiedKFold",
                            "ShuffleSplit",
                            "StratifiedKFold",
                            "StratifiedShuffleSplit",
                            "StratifiedGroupKFold",
                            "TimeSeriesSplit",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "check_cv",
                            "train_test_split",
                        ],
                    },
                ],
            },
            {
                "title": "Hyper-parameter optimizers",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "GridSearchCV",
                            "HalvingGridSearchCV",
                            "ParameterGrid",
                            "ParameterSampler",
                            "RandomizedSearchCV",
                            "HalvingRandomSearchCV",
                        ],
                    },
                ],
            },
            {
                "title": "Model validation",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "cross_validate",
                            "cross_val_predict",
                            "cross_val_score",
                            "learning_curve",
                            "permutation_test_score",
                            "validation_curve",
                        ],
                    },
                ],
            },
            {
                "title": "Visualization",
                "autosummary": [
                    {
                        "template": "display_only_from_estimator",
                        "entries": [
                            "LearningCurveDisplay",
                            "ValidationCurveDisplay",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.multiclass": {
        "short_summary": "Multiclass classification.",
        "description": _get_guide("multiclass_classification"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "OneVsRestClassifier",
                            "OneVsOneClassifier",
                            "OutputCodeClassifier",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.multioutput": {
        "short_summary": "Multioutput regression and classification.",
        "description": _get_guide(
            "multilabel_classification",
            "multiclass_multioutput_classification",
            "multioutput_regression",
        ),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "ClassifierChain",
                            "MultiOutputRegressor",
                            "MultiOutputClassifier",
                            "RegressorChain",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.naive_bayes": {
        "short_summary": "Naive Bayes.",
        "description": _get_guide("naive_bayes"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "BernoulliNB",
                            "CategoricalNB",
                            "ComplementNB",
                            "GaussianNB",
                            "MultinomialNB",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.neighbors": {
        "short_summary": "Nearest neighbors.",
        "description": _get_guide("neighbors"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "BallTree",
                            "KDTree",
                            "KernelDensity",
                            "KNeighborsClassifier",
                            "KNeighborsRegressor",
                            "KNeighborsTransformer",
                            "LocalOutlierFactor",
                            "RadiusNeighborsClassifier",
                            "RadiusNeighborsRegressor",
                            "RadiusNeighborsTransformer",
                            "NearestCentroid",
                            "NearestNeighbors",
                            "NeighborhoodComponentsAnalysis",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "kneighbors_graph",
                            "radius_neighbors_graph",
                            "sort_graph_by_row_values",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.neural_network": {
        "short_summary": "Neural network models.",
        "description": _get_guide(
            "neural_networks_supervised", "neural_networks_unsupervised"
        ),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "BernoulliRBM",
                            "MLPClassifier",
                            "MLPRegressor",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.pipeline": {
        "short_summary": "Pipeline.",
        "description": _get_guide("combining_estimators"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "FeatureUnion",
                            "Pipeline",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "make_pipeline",
                            "make_union",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.preprocessing": {
        "short_summary": "Preprocessing and normalization.",
        "description": _get_guide("preprocessing"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "Binarizer",
                            "FunctionTransformer",
                            "KBinsDiscretizer",
                            "KernelCenterer",
                            "LabelBinarizer",
                            "LabelEncoder",
                            "MultiLabelBinarizer",
                            "MaxAbsScaler",
                            "MinMaxScaler",
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
                    },
                    {
                        "template": "function",
                        "entries": [
                            "add_dummy_feature",
                            "binarize",
                            "label_binarize",
                            "maxabs_scale",
                            "minmax_scale",
                            "normalize",
                            "quantile_transform",
                            "robust_scale",
                            "scale",
                            "power_transform",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.random_projection": {
        "short_summary": "Random projection.",
        "description": _get_guide("random_projection"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "GaussianRandomProjection",
                            "SparseRandomProjection",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": ["johnson_lindenstrauss_min_dim"],
                    },
                ],
            },
        ],
    },
    "sklearn.semi_supervised": {
        "short_summary": "Semi-supervised learning.",
        "description": _get_guide("semi_supervised"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "LabelPropagation",
                            "LabelSpreading",
                            "SelfTrainingClassifier",
                        ],
                    },
                ],
            },
        ],
    },
    "sklearn.svm": {
        "short_summary": "Support vector machines.",
        "description": _get_guide("svm"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "LinearSVC",
                            "LinearSVR",
                            "NuSVC",
                            "NuSVR",
                            "OneClassSVM",
                            "SVC",
                            "SVR",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": ["l1_min_c"],
                    },
                ],
            },
        ],
    },
    "sklearn.tree": {
        "short_summary": "Decision trees.",
        "description": _get_guide("tree"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "DecisionTreeClassifier",
                            "DecisionTreeRegressor",
                            "ExtraTreeClassifier",
                            "ExtraTreeRegressor",
                        ],
                    },
                ],
            },
            {
                "title": "Exporting",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "export_graphviz",
                            "export_text",
                        ],
                    },
                ],
            },
            {
                "title": "Plotting",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["plot_tree"],
                    },
                ],
            },
        ],
    },
    "sklearn.utils": {
        "short_summary": "Utilities.",
        "description": _get_guide("developers-utils", is_developer=True),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["Bunch"],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "as_float_array",
                            "assert_all_finite",
                            "deprecated",
                            "estimator_html_repr",
                            "gen_batches",
                            "gen_even_slices",
                            "indexable",
                            "murmurhash3_32",
                            "resample",
                            "_safe_indexing",
                            "safe_mask",
                            "safe_sqr",
                            "shuffle",
                        ],
                    },
                ],
            },
            {
                "title": "Input and parameter validation",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "check_X_y",
                            "check_array",
                            "check_scalar",
                            "check_consistent_length",
                            "check_random_state",
                            "validation.check_is_fitted",
                            "validation.check_memory",
                            "validation.check_symmetric",
                            "validation.column_or_1d",
                            "validation.has_fit_parameter",
                        ],
                    },
                ],
            },
            {
                "title": "Meta-estimators",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["metaestimators.available_if"],
                    },
                ],
            },
            {
                "title": "Weight handling based on class labels",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "class_weight.compute_class_weight",
                            "class_weight.compute_sample_weight",
                        ],
                    },
                ],
            },
            {
                "title": "Dealing with multiclass target in classifiers",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "multiclass.type_of_target",
                            "multiclass.is_multilabel",
                            "multiclass.unique_labels",
                        ],
                    },
                ],
            },
            {
                "title": "Optimal mathematical operations",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "extmath.safe_sparse_dot",
                            "extmath.randomized_range_finder",
                            "extmath.randomized_svd",
                            "extmath.fast_logdet",
                            "extmath.density",
                            "extmath.weighted_mode",
                        ],
                    },
                ],
            },
            {
                "title": "Working with sparse matrices and arrays",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "sparsefuncs.incr_mean_variance_axis",
                            "sparsefuncs.inplace_column_scale",
                            "sparsefuncs.inplace_row_scale",
                            "sparsefuncs.inplace_swap_row",
                            "sparsefuncs.inplace_swap_column",
                            "sparsefuncs.mean_variance_axis",
                            "sparsefuncs.inplace_csr_column_scale",
                            "sparsefuncs_fast.inplace_csr_row_normalize_l1",
                            "sparsefuncs_fast.inplace_csr_row_normalize_l2",
                        ],
                    },
                ],
            },
            {
                "title": "Working with graphs",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["graph.single_source_shortest_path_length"],
                    },
                ],
            },
            {
                "title": "Random sampling",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["random.sample_without_replacement"],
                    },
                ],
            },
            {
                "title": "Auxiliary functions that operate on arrays",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": ["arrayfuncs.min_pos"],
                    },
                ],
            },
            {
                "title": "Metadata routing",
                "description": _get_guide("metadata_routing"),
                "autosummary": [
                    {
                        "template": "class",
                        "entries": [
                            "metadata_routing.MetadataRouter",
                            "metadata_routing.MetadataRequest",
                            "metadata_routing.MethodMapping",
                        ],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "metadata_routing.get_routing_for_object",
                            "metadata_routing.process_routing",
                        ],
                    },
                ],
            },
            {
                "title": "Discovering scikit-learn objects",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "discovery.all_estimators",
                            "discovery.all_displays",
                            "discovery.all_functions",
                        ],
                    },
                ],
            },
            {
                "title": "API compatibility checkers",
                "autosummary": [
                    {
                        "template": "function",
                        "entries": [
                            "estimator_checks.check_estimator",
                            "estimator_checks.parametrize_with_checks",
                        ],
                    },
                ],
            },
            {
                "title": "Parallel computing",
                "autosummary": [
                    {
                        "template": "class",
                        "entries": ["parallel.Parallel"],
                    },
                    {
                        "template": "function",
                        "entries": [
                            "parallel.delayed",
                            "parallel_backend",
                            "register_parallel_backend",
                        ],
                    },
                ],
            },
        ],
    },
}


DEPRECATED_API_REFERENCE = {}  # type: ignore


def _write_autosummary_rst(autosummary, f):
    """Write the autosummary rst to a text stream."""
    f.write(".. autosummary::\n")
    f.write("   :toctree: ../modules/generated/\n")
    for autosummary_item in autosummary:
        if autosummary_item["template"] is not None:
            f.write(f"   :template: {autosummary_item['template']}.rst\n")
        f.write("\n")
        # Sort the entries in alphabetical order
        for entry in sorted(autosummary_item["entries"]):
            f.write(f"   {entry}\n")
        f.write("\n")


def get_api_reference_rst(module_name):
    """Get the API reference rst for a module."""
    output = StringIO()

    module_info = API_REFERENCE[module_name]
    if module_name == "sklearn":
        subname = None
    else:
        assert module_name.startswith("sklearn.")
        subname = module_name[8:]

    # Print the cross-reference hook
    if subname is not None:
        output.write(f".. _{subname}_ref:\n\n")

    # Print the top-level heading
    output.write(f":mod:`{module_name}`\n")
    output.write("=" * (len(module_name) + 7) + "\n\n")

    # Print the module docstring
    output.write(f".. automodule:: {module_name}\n")
    output.write("   :no-members:\n")
    output.write("   :no-inherited-members:\n\n")

    # Print the additional description if it exists
    if module_info["description"] is not None:
        output.write(module_info["description"] + "\n\n")

    for section in module_info["sections"]:
        # Print the cross-reference hook
        section_title = section["title"]
        if section_title is not None and subname is not None:
            section_refname = section_title.lower().replace(" ", "-")
            output.write(f".. _{subname}_ref-{section_refname}:\n\n")

        # Print the title if it exists
        if section_title is not None:
            output.write(section_title + "\n")
            output.write("-" * len(section_title) + "\n\n")

        # Print the additional description if it exists
        section_description = section.get("description", None)
        if section_description is not None:
            output.write(section_description + "\n\n")

        # Print the autosummary
        _write_autosummary_rst(section["autosummary"], output)

    return output.getvalue()


def get_deprecated_api_reference_rst(version):
    """Print the deprecated API reference for a version."""
    output = StringIO()

    output.write(f"To be removed in {version}\n")
    output.write("-" * (len(version) + 17) + "\n\n")

    # Print the autosummary
    _write_autosummary_rst(DEPRECATED_API_REFERENCE[version], output)

    return output.getvalue()
