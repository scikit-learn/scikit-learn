"""Configuration for the API reference documentation."""


def _get_guide(*refs, is_developer=False):
    """Get the rst to refer to user/developer guide.

    `refs` is several references that can be used in the :ref:`...` directive.
    """
    if len(refs) == 1:
        ref_desc = f":ref:`{refs[0]}` section"
    elif len(refs) == 2:
        ref_desc = f":ref:`{refs[0]}` and :ref:`{refs[1]}` sections"
    else:
        ref_desc = ", ".join(f":ref:`{ref}`" for ref in refs[:-1])
        ref_desc += f", and :ref:`{refs[-1]}` sections"

    guide_name = "Developer" if is_developer else "User"
    return f"**{guide_name} guide.** See the {ref_desc} for further details."


def _get_submodule(module_name, submodule_name):
    """Get the submodule docstring and automatically add the hook.

    `module_name` is e.g. `sklearn.feature_extraction`, and `submodule_name` is e.g.
    `image`, so we get the docstring and hook for `sklearn.feature_extraction.image`
    submodule. `module_name` is used to reset the current module because autosummary
    automatically changes the current module.
    """
    lines = [
        f".. automodule:: {module_name}.{submodule_name}",
        f".. currentmodule:: {module_name}",
    ]
    return "\n\n".join(lines)


"""
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

Hooks will be automatically generated for each module and each section. For a module,
e.g., `sklearn.feature_extraction`, the hook would be `feature_extraction_ref`; for a
section, e.g., "From text" under `sklearn.feature_extraction`, the hook would be
`feature_extraction_ref-from-text`. However, note that a better way is to refer using
the :mod: directive, e.g., :mod:`sklearn.feature_extraction` for the module and
:mod:`sklearn.feature_extraction.text` for the section. Only in case that a section
is not a particular submodule does the hook become useful, e.g., the "Loaders" section
under `sklearn.datasets`.
"""

API_REFERENCE = {
    "sklearn": {
        "short_summary": "Settings and information tools.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": [
                    "config_context",
                    "get_config",
                    "set_config",
                    "show_versions",
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
                    "BaseEstimator",
                    "BiclusterMixin",
                    "ClassNamePrefixFeaturesOutMixin",
                    "ClassifierMixin",
                    "ClusterMixin",
                    "DensityMixin",
                    "MetaEstimatorMixin",
                    "OneToOneFeatureMixin",
                    "OutlierMixin",
                    "RegressorMixin",
                    "TransformerMixin",
                    "clone",
                    "is_classifier",
                    "is_clusterer",
                    "is_regressor",
                    "is_outlier_detector",
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
                "autosummary": ["CalibratedClassifierCV", "calibration_curve"],
            },
            {
                "title": "Visualization",
                "autosummary": ["CalibrationDisplay"],
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
    "sklearn.compose": {
        "short_summary": "Composite estimators.",
        "description": _get_guide("combining_estimators"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    "ColumnTransformer",
                    "TransformedTargetRegressor",
                    "make_column_selector",
                    "make_column_transformer",
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
                    "EllipticEnvelope",
                    "EmpiricalCovariance",
                    "GraphicalLasso",
                    "GraphicalLassoCV",
                    "LedoitWolf",
                    "MinCovDet",
                    "OAS",
                    "ShrunkCovariance",
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
    "sklearn.cross_decomposition": {
        "short_summary": "Cross decomposition.",
        "description": _get_guide("cross_decomposition"),
        "sections": [
            {
                "title": None,
                "autosummary": ["CCA", "PLSCanonical", "PLSRegression", "PLSSVD"],
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
                    "clear_data_home",
                    "dump_svmlight_file",
                    "fetch_20newsgroups",
                    "fetch_20newsgroups_vectorized",
                    "fetch_california_housing",
                    "fetch_covtype",
                    "fetch_file",
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
            {
                "title": "Sample generators",
                "autosummary": [
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
    "sklearn.decomposition": {
        "short_summary": "Matrix decomposition.",
        "description": _get_guide("decompositions"),
        "sections": [
            {
                "title": None,
                "autosummary": [
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
                    "dict_learning",
                    "dict_learning_online",
                    "fastica",
                    "non_negative_factorization",
                    "sparse_encode",
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
                    "LinearDiscriminantAnalysis",
                    "QuadraticDiscriminantAnalysis",
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
                "autosummary": ["DummyClassifier", "DummyRegressor"],
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
                    "ConvergenceWarning",
                    "DataConversionWarning",
                    "DataDimensionalityWarning",
                    "EfficiencyWarning",
                    "FitFailedWarning",
                    "InconsistentVersionWarning",
                    "NotFittedError",
                    "UndefinedMetricWarning",
                    "EstimatorCheckFailedWarning",
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
                "autosummary": ["enable_halving_search_cv", "enable_iterative_imputer"],
            },
        ],
    },
    "sklearn.feature_extraction": {
        "short_summary": "Feature extraction.",
        "description": _get_guide("feature_extraction"),
        "sections": [
            {
                "title": None,
                "autosummary": ["DictVectorizer", "FeatureHasher"],
            },
            {
                "title": "From images",
                "description": _get_submodule("sklearn.feature_extraction", "image"),
                "autosummary": [
                    "image.PatchExtractor",
                    "image.extract_patches_2d",
                    "image.grid_to_graph",
                    "image.img_to_graph",
                    "image.reconstruct_from_patches_2d",
                ],
            },
            {
                "title": "From text",
                "description": _get_submodule("sklearn.feature_extraction", "text"),
                "autosummary": [
                    "text.CountVectorizer",
                    "text.HashingVectorizer",
                    "text.TfidfTransformer",
                    "text.TfidfVectorizer",
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
                    "GenericUnivariateSelect",
                    "RFE",
                    "RFECV",
                    "SelectFdr",
                    "SelectFpr",
                    "SelectFromModel",
                    "SelectFwe",
                    "SelectKBest",
                    "SelectPercentile",
                    "SelectorMixin",
                    "SequentialFeatureSelector",
                    "VarianceThreshold",
                    "chi2",
                    "f_classif",
                    "f_regression",
                    "mutual_info_classif",
                    "mutual_info_regression",
                    "r_regression",
                ],
            },
        ],
    },
    "sklearn.frozen": {
        "short_summary": "Frozen estimators.",
        "description": None,
        "sections": [
            {
                "title": None,
                "autosummary": ["FrozenEstimator"],
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
                    "GaussianProcessClassifier",
                    "GaussianProcessRegressor",
                ],
            },
            {
                "title": "Kernels",
                "description": _get_submodule("sklearn.gaussian_process", "kernels"),
                "autosummary": [
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
    "sklearn.impute": {
        "short_summary": "Imputation.",
        "description": _get_guide("impute"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    "IterativeImputer",
                    "KNNImputer",
                    "MissingIndicator",
                    "SimpleImputer",
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
                "autosummary": ["partial_dependence", "permutation_importance"],
            },
            {
                "title": "Plotting",
                "autosummary": ["DecisionBoundaryDisplay", "PartialDependenceDisplay"],
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
                    "IsotonicRegression",
                    "check_increasing",
                    "isotonic_regression",
                ],
            },
        ],
    },
    "sklearn.kernel_approximation": {
        "short_summary": "Kernel approximation.",
        "description": _get_guide("kernel_approximation"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    "AdditiveChi2Sampler",
                    "Nystroem",
                    "PolynomialCountSketch",
                    "RBFSampler",
                    "SkewedChi2Sampler",
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
                "autosummary": ["KernelRidge"],
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
            {
                "title": "Classical linear regressors",
                "autosummary": ["LinearRegression", "Ridge", "RidgeCV", "SGDRegressor"],
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
            {
                "title": "Bayesian regressors",
                "autosummary": ["ARDRegression", "BayesianRidge"],
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
                    "MultiTaskElasticNet",
                    "MultiTaskElasticNetCV",
                    "MultiTaskLasso",
                    "MultiTaskLassoCV",
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
                    "HuberRegressor",
                    "QuantileRegressor",
                    "RANSACRegressor",
                    "TheilSenRegressor",
                ],
            },
            {
                "title": "Generalized linear models (GLM) for regression",
                "description": (
                    "These models allow for response variables to have error "
                    "distributions other than a normal distribution."
                ),
                "autosummary": [
                    "GammaRegressor",
                    "PoissonRegressor",
                    "TweedieRegressor",
                ],
            },
            {
                "title": "Miscellaneous",
                "autosummary": [
                    "PassiveAggressiveRegressor",
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
    "sklearn.manifold": {
        "short_summary": "Manifold learning.",
        "description": _get_guide("manifold"),
        "sections": [
            {
                "title": None,
                "autosummary": [
                    "Isomap",
                    "LocallyLinearEmbedding",
                    "MDS",
                    "SpectralEmbedding",
                    "TSNE",
                    "locally_linear_embedding",
                    "smacof",
                    "spectral_embedding",
                    "trustworthiness",
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
                    "check_scoring",
                    "get_scorer",
                    "get_scorer_names",
                    "make_scorer",
                ],
            },
            {
                "title": "Classification metrics",
                "description": _get_guide("classification_metrics"),
                "autosummary": [
                    "accuracy_score",
                    "auc",
                    "average_precision_score",
                    "balanced_accuracy_score",
                    "brier_score_loss",
                    "class_likelihood_ratios",
                    "classification_report",
                    "cohen_kappa_score",
                    "confusion_matrix",
                    "d2_brier_score",
                    "d2_log_loss_score",
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
            {
                "title": "Regression metrics",
                "description": _get_guide("regression_metrics"),
                "autosummary": [
                    "d2_absolute_error_score",
                    "d2_pinball_score",
                    "d2_tweedie_score",
                    "explained_variance_score",
                    "max_error",
                    "mean_absolute_error",
                    "mean_absolute_percentage_error",
                    "mean_gamma_deviance",
                    "mean_pinball_loss",
                    "mean_poisson_deviance",
                    "mean_squared_error",
                    "mean_squared_log_error",
                    "mean_tweedie_deviance",
                    "median_absolute_error",
                    "r2_score",
                    "root_mean_squared_error",
                    "root_mean_squared_log_error",
                ],
            },
            {
                "title": "Multilabel ranking metrics",
                "description": _get_guide("multilabel_ranking_metrics"),
                "autosummary": [
                    "coverage_error",
                    "label_ranking_average_precision_score",
                    "label_ranking_loss",
                ],
            },
            {
                "title": "Clustering metrics",
                "description": (
                    _get_submodule("sklearn.metrics", "cluster")
                    + "\n\n"
                    + _get_guide("clustering_evaluation")
                ),
                "autosummary": [
                    "adjusted_mutual_info_score",
                    "adjusted_rand_score",
                    "calinski_harabasz_score",
                    "cluster.contingency_matrix",
                    "cluster.pair_confusion_matrix",
                    "completeness_score",
                    "davies_bouldin_score",
                    "fowlkes_mallows_score",
                    "homogeneity_completeness_v_measure",
                    "homogeneity_score",
                    "mutual_info_score",
                    "normalized_mutual_info_score",
                    "rand_score",
                    "silhouette_samples",
                    "silhouette_score",
                    "v_measure_score",
                ],
            },
            {
                "title": "Biclustering metrics",
                "description": _get_guide("biclustering_evaluation"),
                "autosummary": ["consensus_score"],
            },
            {
                "title": "Distance metrics",
                "autosummary": ["DistanceMetric"],
            },
            {
                "title": "Pairwise metrics",
                "description": (
                    _get_submodule("sklearn.metrics", "pairwise")
                    + "\n\n"
                    + _get_guide("metrics")
                ),
                "autosummary": [
                    "pairwise.additive_chi2_kernel",
                    "pairwise.chi2_kernel",
                    "pairwise.cosine_distances",
                    "pairwise.cosine_similarity",
                    "pairwise.distance_metrics",
                    "pairwise.euclidean_distances",
                    "pairwise.haversine_distances",
                    "pairwise.kernel_metrics",
                    "pairwise.laplacian_kernel",
                    "pairwise.linear_kernel",
                    "pairwise.manhattan_distances",
                    "pairwise.nan_euclidean_distances",
                    "pairwise.paired_cosine_distances",
                    "pairwise.paired_distances",
                    "pairwise.paired_euclidean_distances",
                    "pairwise.paired_manhattan_distances",
                    "pairwise.pairwise_kernels",
                    "pairwise.polynomial_kernel",
                    "pairwise.rbf_kernel",
                    "pairwise.sigmoid_kernel",
                    "pairwise_distances",
                    "pairwise_distances_argmin",
                    "pairwise_distances_argmin_min",
                    "pairwise_distances_chunked",
                ],
            },
            {
                "title": "Plotting",
                "description": _get_guide("visualizations"),
                "autosummary": [
                    "ConfusionMatrixDisplay",
                    "DetCurveDisplay",
                    "PrecisionRecallDisplay",
                    "PredictionErrorDisplay",
                    "RocCurveDisplay",
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
                "autosummary": ["BayesianGaussianMixture", "GaussianMixture"],
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
                    "GroupKFold",
                    "GroupShuffleSplit",
                    "KFold",
                    "LeaveOneGroupOut",
                    "LeaveOneOut",
                    "LeavePGroupsOut",
                    "LeavePOut",
                    "PredefinedSplit",
                    "RepeatedKFold",
                    "RepeatedStratifiedKFold",
                    "ShuffleSplit",
                    "StratifiedGroupKFold",
                    "StratifiedKFold",
                    "StratifiedShuffleSplit",
                    "TimeSeriesSplit",
                    "check_cv",
                    "train_test_split",
                ],
            },
            {
                "title": "Hyper-parameter optimizers",
                "autosummary": [
                    "GridSearchCV",
                    "HalvingGridSearchCV",
                    "HalvingRandomSearchCV",
                    "ParameterGrid",
                    "ParameterSampler",
                    "RandomizedSearchCV",
                ],
            },
            {
                "title": "Post-fit model tuning",
                "autosummary": [
                    "FixedThresholdClassifier",
                    "TunedThresholdClassifierCV",
                ],
            },
            {
                "title": "Model validation",
                "autosummary": [
                    "cross_val_predict",
                    "cross_val_score",
                    "cross_validate",
                    "learning_curve",
                    "permutation_test_score",
                    "validation_curve",
                ],
            },
            {
                "title": "Visualization",
                "autosummary": ["LearningCurveDisplay", "ValidationCurveDisplay"],
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
                    "OneVsOneClassifier",
                    "OneVsRestClassifier",
                    "OutputCodeClassifier",
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
                    "ClassifierChain",
                    "MultiOutputClassifier",
                    "MultiOutputRegressor",
                    "RegressorChain",
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
                    "BernoulliNB",
                    "CategoricalNB",
                    "ComplementNB",
                    "GaussianNB",
                    "MultinomialNB",
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
                    "BallTree",
                    "KDTree",
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
                    "kneighbors_graph",
                    "radius_neighbors_graph",
                    "sort_graph_by_row_values",
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
                "autosummary": ["BernoulliRBM", "MLPClassifier", "MLPRegressor"],
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
                    "FeatureUnion",
                    "Pipeline",
                    "make_pipeline",
                    "make_union",
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
                    "add_dummy_feature",
                    "binarize",
                    "label_binarize",
                    "maxabs_scale",
                    "minmax_scale",
                    "normalize",
                    "power_transform",
                    "quantile_transform",
                    "robust_scale",
                    "scale",
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
                    "GaussianRandomProjection",
                    "SparseRandomProjection",
                    "johnson_lindenstrauss_min_dim",
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
                    "LabelPropagation",
                    "LabelSpreading",
                    "SelfTrainingClassifier",
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
                    "LinearSVC",
                    "LinearSVR",
                    "NuSVC",
                    "NuSVR",
                    "OneClassSVM",
                    "SVC",
                    "SVR",
                    "l1_min_c",
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
                    "DecisionTreeClassifier",
                    "DecisionTreeRegressor",
                    "ExtraTreeClassifier",
                    "ExtraTreeRegressor",
                ],
            },
            {
                "title": "Exporting",
                "autosummary": ["export_graphviz", "export_text"],
            },
            {
                "title": "Plotting",
                "autosummary": ["plot_tree"],
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
                    "Bunch",
                    "_safe_indexing",
                    "as_float_array",
                    "assert_all_finite",
                    "deprecated",
                    "estimator_html_repr",
                    "gen_batches",
                    "gen_even_slices",
                    "indexable",
                    "murmurhash3_32",
                    "resample",
                    "safe_mask",
                    "safe_sqr",
                    "shuffle",
                    "Tags",
                    "InputTags",
                    "TargetTags",
                    "ClassifierTags",
                    "RegressorTags",
                    "TransformerTags",
                    "get_tags",
                ],
            },
            {
                "title": "Input and parameter validation",
                "description": _get_submodule("sklearn.utils", "validation"),
                "autosummary": [
                    "check_X_y",
                    "check_array",
                    "check_consistent_length",
                    "check_random_state",
                    "check_scalar",
                    "validation.check_is_fitted",
                    "validation.check_memory",
                    "validation.check_symmetric",
                    "validation.column_or_1d",
                    "validation.has_fit_parameter",
                    "validation.validate_data",
                ],
            },
            {
                "title": "Meta-estimators",
                "description": _get_submodule("sklearn.utils", "metaestimators"),
                "autosummary": ["metaestimators.available_if"],
            },
            {
                "title": "Weight handling based on class labels",
                "description": _get_submodule("sklearn.utils", "class_weight"),
                "autosummary": [
                    "class_weight.compute_class_weight",
                    "class_weight.compute_sample_weight",
                ],
            },
            {
                "title": "Dealing with multiclass target in classifiers",
                "description": _get_submodule("sklearn.utils", "multiclass"),
                "autosummary": [
                    "multiclass.is_multilabel",
                    "multiclass.type_of_target",
                    "multiclass.unique_labels",
                ],
            },
            {
                "title": "Optimal mathematical operations",
                "description": _get_submodule("sklearn.utils", "extmath"),
                "autosummary": [
                    "extmath.density",
                    "extmath.fast_logdet",
                    "extmath.randomized_range_finder",
                    "extmath.randomized_svd",
                    "extmath.safe_sparse_dot",
                    "extmath.weighted_mode",
                ],
            },
            {
                "title": "Working with sparse matrices and arrays",
                "description": _get_submodule("sklearn.utils", "sparsefuncs"),
                "autosummary": [
                    "sparsefuncs.incr_mean_variance_axis",
                    "sparsefuncs.inplace_column_scale",
                    "sparsefuncs.inplace_csr_column_scale",
                    "sparsefuncs.inplace_row_scale",
                    "sparsefuncs.inplace_swap_column",
                    "sparsefuncs.inplace_swap_row",
                    "sparsefuncs.mean_variance_axis",
                ],
            },
            {
                "title": None,
                "description": _get_submodule("sklearn.utils", "sparsefuncs_fast"),
                "autosummary": [
                    "sparsefuncs_fast.inplace_csr_row_normalize_l1",
                    "sparsefuncs_fast.inplace_csr_row_normalize_l2",
                ],
            },
            {
                "title": "Working with graphs",
                "description": _get_submodule("sklearn.utils", "graph"),
                "autosummary": ["graph.single_source_shortest_path_length"],
            },
            {
                "title": "Random sampling",
                "description": _get_submodule("sklearn.utils", "random"),
                "autosummary": ["random.sample_without_replacement"],
            },
            {
                "title": "Auxiliary functions that operate on arrays",
                "description": _get_submodule("sklearn.utils", "arrayfuncs"),
                "autosummary": ["arrayfuncs.min_pos"],
            },
            {
                "title": "Metadata routing",
                "description": (
                    _get_submodule("sklearn.utils", "metadata_routing")
                    + "\n\n"
                    + _get_guide("metadata_routing")
                ),
                "autosummary": [
                    "metadata_routing.MetadataRequest",
                    "metadata_routing.MetadataRouter",
                    "metadata_routing.MethodMapping",
                    "metadata_routing.get_routing_for_object",
                    "metadata_routing.process_routing",
                ],
            },
            {
                "title": "Discovering scikit-learn objects",
                "description": _get_submodule("sklearn.utils", "discovery"),
                "autosummary": [
                    "discovery.all_displays",
                    "discovery.all_estimators",
                    "discovery.all_functions",
                ],
            },
            {
                "title": "API compatibility checkers",
                "description": _get_submodule("sklearn.utils", "estimator_checks"),
                "autosummary": [
                    "estimator_checks.check_estimator",
                    "estimator_checks.parametrize_with_checks",
                    "estimator_checks.estimator_checks_generator",
                ],
            },
            {
                "title": "Parallel computing",
                "description": _get_submodule("sklearn.utils", "parallel"),
                "autosummary": [
                    "parallel.Parallel",
                    "parallel.delayed",
                ],
            },
        ],
    },
}


"""
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

Example:

DEPRECATED_API_REFERENCE = {
    "0.24": [
        "model_selection.fit_grid_point",
        "utils.safe_indexing",
    ],
}
"""

DEPRECATED_API_REFERENCE = {}  # type: ignore[var-annotated]
