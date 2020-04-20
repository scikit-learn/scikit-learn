"""Generates submodule to allow deprecation of submodules and keeping git
blame."""
from pathlib import Path
from contextlib import suppress

# TODO: Remove the whole file in 0.24

# This is a set of 4-tuples consisting of
# (new_module_name, deprecated_path, correct_import_path, importee)
# importee is used by test_import_deprecations to check for DeprecationWarnings
_DEPRECATED_MODULES = [
    ('_mocking', 'sklearn.utils.mocking', 'sklearn.utils',
     'MockDataFrame'),

    ('_bagging', 'sklearn.ensemble.bagging', 'sklearn.ensemble',
     'BaggingClassifier'),
    ('_base', 'sklearn.ensemble.base', 'sklearn.ensemble',
     'BaseEnsemble'),
    ('_forest', 'sklearn.ensemble.forest', 'sklearn.ensemble',
     'RandomForestClassifier'),
    ('_gb', 'sklearn.ensemble.gradient_boosting', 'sklearn.ensemble',
     'GradientBoostingClassifier'),
    ('_iforest', 'sklearn.ensemble.iforest', 'sklearn.ensemble',
     'IsolationForest'),
    ('_voting', 'sklearn.ensemble.voting', 'sklearn.ensemble',
     'VotingClassifier'),
    ('_weight_boosting', 'sklearn.ensemble.weight_boosting',
     'sklearn.ensemble', 'AdaBoostClassifier'),
    ('_classes', 'sklearn.tree.tree', 'sklearn.tree',
     'DecisionTreeClassifier'),
    ('_export', 'sklearn.tree.export', 'sklearn.tree', 'export_graphviz'),

    ('_rbm', 'sklearn.neural_network.rbm', 'sklearn.neural_network',
     'BernoulliRBM'),
    ('_multilayer_perceptron', 'sklearn.neural_network.multilayer_perceptron',
     'sklearn.neural_network', 'MLPClassifier'),

    ('_weight_vector', 'sklearn.utils.weight_vector', 'sklearn.utils',
     'WeightVector'),
    ('_seq_dataset', 'sklearn.utils.seq_dataset', 'sklearn.utils',
     'ArrayDataset32'),
    ('_fast_dict', 'sklearn.utils.fast_dict', 'sklearn.utils', 'IntFloatDict'),

    ('_affinity_propagation', 'sklearn.cluster.affinity_propagation_',
     'sklearn.cluster', 'AffinityPropagation'),
    ('_bicluster', 'sklearn.cluster.bicluster', 'sklearn.cluster',
     'SpectralBiclustering'),
    ('_birch', 'sklearn.cluster.birch', 'sklearn.cluster', 'Birch'),
    ('_dbscan', 'sklearn.cluster.dbscan_', 'sklearn.cluster', 'DBSCAN'),
    ('_agglomerative', 'sklearn.cluster.hierarchical', 'sklearn.cluster',
     'FeatureAgglomeration'),
    ('_kmeans', 'sklearn.cluster.k_means_', 'sklearn.cluster', 'KMeans'),
    ('_mean_shift', 'sklearn.cluster.mean_shift_', 'sklearn.cluster',
     'MeanShift'),
    ('_optics', 'sklearn.cluster.optics_', 'sklearn.cluster', 'OPTICS'),
    ('_spectral', 'sklearn.cluster.spectral', 'sklearn.cluster',
     'SpectralClustering'),

    ('_base', 'sklearn.mixture.base', 'sklearn.mixture', 'BaseMixture'),
    ('_gaussian_mixture', 'sklearn.mixture.gaussian_mixture',
     'sklearn.mixture', 'GaussianMixture'),
    ('_bayesian_mixture', 'sklearn.mixture.bayesian_mixture',
     'sklearn.mixture', 'BayesianGaussianMixture'),

    ('_empirical_covariance', 'sklearn.covariance.empirical_covariance_',
     'sklearn.covariance', 'EmpiricalCovariance'),
    ('_shrunk_covariance', 'sklearn.covariance.shrunk_covariance_',
     'sklearn.covariance', 'ShrunkCovariance'),
    ('_robust_covariance', 'sklearn.covariance.robust_covariance',
     'sklearn.covariance', 'MinCovDet'),
    ('_graph_lasso', 'sklearn.covariance.graph_lasso_',
     'sklearn.covariance', 'GraphicalLasso'),
    ('_elliptic_envelope', 'sklearn.covariance.elliptic_envelope',
     'sklearn.covariance', 'EllipticEnvelope'),

    ('_cca', 'sklearn.cross_decomposition.cca_',
     'sklearn.cross_decomposition', 'CCA'),
    ('_pls', 'sklearn.cross_decomposition.pls_',
     'sklearn.cross_decomposition', 'PLSSVD'),

    ('_base', 'sklearn.svm.base', 'sklearn.svm', 'BaseLibSVM'),
    ('_bounds', 'sklearn.svm.bounds', 'sklearn.svm', 'l1_min_c'),
    ('_classes', 'sklearn.svm.classes', 'sklearn.svm', 'SVR'),
    ('_libsvm', 'sklearn.svm.libsvm', 'sklearn.svm', 'fit'),
    ('_libsvm_sparse', 'sklearn.svm.libsvm_sparse', 'sklearn.svm',
     'set_verbosity_wrap'),
    ('_liblinear', 'sklearn.svm.liblinear', 'sklearn.svm', 'train_wrap'),

    ('_base', 'sklearn.decomposition.base', 'sklearn.decomposition',
     'BaseEstimator'),
    ('_dict_learning', 'sklearn.decomposition.dict_learning',
     'sklearn.decomposition', 'MiniBatchDictionaryLearning'),
    ('_cdnmf_fast', 'sklearn.decomposition.cdnmf_fast',
     'sklearn.decomposition', '__dict__'),
    ('_factor_analysis', 'sklearn.decomposition.factor_analysis',
     'sklearn.decomposition', 'FactorAnalysis'),
    ('_fastica', 'sklearn.decomposition.fastica_', 'sklearn.decomposition',
     'FastICA'),
    ('_incremental_pca', 'sklearn.decomposition.incremental_pca',
     'sklearn.decomposition', 'IncrementalPCA'),
    ('_kernel_pca', 'sklearn.decomposition.kernel_pca',
     'sklearn.decomposition', 'KernelPCA'),
    ('_nmf', 'sklearn.decomposition.nmf', 'sklearn.decomposition', 'NMF'),
    ('_lda', 'sklearn.decomposition.online_lda',
     'sklearn.decomposition', 'LatentDirichletAllocation'),
    ('_online_lda_fast', 'sklearn.decomposition.online_lda_fast',
     'sklearn.decomposition', 'mean_change'),
    ('_pca', 'sklearn.decomposition.pca', 'sklearn.decomposition', 'PCA'),
    ('_sparse_pca', 'sklearn.decomposition.sparse_pca',
     'sklearn.decomposition', 'SparsePCA'),
    ('_truncated_svd', 'sklearn.decomposition.truncated_svd',
     'sklearn.decomposition', 'TruncatedSVD'),

    ('_gpr', 'sklearn.gaussian_process.gpr', 'sklearn.gaussian_process',
     'GaussianProcessRegressor'),
    ('_gpc', 'sklearn.gaussian_process.gpc', 'sklearn.gaussian_process',
     'GaussianProcessClassifier'),

    ('_base', 'sklearn.datasets.base', 'sklearn.datasets', 'get_data_home'),
    ('_california_housing', 'sklearn.datasets.california_housing',
     'sklearn.datasets', 'fetch_california_housing'),
    ('_covtype', 'sklearn.datasets.covtype', 'sklearn.datasets',
     'fetch_covtype'),
    ('_kddcup99', 'sklearn.datasets.kddcup99', 'sklearn.datasets',
     'fetch_kddcup99'),
    ('_lfw', 'sklearn.datasets.lfw', 'sklearn.datasets',
     'fetch_lfw_people'),
    ('_olivetti_faces', 'sklearn.datasets.olivetti_faces', 'sklearn.datasets',
     'fetch_olivetti_faces'),
    ('_openml', 'sklearn.datasets.openml', 'sklearn.datasets', 'fetch_openml'),
    ('_rcv1', 'sklearn.datasets.rcv1', 'sklearn.datasets', 'fetch_rcv1'),
    ('_samples_generator', 'sklearn.datasets.samples_generator',
     'sklearn.datasets', 'make_classification'),
    ('_species_distributions', 'sklearn.datasets.species_distributions',
     'sklearn.datasets', 'fetch_species_distributions'),
    ('_svmlight_format_io', 'sklearn.datasets.svmlight_format',
     'sklearn.datasets', 'load_svmlight_file'),
    ('_twenty_newsgroups', 'sklearn.datasets.twenty_newsgroups',
     'sklearn.datasets', 'strip_newsgroup_header'),

    ('_dict_vectorizer', 'sklearn.feature_extraction.dict_vectorizer',
     'sklearn.feature_extraction', 'DictVectorizer'),
    ('_hash', 'sklearn.feature_extraction.hashing',
     'sklearn.feature_extraction', 'FeatureHasher'),
    ('_stop_words', 'sklearn.feature_extraction.stop_words',
     'sklearn.feature_extraction.text', 'ENGLISH_STOP_WORDS'),

    ('_base', 'sklearn.linear_model.base', 'sklearn.linear_model',
     'LinearRegression'),
    ('_cd_fast', 'sklearn.linear_model.cd_fast', 'sklearn.linear_model',
     'sparse_enet_coordinate_descent'),
    ('_bayes', 'sklearn.linear_model.bayes', 'sklearn.linear_model',
     'BayesianRidge'),
    ('_coordinate_descent', 'sklearn.linear_model.coordinate_descent',
     'sklearn.linear_model', 'Lasso'),
    ('_huber', 'sklearn.linear_model.huber', 'sklearn.linear_model',
     'HuberRegressor'),
    ('_least_angle', 'sklearn.linear_model.least_angle',
     'sklearn.linear_model', 'LassoLarsCV'),
    ('_logistic', 'sklearn.linear_model.logistic', 'sklearn.linear_model',
     'LogisticRegression'),
    ('_omp', 'sklearn.linear_model.omp', 'sklearn.linear_model',
     'OrthogonalMatchingPursuit'),
    ('_passive_aggressive', 'sklearn.linear_model.passive_aggressive',
     'sklearn.linear_model', 'PassiveAggressiveClassifier'),
    ('_perceptron', 'sklearn.linear_model.perceptron', 'sklearn.linear_model',
     'Perceptron'),
    ('_ransac', 'sklearn.linear_model.ransac', 'sklearn.linear_model',
     'RANSACRegressor'),
    ('_ridge', 'sklearn.linear_model.ridge', 'sklearn.linear_model',
     'Ridge'),
    ('_sag', 'sklearn.linear_model.sag', 'sklearn.linear_model',
     'get_auto_step_size'),
    ('_sag_fast', 'sklearn.linear_model.sag_fast', 'sklearn.linear_model',
     'MultinomialLogLoss64'),
    ('_sgd_fast', 'sklearn.linear_model.sgd_fast', 'sklearn.linear_model',
     'Hinge'),
    ('_stochastic_gradient', 'sklearn.linear_model.stochastic_gradient',
     'sklearn.linear_model', 'SGDClassifier'),
    ('_theil_sen', 'sklearn.linear_model.theil_sen', 'sklearn.linear_model',
     'TheilSenRegressor'),

    ('_bicluster', 'sklearn.metrics.cluster.bicluster',
     'sklearn.metrics.cluster', 'consensus_score'),
    ('_supervised', 'sklearn.metrics.cluster.supervised',
     'sklearn.metrics.cluster', 'entropy'),
    ('_unsupervised', 'sklearn.metrics.cluster.unsupervised',
     'sklearn.metrics.cluster', 'silhouette_score'),
    ('_expected_mutual_info_fast',
     'sklearn.metrics.cluster.expected_mutual_info_fast',
     'sklearn.metrics.cluster', 'expected_mutual_information'),

    ('_base', 'sklearn.metrics.base', 'sklearn.metrics', 'combinations'),
    ('_classification', 'sklearn.metrics.classification', 'sklearn.metrics',
     'accuracy_score'),
    ('_regression', 'sklearn.metrics.regression', 'sklearn.metrics',
     'max_error'),
    ('_ranking', 'sklearn.metrics.ranking', 'sklearn.metrics', 'roc_curve'),
    ('_pairwise_fast', 'sklearn.metrics.pairwise_fast', 'sklearn.metrics',
     'np'),
    ('_scorer', 'sklearn.metrics.scorer', 'sklearn.metrics', 'get_scorer'),

    ('_partial_dependence', 'sklearn.inspection.partial_dependence',
     'sklearn.inspection', 'partial_dependence'),

    ('_ball_tree', 'sklearn.neighbors.ball_tree', 'sklearn.neighbors',
     'BallTree'),
    ('_base', 'sklearn.neighbors.base', 'sklearn.neighbors',
     'VALID_METRICS'),
    ('_classification', 'sklearn.neighbors.classification',
     'sklearn.neighbors', 'KNeighborsClassifier'),
    ('_dist_metrics', 'sklearn.neighbors.dist_metrics', 'sklearn.neighbors',
     'DistanceMetric'),
    ('_graph', 'sklearn.neighbors.graph', 'sklearn.neighbors',
     'KNeighborsTransformer'),
    ('_kd_tree', 'sklearn.neighbors.kd_tree', 'sklearn.neighbors',
     'KDTree'),
    ('_kde', 'sklearn.neighbors.kde', 'sklearn.neighbors',
     'KernelDensity'),
    ('_lof', 'sklearn.neighbors.lof', 'sklearn.neighbors',
     'LocalOutlierFactor'),
    ('_nca', 'sklearn.neighbors.nca', 'sklearn.neighbors',
     'NeighborhoodComponentsAnalysis'),
    ('_nearest_centroid', 'sklearn.neighbors.nearest_centroid',
     'sklearn.neighbors', 'NearestCentroid'),
    ('_quad_tree', 'sklearn.neighbors.quad_tree', 'sklearn.neighbors',
     'CELL_DTYPE'),
    ('_regression', 'sklearn.neighbors.regression', 'sklearn.neighbors',
     'KNeighborsRegressor'),
    ('_typedefs', 'sklearn.neighbors.typedefs', 'sklearn.neighbors',
     'DTYPE'),
    ('_unsupervised', 'sklearn.neighbors.unsupervised', 'sklearn.neighbors',
     'NearestNeighbors'),

    ('_isomap', 'sklearn.manifold.isomap', 'sklearn.manifold', 'Isomap'),
    ('_locally_linear', 'sklearn.manifold.locally_linear', 'sklearn.manifold',
     'LocallyLinearEmbedding'),
    ('_mds', 'sklearn.manifold.mds', 'sklearn.manifold', 'MDS'),
    ('_spectral_embedding', 'sklearn.manifold.spectral_embedding_',
     'sklearn.manifold', 'SpectralEmbedding'),
    ('_t_sne', 'sklearn.manifold.t_sne', 'sklearn.manifold', 'TSNE'),

    ('_label_propagation', 'sklearn.semi_supervised.label_propagation',
     'sklearn.semi_supervised', 'LabelPropagation'),

    ('_data', 'sklearn.preprocessing.data', 'sklearn.preprocessing',
     'Binarizer'),
    ('_label', 'sklearn.preprocessing.label', 'sklearn.preprocessing',
     'LabelEncoder'),

    ('_base', 'sklearn.feature_selection.base', 'sklearn.feature_selection',
     'SelectorMixin'),
    ('_from_model', 'sklearn.feature_selection.from_model',
     'sklearn.feature_selection', 'SelectFromModel'),
    ('_mutual_info', 'sklearn.feature_selection.mutual_info',
     'sklearn.feature_selection', 'mutual_info_regression'),
    ('_rfe', 'sklearn.feature_selection.rfe',
     'sklearn.feature_selection.rfe', 'RFE'),
    ('_univariate_selection',
     'sklearn.feature_selection.univariate_selection',
     'sklearn.feature_selection', 'chi2'),
    ('_variance_threshold',
     'sklearn.feature_selection.variance_threshold',
     'sklearn.feature_selection', 'VarianceThreshold'),

    ('_testing', 'sklearn.utils.testing', 'sklearn.utils',
     'all_estimators'),
]


_FILE_CONTENT_TEMPLATE = """
# THIS FILE WAS AUTOMATICALLY GENERATED BY deprecated_modules.py
import sys
# mypy error: Module X has no attribute y (typically for C extensions)
from . import {new_module_name}  # type: ignore
from {relative_dots}externals._pep562 import Pep562
from {relative_dots}utils.deprecation import _raise_dep_warning_if_not_pytest

deprecated_path = '{deprecated_path}'
correct_import_path = '{correct_import_path}'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_import_path)

def __getattr__(name):
    return getattr({new_module_name}, name)

if not sys.version_info >= (3, 7):
    Pep562(__name__)
"""


def _get_deprecated_path(deprecated_path):
    deprecated_parts = deprecated_path.split(".")
    deprecated_parts[-1] = deprecated_parts[-1] + ".py"
    return Path(*deprecated_parts)


def _create_deprecated_modules_files():
    """Add submodules that will be deprecated. A file is created based
    on the deprecated submodule's name. When this submodule is imported a
    deprecation warning will be raised.
    """
    for (new_module_name, deprecated_path,
         correct_import_path, _) in _DEPRECATED_MODULES:
        relative_dots = deprecated_path.count(".") * "."
        deprecated_content = _FILE_CONTENT_TEMPLATE.format(
            new_module_name=new_module_name,
            relative_dots=relative_dots,
            deprecated_path=deprecated_path,
            correct_import_path=correct_import_path)

        with _get_deprecated_path(deprecated_path).open('w') as f:
            f.write(deprecated_content)


def _clean_deprecated_modules_files():
    """Removes submodules created by _create_deprecated_modules_files."""
    for _, deprecated_path, _, _ in _DEPRECATED_MODULES:
        with suppress(FileNotFoundError):
            _get_deprecated_path(deprecated_path).unlink()


if __name__ == "__main__":
    _clean_deprecated_modules_files()
