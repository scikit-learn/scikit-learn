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
    ('_hierarchical', 'sklearn.cluster.hierarchical', 'sklearn.cluster',
     'FeatureAgglomeration'),
    ('_k_means', 'sklearn.cluster.k_means_', 'sklearn.cluster', 'KMeans'),
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

    ('_empirical_covariance_', 'sklearn.covariance.empirical_covariance_',
     'sklearn.covariance', 'EmpiricalCovariance'),
    ('_shrunk_covariance_', 'sklearn.covariance.shrunk_covariance_',
     'sklearn.covariance', 'ShrunkCovariance'),
    ('_robust_covariance', 'sklearn.covariance.robust_covariance',
     'sklearn.covariance', 'MinCovDet'),
    ('_graph_lasso_', 'sklearn.covariance.graph_lasso_',
     'sklearn.covariance', 'GraphicalLasso'),
    ('_elliptic_envelope', 'sklearn.covariance.elliptic_envelope',
     'sklearn.covariance', 'EllipticEnvelope'),

    ('_cca_', 'sklearn.cross_decomposition.cca_',
     'sklearn.cross_decomposition', 'CCA'),
    ('_pls_', 'sklearn.cross_decomposition.pls_',
     'sklearn.cross_decomposition', 'PLSSVD'),

    ('_base', 'sklearn.svm.base', 'sklearn.svm', 'BaseLibSVM'),
    ('_bounds', 'sklearn.svm.bounds', 'sklearn.svm', 'l1_min_c'),
    ('_classes', 'sklearn.svm.classes', 'sklearn.svm', 'SVR'),
    ('_libsvm', 'sklearn.svm.libsvm', 'sklearn.svm', 'fit'),
    ('_libsvm_sparse', 'sklearn.svm.libsvm_sparse', 'sklearn.svm',
     'set_verbosity_wrap'),
    ('_liblinear', 'sklearn.svm.liblinear', 'sklearn.svm', 'train_wrap'),

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
    ('_permutation_importance', 'sklearn.inspection.permutation_importance',
     'sklearn.inspection', 'permutation_importance'),

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
    ('_spectral_embedding_', 'sklearn.manifold.spectral_embedding_',
     'sklearn.manifold', 'SpectralEmbedding'),
    ('_t_sne', 'sklearn.manifold.t_sne', 'sklearn.manifold', 'TSNE'),

    ('_label_propagation', 'sklearn.semi_supervised.label_propagation',
     'sklearn.semi_supervised', 'LabelPropagation'),

    ('_data', 'sklearn.preprocessing.data', 'sklearn.preprocessing',
     'Binarizer'),
    ('_label', 'sklearn.preprocessing.label', 'sklearn.preprocessing',
     'LabelEncoder'),
]


_FILE_CONTENT_TEMPLATE = """
# THIS FILE WAS AUTOMATICALLY GENERATED BY deprecated_modules.py

from .{new_module_name} import *  # noqa
from {relative_dots}utils.deprecation import _raise_dep_warning_if_not_pytest

deprecated_path = '{deprecated_path}'
correct_import_path = '{correct_import_path}'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_import_path)
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
