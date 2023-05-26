"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""

from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_classification": [
            "accuracy_score",
            "balanced_accuracy_score",
            "class_likelihood_ratios",
            "classification_report",
            "cohen_kappa_score",
            "confusion_matrix",
            "f1_score",
            "fbeta_score",
            "hamming_loss",
            "hinge_loss",
            "jaccard_score",
            "log_loss",
            "matthews_corrcoef",
            "multilabel_confusion_matrix",
            "precision_recall_fscore_support",
            "precision_score",
            "recall_score",
            "zero_one_loss",
            "brier_score_loss",
        ],
        "_dist_metrics": ["DistanceMetric"],
        "_plot.confusion_matrix": ["ConfusionMatrixDisplay"],
        "_plot.det_curve": ["DetCurveDisplay"],
        "_plot.precision_recall_curve": ["PrecisionRecallDisplay"],
        "_plot.roc_curve": ["RocCurveDisplay"],
        "_plot.regression": ["PredictionErrorDisplay"],
        "_ranking": [
            "auc",
            "average_precision_score",
            "coverage_error",
            "dcg_score",
            "det_curve",
            "label_ranking_average_precision_score",
            "label_ranking_loss",
            "ndcg_score",
            "precision_recall_curve",
            "roc_auc_score",
            "roc_curve",
            "top_k_accuracy_score",
        ],
        "_regression": [
            "d2_tweedie_score",
            "d2_absolute_error_score",
            "d2_pinball_score",
            "explained_variance_score",
            "max_error",
            "mean_absolute_error",
            "mean_squared_error",
            "mean_squared_log_error",
            "mean_pinball_loss",
            "mean_poisson_deviance",
            "mean_gamma_deviance",
            "mean_tweedie_deviance",
            "median_absolute_error",
            "mean_absolute_percentage_error",
            "r2_score",
        ],
        "_scorer": [
            "check_scoring",
            "get_scorer",
            "make_scorer",
            "get_scorer_names",
        ],
        "cluster._bicluster": ["consensus_score"],
        "cluster._supervised": [
            "adjusted_mutual_info_score",
            "adjusted_rand_score",
            "completeness_score",
            "fowlkes_mallows_score",
            "homogeneity_completeness_v_measure",
            "homogeneity_score",
            "mutual_info_score",
            "normalized_mutual_info_score",
            "pair_confusion_matrix",
            "rand_score",
            "v_measure_score",
        ],
        "cluster._unsupervised": [
            "calinski_harabasz_score",
            "davies_bouldin_score",
            "silhouette_samples",
            "silhouette_score",
        ],
        "pairwise": [
            "euclidean_distances",
            "nan_euclidean_distances",
            "pairwise_distances",
            "pairwise_distances_argmin",
            "pairwise_distances_argmin_min",
            "pairwise_distances_chunked",
            "pairwise_kernels",
        ],
    },
)
