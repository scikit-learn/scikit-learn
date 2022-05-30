import re
from inspect import signature
import pkgutil
import inspect
import importlib
from typing import Optional

import pytest
from sklearn.utils import all_estimators
import sklearn

numpydoc_validation = pytest.importorskip("numpydoc.validate")

FUNCTION_DOCSTRING_IGNORE_LIST = [
    "sklearn.datasets._kddcup99.fetch_kddcup99",
    "sklearn.datasets._lfw.fetch_lfw_pairs",
    "sklearn.datasets._lfw.fetch_lfw_people",
    "sklearn.datasets._samples_generator.make_gaussian_quantiles",
    "sklearn.datasets._samples_generator.make_spd_matrix",
    "sklearn.datasets._species_distributions.fetch_species_distributions",
    "sklearn.datasets._svmlight_format_io.load_svmlight_file",
    "sklearn.datasets._svmlight_format_io.load_svmlight_files",
    "sklearn.decomposition._dict_learning.dict_learning",
    "sklearn.decomposition._dict_learning.dict_learning_online",
    "sklearn.decomposition._nmf.non_negative_factorization",
    "sklearn.externals._packaging.version.parse",
    "sklearn.feature_extraction.image.extract_patches_2d",
    "sklearn.feature_extraction.text.strip_accents_unicode",
    "sklearn.feature_selection._univariate_selection.chi2",
    "sklearn.feature_selection._univariate_selection.f_oneway",
    "sklearn.inspection._partial_dependence.partial_dependence",
    "sklearn.inspection._plot.partial_dependence.plot_partial_dependence",
    "sklearn.linear_model._least_angle.lars_path_gram",
    "sklearn.linear_model._omp.orthogonal_mp_gram",
    "sklearn.manifold._locally_linear.locally_linear_embedding",
    "sklearn.manifold._t_sne.trustworthiness",
    "sklearn.metrics._classification.brier_score_loss",
    "sklearn.metrics._classification.cohen_kappa_score",
    "sklearn.metrics._classification.fbeta_score",
    "sklearn.metrics._classification.jaccard_score",
    "sklearn.metrics._classification.log_loss",
    "sklearn.metrics._plot.det_curve.plot_det_curve",
    "sklearn.metrics._plot.precision_recall_curve.plot_precision_recall_curve",
    "sklearn.metrics._ranking.coverage_error",
    "sklearn.metrics._ranking.dcg_score",
    "sklearn.metrics._ranking.label_ranking_average_precision_score",
    "sklearn.metrics._ranking.roc_auc_score",
    "sklearn.metrics._ranking.roc_curve",
    "sklearn.metrics._ranking.top_k_accuracy_score",
    "sklearn.metrics._regression.mean_pinball_loss",
    "sklearn.metrics.cluster._bicluster.consensus_score",
    "sklearn.metrics.cluster._supervised.adjusted_mutual_info_score",
    "sklearn.metrics.cluster._supervised.adjusted_rand_score",
    "sklearn.metrics.cluster._supervised.entropy",
    "sklearn.metrics.cluster._supervised.fowlkes_mallows_score",
    "sklearn.metrics.cluster._supervised.homogeneity_completeness_v_measure",
    "sklearn.metrics.cluster._supervised.mutual_info_score",
    "sklearn.metrics.cluster._supervised.normalized_mutual_info_score",
    "sklearn.metrics.cluster._supervised.pair_confusion_matrix",
    "sklearn.metrics.cluster._supervised.rand_score",
    "sklearn.metrics.cluster._supervised.v_measure_score",
    "sklearn.metrics.pairwise.additive_chi2_kernel",
    "sklearn.metrics.pairwise.check_paired_arrays",
    "sklearn.metrics.pairwise.check_pairwise_arrays",
    "sklearn.metrics.pairwise.chi2_kernel",
    "sklearn.metrics.pairwise.cosine_distances",
    "sklearn.metrics.pairwise.cosine_similarity",
    "sklearn.metrics.pairwise.distance_metrics",
    "sklearn.metrics.pairwise.kernel_metrics",
    "sklearn.metrics.pairwise.paired_manhattan_distances",
    "sklearn.metrics.pairwise.pairwise_distances_argmin",
    "sklearn.metrics.pairwise.pairwise_distances_argmin_min",
    "sklearn.metrics.pairwise.pairwise_distances_chunked",
    "sklearn.metrics.pairwise.pairwise_kernels",
    "sklearn.metrics.pairwise.polynomial_kernel",
    "sklearn.metrics.pairwise.rbf_kernel",
    "sklearn.metrics.pairwise.sigmoid_kernel",
    "sklearn.model_selection._validation.learning_curve",
    "sklearn.model_selection._validation.permutation_test_score",
    "sklearn.model_selection._validation.validation_curve",
    "sklearn.pipeline.make_union",
    "sklearn.preprocessing._data.maxabs_scale",
    "sklearn.preprocessing._data.robust_scale",
    "sklearn.preprocessing._data.scale",
    "sklearn.preprocessing._label.label_binarize",
    "sklearn.random_projection.johnson_lindenstrauss_min_dim",
    "sklearn.svm._bounds.l1_min_c",
    "sklearn.tree._export.plot_tree",
    "sklearn.utils.axis0_safe_slice",
    "sklearn.utils.extmath.density",
    "sklearn.utils.extmath.fast_logdet",
    "sklearn.utils.extmath.randomized_svd",
    "sklearn.utils.extmath.safe_sparse_dot",
    "sklearn.utils.extmath.squared_norm",
    "sklearn.utils.extmath.stable_cumsum",
    "sklearn.utils.extmath.svd_flip",
    "sklearn.utils.extmath.weighted_mode",
    "sklearn.utils.fixes.delayed",
    "sklearn.utils.fixes.linspace",
    # To be fixed in upstream issue:
    # https://github.com/joblib/threadpoolctl/issues/108
    "sklearn.utils.fixes.threadpool_info",
    "sklearn.utils.fixes.threadpool_limits",
    "sklearn.utils.gen_batches",
    "sklearn.utils.gen_even_slices",
    "sklearn.utils.graph.graph_shortest_path",
    "sklearn.utils.graph.single_source_shortest_path_length",
    "sklearn.utils.is_scalar_nan",
    "sklearn.utils.metaestimators.available_if",
    "sklearn.utils.metaestimators.if_delegate_has_method",
    "sklearn.utils.multiclass.class_distribution",
    "sklearn.utils.multiclass.type_of_target",
    "sklearn.utils.multiclass.unique_labels",
    "sklearn.utils.resample",
    "sklearn.utils.safe_mask",
    "sklearn.utils.safe_sqr",
    "sklearn.utils.shuffle",
    "sklearn.utils.sparsefuncs.count_nonzero",
    "sklearn.utils.sparsefuncs.csc_median_axis_0",
    "sklearn.utils.sparsefuncs.incr_mean_variance_axis",
    "sklearn.utils.sparsefuncs.inplace_swap_column",
    "sklearn.utils.sparsefuncs.inplace_swap_row",
    "sklearn.utils.sparsefuncs.inplace_swap_row_csc",
    "sklearn.utils.sparsefuncs.inplace_swap_row_csr",
    "sklearn.utils.sparsefuncs.mean_variance_axis",
    "sklearn.utils.validation.check_is_fitted",
]
FUNCTION_DOCSTRING_IGNORE_LIST = set(FUNCTION_DOCSTRING_IGNORE_LIST)


def get_all_methods():
    estimators = all_estimators()
    for name, Estimator in estimators:
        if name.startswith("_"):
            # skip private classes
            continue
        methods = []
        for name in dir(Estimator):
            if name.startswith("_"):
                continue
            method_obj = getattr(Estimator, name)
            if hasattr(method_obj, "__call__") or isinstance(method_obj, property):
                methods.append(name)
        methods.append(None)

        for method in sorted(methods, key=str):
            yield Estimator, method


def _is_checked_function(item):
    if not inspect.isfunction(item):
        return False

    if item.__name__.startswith("_"):
        return False

    mod = item.__module__
    if not mod.startswith("sklearn.") or mod.endswith("estimator_checks"):
        return False

    return True


def get_all_functions_names():
    """Get all public functions define in the sklearn module"""
    modules_to_ignore = {
        "tests",
        "externals",
        "setup",
        "conftest",
        "experimental",
        "estimator_checks",
    }

    all_functions_names = set()
    for module_finder, module_name, ispkg in pkgutil.walk_packages(
        path=sklearn.__path__, prefix="sklearn."
    ):
        module_parts = module_name.split(".")
        if (
            any(part in modules_to_ignore for part in module_parts)
            or "._" in module_name
        ):
            continue

        module = importlib.import_module(module_name)
        functions = inspect.getmembers(module, _is_checked_function)
        for name, func in functions:
            full_name = f"{func.__module__}.{func.__name__}"
            all_functions_names.add(full_name)

    return sorted(all_functions_names)


def filter_errors(errors, method, Estimator=None):
    """
    Ignore some errors based on the method type.

    These rules are specific for scikit-learn."""
    for code, message in errors:
        # We ignore following error code,
        #  - RT02: The first line of the Returns section
        #    should contain only the type, ..
        #   (as we may need refer to the name of the returned
        #    object)
        #  - GL01: Docstring text (summary) should start in the line
        #    immediately after the opening quotes (not in the same line,
        #    or leaving a blank line in between)
        #  - GL02: If there's a blank line, it should be before the
        #    first line of the Returns section, not after (it allows to have
        #    short docstrings for properties).

        if code in ["RT02", "GL01", "GL02"]:
            continue

        # Ignore PR02: Unknown parameters for properties. We sometimes use
        # properties for ducktyping, i.e. SGDClassifier.predict_proba
        if code == "PR02" and Estimator is not None and method is not None:
            method_obj = getattr(Estimator, method)
            if isinstance(method_obj, property):
                continue

        # Following codes are only taken into account for the
        # top level class docstrings:
        #  - ES01: No extended summary found
        #  - SA01: See Also section not found
        #  - EX01: No examples section found

        if method is not None and code in ["EX01", "SA01", "ES01"]:
            continue
        yield code, message


def repr_errors(res, estimator=None, method: Optional[str] = None) -> str:
    """Pretty print original docstring and the obtained errors

    Parameters
    ----------
    res : dict
        result of numpydoc.validate.validate
    estimator : {estimator, None}
        estimator object or None
    method : str
        if estimator is not None, either the method name or None.

    Returns
    -------
    str
       String representation of the error.
    """
    if method is None:
        if hasattr(estimator, "__init__"):
            method = "__init__"
        elif estimator is None:
            raise ValueError("At least one of estimator, method should be provided")
        else:
            raise NotImplementedError

    if estimator is not None:
        obj = getattr(estimator, method)
        try:
            obj_signature = str(signature(obj))
        except TypeError:
            # In particular we can't parse the signature of properties
            obj_signature = (
                "\nParsing of the method signature failed, "
                "possibly because this is a property."
            )

        obj_name = estimator.__name__ + "." + method
    else:
        obj_signature = ""
        obj_name = method

    msg = "\n\n" + "\n\n".join(
        [
            str(res["file"]),
            obj_name + obj_signature,
            res["docstring"],
            "# Errors",
            "\n".join(
                " - {}: {}".format(code, message) for code, message in res["errors"]
            ),
        ]
    )
    return msg


@pytest.mark.parametrize("function_name", get_all_functions_names())
def test_function_docstring(function_name, request):
    """Check function docstrings using numpydoc."""
    if function_name in FUNCTION_DOCSTRING_IGNORE_LIST:
        request.applymarker(
            pytest.mark.xfail(run=False, reason="TODO pass numpydoc validation")
        )

    res = numpydoc_validation.validate(function_name)

    res["errors"] = list(filter_errors(res["errors"], method="function"))

    if res["errors"]:
        msg = repr_errors(res, method=f"Tested function: {function_name}")

        raise ValueError(msg)


@pytest.mark.parametrize("Estimator, method", get_all_methods())
def test_docstring(Estimator, method, request):
    base_import_path = Estimator.__module__
    import_path = [base_import_path, Estimator.__name__]
    if method is not None:
        import_path.append(method)

    import_path = ".".join(import_path)

    res = numpydoc_validation.validate(import_path)

    res["errors"] = list(filter_errors(res["errors"], method, Estimator=Estimator))

    if res["errors"]:
        msg = repr_errors(res, Estimator, method)

        raise ValueError(msg)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description="Validate docstring with numpydoc.")
    parser.add_argument("import_path", help="Import path to validate")

    args = parser.parse_args()

    res = numpydoc_validation.validate(args.import_path)

    import_path_sections = args.import_path.split(".")
    # When applied to classes, detect class method. For functions
    # method = None.
    # TODO: this detection can be improved. Currently we assume that we have
    # class # methods if the second path element before last is in camel case.
    if len(import_path_sections) >= 2 and re.match(
        r"(?:[A-Z][a-z]*)+", import_path_sections[-2]
    ):
        method = import_path_sections[-1]
    else:
        method = None

    res["errors"] = list(filter_errors(res["errors"], method))

    if res["errors"]:
        msg = repr_errors(res, method=args.import_path)

        print(msg)
        sys.exit(1)
    else:
        print("All docstring checks passed for {}!".format(args.import_path))


def get_all_functions():
    """Get all public functions defined in the sklearn module"""
    modules_to_ignore = {
        "tests",
        "externals",
        "setup",
        "conftest",
        "experimental",
        "estimator_checks",
    }

    all_functions = set()
    for module_finder, module_name, ispkg in pkgutil.walk_packages(
        path=sklearn.__path__, prefix="sklearn."
    ):
        module_parts = module_name.split(".")
        if (
            any(part in modules_to_ignore for part in module_parts)
            or "._" in module_name
        ):
            continue

        module = importlib.import_module(module_name)
        functions = inspect.getmembers(module, _is_checked_function)
        for name, func in functions:
            all_functions.add((func.__module__, name))

    return sorted(all_functions)


FUNCTION_PARAM_VALIDATION_IGNORE_LIST = []


@pytest.mark.parametrize("module, function_name", get_all_functions())
def test_function_param_validation(module, function_name, request):
    """TODO"""
    if function_name in FUNCTION_PARAM_VALIDATION_IGNORE_LIST:
        pytest.skip("TODO use validate_params.")

    return



   base.clone
   base.is_classifier
   base.is_regressor
   config_context
   get_config
   set_config
   show_versions
   calibration.calibration_curve
   cluster.affinity_propagation
   cluster.cluster_optics_dbscan
   cluster.cluster_optics_xi
   cluster.compute_optics_graph
   cluster.dbscan
   cluster.estimate_bandwidth
   cluster.k_means
   cluster.kmeans_plusplus
   cluster.mean_shift
   cluster.spectral_clustering
   cluster.ward_tree
   compose.make_column_transformer
   compose.make_column_selector
   covariance.empirical_covariance
   covariance.graphical_lasso
   covariance.ledoit_wolf
   covariance.oas
   covariance.shrunk_covariance
   datasets.clear_data_home
   datasets.dump_svmlight_file
   datasets.fetch_20newsgroups
   datasets.fetch_20newsgroups_vectorized
   datasets.fetch_california_housing
   datasets.fetch_covtype
   datasets.fetch_kddcup99
   datasets.fetch_lfw_pairs
   datasets.fetch_lfw_people
   datasets.fetch_olivetti_faces
   datasets.fetch_openml
   datasets.fetch_rcv1
   datasets.fetch_species_distributions
   datasets.get_data_home
   datasets.load_boston
   datasets.load_breast_cancer
   datasets.load_diabetes
   datasets.load_digits
   datasets.load_files
   datasets.load_iris
   datasets.load_linnerud
   datasets.load_sample_image
   datasets.load_sample_images
   datasets.load_svmlight_file
   datasets.load_svmlight_files
   datasets.load_wine
   datasets.make_biclusters
   datasets.make_blobs
   datasets.make_checkerboard
   datasets.make_circles
   datasets.make_classification
   datasets.make_friedman1
   datasets.make_friedman2
   datasets.make_friedman3
   datasets.make_gaussian_quantiles
   datasets.make_hastie_10_2
   datasets.make_low_rank_matrix
   datasets.make_moons
   datasets.make_multilabel_classification
   datasets.make_regression
   datasets.make_s_curve
   datasets.make_sparse_coded_signal
   datasets.make_sparse_spd_matrix
   datasets.make_sparse_uncorrelated
   datasets.make_spd_matrix
   datasets.make_swiss_roll
   decomposition.dict_learning
   decomposition.dict_learning_online
   decomposition.fastica
   decomposition.non_negative_factorization
   decomposition.sparse_encode
   feature_extraction.image.extract_patches_2d
   feature_extraction.image.grid_to_graph
   feature_extraction.image.img_to_graph
   feature_extraction.image.reconstruct_from_patches_2d
   feature_selection.chi2
   feature_selection.f_classif
   feature_selection.f_regression
   feature_selection.r_regression
   feature_selection.mutual_info_classif
   feature_selection.mutual_info_regression
   inspection.partial_dependence
   inspection.permutation_importance
   inspection.plot_partial_dependence
   isotonic.check_increasing
   isotonic.isotonic_regression
   linear_model.enet_path
   linear_model.lars_path
   linear_model.lars_path_gram
   linear_model.lasso_path
   linear_model.orthogonal_mp
   linear_model.orthogonal_mp_gram
   linear_model.ridge_regression
   manifold.locally_linear_embedding
   manifold.smacof
   manifold.spectral_embedding
   manifold.trustworthiness
   metrics.check_scoring
   metrics.get_scorer
   metrics.get_scorer_names
   metrics.make_scorer
   metrics.accuracy_score
   metrics.auc
   metrics.average_precision_score
   metrics.balanced_accuracy_score
   metrics.brier_score_loss
   metrics.class_likelihood_ratios
   metrics.classification_report
   metrics.cohen_kappa_score
   metrics.confusion_matrix
   metrics.dcg_score
   metrics.det_curve
   metrics.f1_score
   metrics.fbeta_score
   metrics.hamming_loss
   metrics.hinge_loss
   metrics.jaccard_score
   metrics.log_loss
   metrics.matthews_corrcoef
   metrics.multilabel_confusion_matrix
   metrics.ndcg_score
   metrics.precision_recall_curve
   metrics.precision_recall_fscore_support
   metrics.precision_score
   metrics.recall_score
   metrics.roc_auc_score
   metrics.roc_curve
   metrics.top_k_accuracy_score
   metrics.zero_one_loss
   metrics.explained_variance_score
   metrics.max_error
   metrics.mean_absolute_error
   metrics.mean_squared_error
   metrics.mean_squared_log_error
   metrics.median_absolute_error
   metrics.mean_absolute_percentage_error
   metrics.r2_score
   metrics.mean_poisson_deviance
   metrics.mean_gamma_deviance
   metrics.mean_tweedie_deviance
   metrics.d2_tweedie_score
   metrics.mean_pinball_loss
   metrics.d2_pinball_score
   metrics.d2_absolute_error_score
   metrics.coverage_error
   metrics.label_ranking_average_precision_score
   metrics.label_ranking_loss
   metrics.adjusted_mutual_info_score
   metrics.adjusted_rand_score
   metrics.calinski_harabasz_score
   metrics.davies_bouldin_score
   metrics.completeness_score
   metrics.cluster.contingency_matrix
   metrics.cluster.pair_confusion_matrix
   metrics.fowlkes_mallows_score
   metrics.homogeneity_completeness_v_measure
   metrics.homogeneity_score
   metrics.mutual_info_score
   metrics.normalized_mutual_info_score
   metrics.rand_score
   metrics.silhouette_score
   metrics.silhouette_samples
   metrics.v_measure_score
   metrics.consensus_score
   metrics.pairwise.additive_chi2_kernel
   metrics.pairwise.chi2_kernel
   metrics.pairwise.cosine_similarity
   metrics.pairwise.cosine_distances
   metrics.pairwise.distance_metrics
   metrics.pairwise.euclidean_distances
   metrics.pairwise.haversine_distances
   metrics.pairwise.kernel_metrics
   metrics.pairwise.laplacian_kernel
   metrics.pairwise.linear_kernel
   metrics.pairwise.manhattan_distances
   metrics.pairwise.nan_euclidean_distances
   metrics.pairwise.pairwise_kernels
   metrics.pairwise.polynomial_kernel
   metrics.pairwise.rbf_kernel
   metrics.pairwise.sigmoid_kernel
   metrics.pairwise.paired_euclidean_distances
   metrics.pairwise.paired_manhattan_distances
   metrics.pairwise.paired_cosine_distances
   metrics.pairwise.paired_distances
   metrics.pairwise_distances
   metrics.pairwise_distances_argmin
   metrics.pairwise_distances_argmin_min
   metrics.pairwise_distances_chunked
   metrics.plot_confusion_matrix
   metrics.plot_det_curve
   metrics.plot_precision_recall_curve
   metrics.plot_roc_curve
   model_selection.check_cv
   model_selection.train_test_split
   model_selection.cross_validate
   model_selection.cross_val_predict
   model_selection.cross_val_score
   model_selection.learning_curve
   model_selection.permutation_test_score
   model_selection.validation_curve
   neighbors.kneighbors_graph
   neighbors.radius_neighbors_graph
   neighbors.sort_graph_by_row_values
   pipeline.make_pipeline
   pipeline.make_union
   preprocessing.add_dummy_feature
   preprocessing.binarize
   preprocessing.label_binarize
   preprocessing.maxabs_scale
   preprocessing.minmax_scale
   preprocessing.normalize
   preprocessing.quantile_transform
   preprocessing.robust_scale
   preprocessing.scale
   preprocessing.power_transform
   random_projection.johnson_lindenstrauss_min_dim
   svm.l1_min_c
   tree.export_graphviz
   tree.export_text
   tree.plot_tree
   utils.arrayfuncs.min_pos
   utils.as_float_array
   utils.assert_all_finite
   utils.check_X_y
   utils.check_array
   utils.check_scalar
   utils.check_consistent_length
   utils.check_random_state
   utils.class_weight.compute_class_weight
   utils.class_weight.compute_sample_weight
   utils.deprecated
   utils.estimator_checks.check_estimator
   utils.estimator_checks.parametrize_with_checks
   utils.estimator_html_repr
   utils.extmath.safe_sparse_dot
   utils.extmath.randomized_range_finder
   utils.extmath.randomized_svd
   utils.extmath.fast_logdet
   utils.extmath.density
   utils.extmath.weighted_mode
   utils.gen_batches
   utils.gen_even_slices
   utils.graph.single_source_shortest_path_length
   utils.indexable
   utils.metaestimators.available_if
   utils.multiclass.type_of_target
   utils.multiclass.is_multilabel
   utils.multiclass.unique_labels
   utils.murmurhash3_32
   utils.resample
   utils._safe_indexing
   utils.safe_mask
   utils.safe_sqr
   utils.shuffle
   utils.sparsefuncs.incr_mean_variance_axis
   utils.sparsefuncs.inplace_column_scale
   utils.sparsefuncs.inplace_row_scale
   utils.sparsefuncs.inplace_swap_row
   utils.sparsefuncs.inplace_swap_column
   utils.sparsefuncs.mean_variance_axis
   utils.sparsefuncs.inplace_csr_column_scale
   utils.sparsefuncs_fast.inplace_csr_row_normalize_l1
   utils.sparsefuncs_fast.inplace_csr_row_normalize_l2
   utils.random.sample_without_replacement
   utils.validation.check_is_fitted
   utils.validation.check_memory
   utils.validation.check_symmetric
   utils.validation.column_or_1d
   utils.validation.has_fit_parameter
   utils.all_estimators
   utils.parallel_backend
   utils.register_parallel_backend
   utils.metaestimators.if_delegate_has_method
