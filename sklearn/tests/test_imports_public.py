import importlib
from collections import defaultdict
from importlib import import_module
from pathlib import Path

import pytest

import sklearn
from sklearn.experimental import (  # noqa
    enable_halving_search_cv,
    enable_iterative_imputer,
)

# These functions or classes are added here if:
# - A submodule  (This is okay)
# - Reimported from somewhere else, and is documented there
# - Importable but not documented
# - deprecated
OBJECTS_NOT_IN_API_REFERENCE = {
    "sklearn": [
        "clone",  # reimport: base
    ]
    + sklearn._submodules,
    "sklearn.cluster": [
        "get_bin_seeds",  # undocumented
    ],
    "sklearn.covariance": [
        "log_likelihood",  # undocumented
    ],
    "sklearn.feature_extraction": [
        "image",  # submodule
        "text",  # submodule
        "img_to_graph",  # reimport: image
        "grid_to_graph",  # reimport: image
    ],
    "sklearn.feature_extraction.text": [
        "ENGLISH_STOP_WORDS",  # undocumented
        "strip_accents_ascii",  # undocumented
        "strip_accents_unicode",  # undocumented
        "strip_tags",  # undocumented
    ],
    "sklearn.gaussian_process": [
        "kernels",  # submodule
    ],
    "sklearn.metrics": [
        "cluster",  # submodule
        "euclidean_distances",  # reimport: pairwise
        "nan_euclidean_distances",  # reimport: pairwise
        "pair_confusion_matrix",  # undocumented
        "pairwise_kernels",  # reimport: pairwise
    ],
    "sklearn.metrics.cluster": [
        "fowlkes_mallows_score",  # reimport: metrics
        "v_measure_score",  # reimport: metrics
        "entropy",  # undocumented
        "expected_mutual_information",  # undocumented
        "mutual_info_score",  # reimport: metrics
        "normalized_mutual_info_score",  # reimport: metrics
        "calinski_harabasz_score",  # reimport: metrics
        "adjusted_rand_score",  # reimport: metrics
        "homogeneity_score",  # reimport: metrics
        "consensus_score",  # reimport: metrics
        "davies_bouldin_score",  # reimport: metrics
        "silhouette_score",  # reimport: metrics
        "completeness_score",  # reimport: metrics
        "adjusted_mutual_info_score",  # reimport: metrics
        "rand_score",  # reimport: metrics
        "homogeneity_completeness_v_measure",  # reimport: metrics
        "silhouette_samples",  # reimport: metrics
    ],
    "sklearn.model_selection": [
        "BaseCrossValidator",  # undocumented
        "BaseShuffleSplit",  # undocumented
    ],
    "sklearn.neighbors": [
        "VALID_METRICS",  # undocumented
        "VALID_METRICS_SPARSE",  # undocumented
    ],
    "sklearn.tree": [
        "BaseDecisionTree",  # undocumented
    ],
    "sklearn.utils.multiclass": [
        "check_classification_targets",  # undocumented
        "class_distribution",  # undocumented
    ],
    "sklearn.utils.extmath": [
        "make_nonnegative",  # undocumented
        "svd_flip",  # undocumented
        "row_norms",  # undocumented
        "cartesian",  # undocumented
        "softmax",  # undocumented
        "stable_cumsum",  # undocumented
        "squared_norm",  # undocumented
    ],
    "sklearn.utils.validation": [
        "assert_all_finite",  # reimport: utils
    ],
    "sklearn.utils": [
        "metadata_routing",  # module
        "parallel_backend",  # deprecated
        "check_symmetric",  # reimport: utils.validation
        "compute_sample_weight",  # reimport: utils.class_weight
        "default_tags",  # undocumented
        "column_or_1d",  # reimport: utils.validation
        "all_estimators",  # reimport: discovery.all_estimator
        "compute_class_weight",  # reimport: utils.class_weight
        "tosequence",  # undocumented
        "DataConversionWarning",  # reimport exceptions
    ],
}


def import_module_from_path(path):
    module_name = path.name
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def yield_all_public_apis():
    root_directory = Path(sklearn.__file__).parent.parent
    api_reference_file = root_directory / "doc" / "api_reference.py"

    if not api_reference_file.exists():
        return []

    api_reference = import_module_from_path(api_reference_file)
    module_to_public_names = defaultdict(list)
    module_to_public_names.update(OBJECTS_NOT_IN_API_REFERENCE)

    for module_name, info in api_reference.API_REFERENCE.items():
        for section in info["sections"]:
            for public_name in section["autosummary"]:
                if "." in public_name:
                    assert public_name.count(".") == 1
                    submodule, public_name = public_name.split(".")
                    module_to_public_names[f"{module_name}.{submodule}"].append(
                        public_name
                    )
                else:
                    # Part of the parent module
                    module_to_public_names[module_name].append(public_name)

    for mod_name, public_names in module_to_public_names.items():
        yield mod_name, sorted(public_names)


@pytest.mark.parametrize(
    "module_name, public_names",
    yield_all_public_apis(),
)
def test_public_functions_are_in_all_and_dir(module_name, public_names):
    """Check that public functions and in __all__ and returned by __dir__()."""
    module = import_module(module_name)

    if module_name == "sklearn.experimental":
        pytest.skip(reason="Do not need to run sklearn.experimental")

    public_name_set = set(public_names)
    module_all_set = set(module.__all__)

    in_module_but_not_in_api = module_all_set - public_name_set
    in_api_but_not_in_module = public_name_set - module_all_set

    errors = []
    if in_module_but_not_in_api:
        errors.append(
            f"api_reference.py for {module.__name__} is missing "
            f"{in_module_but_not_in_api}"
        )
    if in_api_but_not_in_module:
        errors.append(
            f"{module.__name__}'s __all__ is missing {in_api_but_not_in_module}"
        )

    assert not in_module_but_not_in_api and not in_api_but_not_in_module, "\n".join(
        errors
    )

    assert public_names == sorted(
        module.__dir__()
    ), f"Expected {module.__name__}'s __dir__() to be:\n{public_names}"
