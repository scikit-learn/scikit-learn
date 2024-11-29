import importlib
from collections import defaultdict
from importlib import import_module
from pathlib import Path

import pytest

import sklearn
from sklearn.experimental import (
    enable_halving_search_cv,  # noqa
    enable_iterative_imputer,  # noqa
)

# These functions or classes are added here if:
# - Reimported from another module
# - Importable but not documented
# - A submodule
EXPORTED_REIMPORTS = {
    "sklearn": ["clone"] + sklearn._submodules,
    "sklearn.feature_extraction": [
        "image",
        "text",
        "img_to_graph",
        "grid_to_graph",
    ],
    "sklearn.gaussian_process": [
        "kernels",
    ],
    "sklearn.metrics": [
        "cluster",
        "euclidean_distances",
        "nan_euclidean_distances",
        "pair_confusion_matrix",
        "pairwise_kernels",
    ],
    "sklearn.metrics.cluster": [
        "fowlkes_mallows_score",
        "v_measure_score",
        "entropy",
        "expected_mutual_information",
        "mutual_info_score",
        "normalized_mutual_info_score",
        "calinski_harabasz_score",
        "adjusted_rand_score",
        "homogeneity_score",
        "consensus_score",
        "davies_bouldin_score",
        "silhouette_score",
        "completeness_score",
        "adjusted_mutual_info_score",
        "rand_score",
        "homogeneity_completeness_v_measure",
        "silhouette_samples",
    ],
    "sklearn.tree": [
        "BaseDecisionTree",
    ],
    "sklearn.utils.multiclass": [
        "check_classification_targets",
        "class_distribution",
    ],
    "sklearn.utils.extmath": [
        "make_nonnegative",
        "svd_flip",
        "row_norms",
        "cartesian",
        "softmax",
        "stable_cumsum",
        "squared_norm",
    ],
    "sklearn.utils.validation": ["assert_all_finite"],
    "sklearn.utils": [
        "parallel_backend",
        "check_symmetric",
        "compute_sample_weight",
        "default_tags",
        "column_or_1d",
        "all_estimators",
        "compute_class_weight",
        "tosequence",
        "metadata_routing",
        "DataConversionWarning",
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
    module_to_public_names.update(EXPORTED_REIMPORTS)
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

    # Remove when https://github.com/scikit-learn/scikit-learn/pull/30368 is
    # merged
    if module_name == "sklearn.base":
        public_names.remove("is_transformer")

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
