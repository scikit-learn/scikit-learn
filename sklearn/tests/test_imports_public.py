import inspect
import itertools
import pkgutil
from importlib import import_module

import pytest

import sklearn

module_part_to_ignore = {
    "tests",
    "externals",
    "setup",
    "conftest",
    "experimental",
    "estimator_checks",
}

# Empty modules
module_name_to_ignore = {
    "sklearn.datasets.data",
    "sklearn.datasets.descr",
    "sklearn.datasets.images",
}

# Functions and classes from these modules are not public
objects_in_module_to_ignore = {
    "sklearn.utils._param_validation",
    "sklearn.utils._array_api",
    "sklearn.utils._unique",
}

# These public names are imported in an __init__.py file, but should also be imported
# from the a specific module
public_names_to_include_in_module = {
    "sklearn.utils.class_weight": {"compute_class_weight", "compute_sample_weight"},
    "sklearn.utils.discovery": {"all_estimators"},
    "sklearn.utils.extmath": {"randomized_svd", "randomized_range_finder"},
}


def yield_module_and_public_names():
    root_path = sklearn.__path__
    module_with_public_names = []

    for _, module_name, _ in pkgutil.walk_packages(path=root_path, prefix="sklearn."):
        if module_name in module_name_to_ignore:
            continue

        module_parts = module_name.split(".")
        if (
            any(part in module_part_to_ignore for part in module_parts)
            or "._" in module_name
        ):
            continue

        module = import_module(module_name)

        # Only check python files
        if not module.__file__.endswith(".py"):
            continue

        def _is_public_class_function(x):
            return (hasattr(x, "__name__") and not x.__name__.startswith("_")) and (
                hasattr(x, "__module__")
                and x.__module__.startswith("sklearn.")
                and x.__module__ not in objects_in_module_to_ignore
            )

        sklearn_objects = inspect.getmembers(module, _is_public_class_function)
        public_names = set([o[0] for o in sklearn_objects])

        module_with_public_names.append((module, public_names))

    # Public names in __init__.py and sklearn.base are generally not imported by
    # other modules
    public_names_to_remove = set(
        itertools.chain.from_iterable(
            public_names
            for module, public_names in module_with_public_names
            if module.__file__.endswith("__init__.py")
            or module.__name__ == "sklearn.base"
        )
    )

    for module, public_names_in_module in module_with_public_names:
        if module.__file__.endswith("__init__.py"):
            yield module, sorted(public_names_in_module)
        else:
            if module.__name__ in public_names_to_include_in_module:
                public_modules_to_remove = (
                    public_names_to_remove
                    - public_names_to_include_in_module[module.__name__]
                )
            else:
                public_modules_to_remove = public_names_to_remove

            yield module, sorted(public_names_in_module - public_modules_to_remove)


@pytest.mark.parametrize(
    "module, public_names",
    yield_module_and_public_names(),
)
def test_public_functions_consistent_with_all_and_dir(module, public_names):
    assert (
        module.__all__ == public_names
    ), f"Expected {module.__name__} __all__ to have:\n{public_names}"

    assert (
        module.__dir__() == public_names
    ), f"Expected {module.__name__} __dir__() to return:\n{public_names}"
