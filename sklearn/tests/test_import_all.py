"""
Smoke tests for various import statements.

Such an import can fail if the content of the __all__ variable is inconsistent
with the package structure.

"""
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import pkgutil
import sklearn


def test_import_all_consistency():
    # Smoke test to check that any name in a __all__ list is actually defined
    # in the namespace of the module or package.
    for importer, modname, ispkg in pkgutil.walk_packages(
        path=sklearn.__path__, prefix='sklearn.', onerror=lambda x: None):
        if ".tests." in modname or not ispkg:
            continue
        package = __import__(modname, fromlist="dummy")
        for name in getattr(package, '__all__', ()):
            if getattr(package, name, None) is None:
                raise AttributeError(
                    "Module '{}' has no attribute '{}'".format(
                        modname, name))
