import re
from inspect import signature

import pytest
from sklearn.utils.testing import all_estimators

numpydoc_validation = pytest.importorskip("numpydoc.validate")

# List of whitelisted modules and methods; regexp is supported.
#
DOCSTRING_WHITELIST = [r"LogisticRegression$"]


def get_all_methods():
    estimators = all_estimators(include_meta_estimators=True)
    for name, estimator in estimators:
        if name.startswith("_"):
            # skip private classes
            continue
        methods = [el for el in dir(estimator) if not el.startswith("_")]
        methods.append(None)

        for method in sorted(methods, key=lambda x: str(x)):
            yield estimator, method


def filter_errors(errors, method):
    """Ignore some errors based on the method type"""
    for code, message in errors:
        # Following codes are only taken into account for the
        # top level class docstrings:
        #  - ES01: No extended summary found
        #  - SA01: See Also section not found
        #  - EX01: No examples section found
        if method is not None and code in [
            "EX01",
            "SA01",
            "ES01",
        ]:
            continue
        yield code, message

def repr_errors(res, estimator, method: str) -> str:
    """Pretty print original docstring and the obtained errors"""
    if method is None:
        if hasattr(estimator, "__init__"):
            method = "__init__"
        else:
            raise NotImplementedError

    obj_signature = signature(getattr(estimator, method))
    obj_name = estimator.__name__ + "." + method

    msg = "\n\n" + "\n\n".join(
        [
            res['file'],
            obj_name + str(obj_signature),
            res["docstring"],
            "# Errors",
            "\n".join(
                " - {}: {}".format(code, message)
                for code, message in res['errors']
            ),
        ]
    )
    return msg


@pytest.mark.parametrize("estimator, method", get_all_methods())
def test_docstring(estimator, method):
    base_import_path = estimator.__module__
    import_path = [base_import_path, estimator.__name__]
    if method is not None:
        import_path.append(method)

    import_path = ".".join(import_path)

    if not any(re.search(regex, import_path) for regex in DOCSTRING_WHITELIST):
        pytest.xfail(reason="TODO")

    res = numpydoc_validation.validate(import_path)

    res['errors'] = list(filter_errors(res["errors"], method))

    if res['errors']:
        msg = repr_errors(res, estimator, method)

        raise ValueError(msg)
