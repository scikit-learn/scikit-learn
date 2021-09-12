import re
from inspect import signature
from typing import Optional

import pytest
from sklearn.utils import all_estimators

numpydoc_validation = pytest.importorskip("numpydoc.validate")

# List of modules ignored when checking for numpydoc validation.
DOCSTRING_IGNORE_LIST = [
    "GridSearchCV",
    "HalvingGridSearchCV",
    "HalvingRandomSearchCV",
    "HashingVectorizer",
    "HuberRegressor",
    "Isomap",
    "IterativeImputer",
    "KNNImputer",
    "LabelPropagation",
    "LabelSpreading",
    "LocalOutlierFactor",
    "LocallyLinearEmbedding",
    "MDS",
    "MeanShift",
    "MiniBatchKMeans",
    "MissingIndicator",
    "MultiLabelBinarizer",
    "MultiTaskElasticNet",
    "MultiTaskElasticNetCV",
    "MultiTaskLasso",
    "MultiTaskLassoCV",
    "NearestCentroid",
    "NeighborhoodComponentsAnalysis",
    "Normalizer",
    "OPTICS",
    "OneVsOneClassifier",
    "OneVsRestClassifier",
    "OrdinalEncoder",
    "OrthogonalMatchingPursuit",
    "OrthogonalMatchingPursuitCV",
    "OutputCodeClassifier",
    "PLSCanonical",
    "PLSRegression",
    "PLSSVD",
    "PassiveAggressiveClassifier",
    "PassiveAggressiveRegressor",
    "PatchExtractor",
    "PolynomialFeatures",
    "PowerTransformer",
    "QuadraticDiscriminantAnalysis",
    "QuantileRegressor",
    "QuantileTransformer",
    "RANSACRegressor",
    "RadiusNeighborsTransformer",
    "RandomizedSearchCV",
    "RobustScaler",
    "SGDOneClassSVM",
    "SGDRegressor",
    "SelfTrainingClassifier",
    "SequentialFeatureSelector",
    "SimpleImputer",
    "SparseRandomProjection",
    "SpectralBiclustering",
    "SpectralClustering",
    "SpectralCoclustering",
    "SpectralEmbedding",
    "SplineTransformer",
    "StackingClassifier",
    "StackingRegressor",
    "TheilSenRegressor",
    "TransformedTargetRegressor",
    "TweedieRegressor",
    "VarianceThreshold",
]


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

        for method in sorted(methods, key=lambda x: str(x)):
            yield Estimator, method


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
            obj_signature = signature(obj)
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
            obj_name + str(obj_signature),
            res["docstring"],
            "# Errors",
            "\n".join(
                " - {}: {}".format(code, message) for code, message in res["errors"]
            ),
        ]
    )
    return msg


@pytest.mark.parametrize("Estimator, method", get_all_methods())
def test_docstring(Estimator, method, request):
    base_import_path = Estimator.__module__
    import_path = [base_import_path, Estimator.__name__]
    if method is not None:
        import_path.append(method)

    import_path = ".".join(import_path)

    if Estimator.__name__ in DOCSTRING_IGNORE_LIST:
        request.applymarker(
            pytest.mark.xfail(run=False, reason="TODO pass numpydoc validation")
        )

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
