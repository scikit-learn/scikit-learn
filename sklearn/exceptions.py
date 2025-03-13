"""Custom warnings and errors used across scikit-learn."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

__all__ = [
    "NotFittedError",
    "ConvergenceWarning",
    "DataConversionWarning",
    "DataDimensionalityWarning",
    "EfficiencyWarning",
    "FitFailedWarning",
    "SkipTestWarning",
    "UndefinedMetricWarning",
    "PositiveSpectrumWarning",
    "UnsetMetadataPassedError",
    "EstimatorCheckFailedWarning",
]


class UnsetMetadataPassedError(ValueError):
    """Exception class to raise if a metadata is passed which is not explicitly \
        requested (metadata=True) or not requested (metadata=False).

    .. versionadded:: 1.3

    Parameters
    ----------
    message : str
        The message

    unrequested_params : dict
        A dictionary of parameters and their values which are provided but not
        requested.

    routed_params : dict
        A dictionary of routed parameters.
    """

    def __init__(self, *, message, unrequested_params, routed_params):
        super().__init__(message)
        self.unrequested_params = unrequested_params
        self.routed_params = routed_params


class NotFittedError(ValueError, AttributeError):
    """Exception class to raise if estimator is used before fitting.

    This class inherits from both ValueError and AttributeError to help with
    exception handling and backward compatibility.

    Examples
    --------
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.exceptions import NotFittedError
    >>> try:
    ...     LinearSVC().predict([[1, 2], [2, 3], [3, 4]])
    ... except NotFittedError as e:
    ...     print(repr(e))
    NotFittedError("This LinearSVC instance is not fitted yet. Call 'fit' with
    appropriate arguments before using this estimator."...)

    .. versionchanged:: 0.18
       Moved from sklearn.utils.validation.
    """


class ConvergenceWarning(UserWarning):
    """Custom warning to capture convergence problems

    .. versionchanged:: 0.18
       Moved from sklearn.utils.
    """


class DataConversionWarning(UserWarning):
    """Warning used to notify implicit data conversions happening in the code.

    This warning occurs when some input data needs to be converted or
    interpreted in a way that may not match the user's expectations.

    For example, this warning may occur when the user
        - passes an integer array to a function which expects float input and
          will convert the input
        - requests a non-copying operation, but a copy is required to meet the
          implementation's data-type expectations;
        - passes an input whose shape can be interpreted ambiguously.

    .. versionchanged:: 0.18
       Moved from sklearn.utils.validation.
    """


class DataDimensionalityWarning(UserWarning):
    """Custom warning to notify potential issues with data dimensionality.

    For example, in random projection, this warning is raised when the
    number of components, which quantifies the dimensionality of the target
    projection space, is higher than the number of features, which quantifies
    the dimensionality of the original source space, to imply that the
    dimensionality of the problem will not be reduced.

    .. versionchanged:: 0.18
       Moved from sklearn.utils.
    """


class EfficiencyWarning(UserWarning):
    """Warning used to notify the user of inefficient computation.

    This warning notifies the user that the efficiency may not be optimal due
    to some reason which may be included as a part of the warning message.
    This may be subclassed into a more specific Warning class.

    .. versionadded:: 0.18
    """


class FitFailedWarning(RuntimeWarning):
    """Warning class used if there is an error while fitting the estimator.

    This Warning is used in meta estimators GridSearchCV and RandomizedSearchCV
    and the cross-validation helper function cross_val_score to warn when there
    is an error while fitting the estimator.

    .. versionchanged:: 0.18
       Moved from sklearn.cross_validation.
    """


class SkipTestWarning(UserWarning):
    """Warning class used to notify the user of a test that was skipped.

    For example, one of the estimator checks requires a pandas import.
    If the pandas package cannot be imported, the test will be skipped rather
    than register as a failure.
    """


class UndefinedMetricWarning(UserWarning):
    """Warning used when the metric is invalid

    .. versionchanged:: 0.18
       Moved from sklearn.base.
    """


class PositiveSpectrumWarning(UserWarning):
    """Warning raised when the eigenvalues of a PSD matrix have issues

    This warning is typically raised by ``_check_psd_eigenvalues`` when the
    eigenvalues of a positive semidefinite (PSD) matrix such as a gram matrix
    (kernel) present significant negative eigenvalues, or bad conditioning i.e.
    very small non-zero eigenvalues compared to the largest eigenvalue.

    .. versionadded:: 0.22
    """


class InconsistentVersionWarning(UserWarning):
    """Warning raised when an estimator is unpickled with a inconsistent version.

    Parameters
    ----------
    estimator_name : str
        Estimator name.

    current_sklearn_version : str
        Current scikit-learn version.

    original_sklearn_version : str
        Original scikit-learn version.
    """

    def __init__(
        self, *, estimator_name, current_sklearn_version, original_sklearn_version
    ):
        self.estimator_name = estimator_name
        self.current_sklearn_version = current_sklearn_version
        self.original_sklearn_version = original_sklearn_version

    def __str__(self):
        return (
            f"Trying to unpickle estimator {self.estimator_name} from version"
            f" {self.original_sklearn_version} when "
            f"using version {self.current_sklearn_version}. This might lead to breaking"
            " code or "
            "invalid results. Use at your own risk. "
            "For more info please refer to:\n"
            "https://scikit-learn.org/stable/model_persistence.html"
            "#security-maintainability-limitations"
        )


class EstimatorCheckFailedWarning(UserWarning):
    """Warning raised when an estimator check from the common tests fails.

    Parameters
    ----------
    estimator : estimator object
        Estimator instance for which the test failed.

    check_name : str
        Name of the check that failed.

    exception : Exception
        Exception raised by the failed check.

    status : str
        Status of the check.

    expected_to_fail : bool
        Whether the check was expected to fail.

    expected_to_fail_reason : str
        Reason for the expected failure.
    """

    def __init__(
        self,
        *,
        estimator,
        check_name: str,
        exception: Exception,
        status: str,
        expected_to_fail: bool,
        expected_to_fail_reason: str,
    ):
        self.estimator = estimator
        self.check_name = check_name
        self.exception = exception
        self.status = status
        self.expected_to_fail = expected_to_fail
        self.expected_to_fail_reason = expected_to_fail_reason

    def __repr__(self):
        expected_to_fail_str = (
            f"Expected to fail: {self.expected_to_fail_reason}"
            if self.expected_to_fail
            else "Not expected to fail"
        )
        return (
            f"Test {self.check_name} failed for estimator {self.estimator!r}.\n"
            f"Expected to fail reason: {expected_to_fail_str}\n"
            f"Exception: {self.exception}"
        )

    def __str__(self):
        return self.__repr__()
