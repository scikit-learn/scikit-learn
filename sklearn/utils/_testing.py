"""Testing utilities."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import atexit
import contextlib
import functools
import importlib
import inspect
import os
import os.path as op
import re
import shutil
import sys
import tempfile
import textwrap
import unittest
import warnings
from collections import defaultdict, namedtuple
from collections.abc import Iterable
from dataclasses import dataclass
from difflib import context_diff
from functools import wraps
from inspect import signature
from itertools import chain, groupby
from subprocess import STDOUT, CalledProcessError, TimeoutExpired, check_output

import joblib
import numpy as np
import scipy as sp
from numpy.testing import assert_allclose as np_assert_allclose
from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
    assert_array_less,
)

from sklearn import __file__ as sklearn_path
from sklearn.utils import (
    ClassifierTags,
    RegressorTags,
    Tags,
    TargetTags,
    TransformerTags,
)
from sklearn.utils._array_api import _check_array_api_dispatch
from sklearn.utils.fixes import (
    _IS_32BIT,
    VisibleDeprecationWarning,
    _in_unstable_openblas_configuration,
)
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.validation import check_array, check_is_fitted, check_X_y

__all__ = [
    "SkipTest",
    "assert_allclose",
    "assert_almost_equal",
    "assert_array_almost_equal",
    "assert_array_equal",
    "assert_array_less",
    "assert_run_python_script_without_output",
]

SkipTest = unittest.case.SkipTest


def ignore_warnings(obj=None, category=Warning):
    """Context manager and decorator to ignore warnings.

    Note: Using this (in both variants) will clear all warnings
    from all python modules loaded. In case you need to test
    cross-module-warning-logging, this is not your tool of choice.

    Parameters
    ----------
    obj : callable, default=None
        callable where you want to ignore the warnings.
    category : warning class, default=Warning
        The category to filter. If Warning, all categories will be muted.

    Examples
    --------
    >>> import warnings
    >>> from sklearn.utils._testing import ignore_warnings
    >>> with ignore_warnings():
    ...     warnings.warn('buhuhuhu')

    >>> def nasty_warn():
    ...     warnings.warn('buhuhuhu')
    ...     print(42)

    >>> ignore_warnings(nasty_warn)()
    42
    """
    if isinstance(obj, type) and issubclass(obj, Warning):
        # Avoid common pitfall of passing category as the first positional
        # argument which result in the test not being run
        warning_name = obj.__name__
        raise ValueError(
            "'obj' should be a callable where you want to ignore warnings. "
            "You passed a warning class instead: 'obj={warning_name}'. "
            "If you want to pass a warning class to ignore_warnings, "
            "you should use 'category={warning_name}'".format(warning_name=warning_name)
        )
    elif callable(obj):
        return _IgnoreWarnings(category=category)(obj)
    else:
        return _IgnoreWarnings(category=category)


class _IgnoreWarnings:
    """Improved and simplified Python warnings context manager and decorator.

    This class allows the user to ignore the warnings raised by a function.
    Copied from Python 2.7.5 and modified as required.

    Parameters
    ----------
    category : tuple of warning class, default=Warning
        The category to filter. By default, all the categories will be muted.

    """

    def __init__(self, category):
        self._record = True
        self._module = sys.modules["warnings"]
        self._entered = False
        self.log = []
        self.category = category

    def __call__(self, fn):
        """Decorator to catch and hide warnings without visual nesting."""

        @wraps(fn)
        def wrapper(*args, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", self.category)
                return fn(*args, **kwargs)

        return wrapper

    def __repr__(self):
        args = []
        if self._record:
            args.append("record=True")
        if self._module is not sys.modules["warnings"]:
            args.append("module=%r" % self._module)
        name = type(self).__name__
        return "%s(%s)" % (name, ", ".join(args))

    def __enter__(self):
        if self._entered:
            raise RuntimeError("Cannot enter %r twice" % self)
        self._entered = True
        self._filters = self._module.filters
        self._module.filters = self._filters[:]
        self._showwarning = self._module.showwarning
        warnings.simplefilter("ignore", self.category)

    def __exit__(self, *exc_info):
        if not self._entered:
            raise RuntimeError("Cannot exit %r without entering first" % self)
        self._module.filters = self._filters
        self._module.showwarning = self._showwarning
        self.log[:] = []


def assert_allclose(
    actual, desired, rtol=None, atol=0.0, equal_nan=True, err_msg="", verbose=True
):
    """dtype-aware variant of numpy.testing.assert_allclose

    This variant introspects the least precise floating point dtype
    in the input argument and automatically sets the relative tolerance
    parameter to 1e-4 float32 and use 1e-7 otherwise (typically float64
    in scikit-learn).

    `atol` is always left to 0. by default. It should be adjusted manually
    to an assertion-specific value in case there are null values expected
    in `desired`.

    The aggregate tolerance is `atol + rtol * abs(desired)`.

    Parameters
    ----------
    actual : array_like
        Array obtained.
    desired : array_like
        Array desired.
    rtol : float, optional, default=None
        Relative tolerance.
        If None, it is set based on the provided arrays' dtypes.
    atol : float, optional, default=0.
        Absolute tolerance.
    equal_nan : bool, optional, default=True
        If True, NaNs will compare equal.
    err_msg : str, optional, default=''
        The error message to be printed in case of failure.
    verbose : bool, optional, default=True
        If True, the conflicting values are appended to the error message.

    Raises
    ------
    AssertionError
        If actual and desired are not equal up to specified precision.

    See Also
    --------
    numpy.testing.assert_allclose

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils._testing import assert_allclose
    >>> x = [1e-5, 1e-3, 1e-1]
    >>> y = np.arccos(np.cos(x))
    >>> assert_allclose(x, y, rtol=1e-5, atol=0)
    >>> a = np.full(shape=10, fill_value=1e-5, dtype=np.float32)
    >>> assert_allclose(a, 1e-5)
    """
    dtypes = []

    actual, desired = np.asanyarray(actual), np.asanyarray(desired)
    dtypes = [actual.dtype, desired.dtype]

    if rtol is None:
        rtols = [1e-4 if dtype == np.float32 else 1e-7 for dtype in dtypes]
        rtol = max(rtols)

    np_assert_allclose(
        actual,
        desired,
        rtol=rtol,
        atol=atol,
        equal_nan=equal_nan,
        err_msg=err_msg,
        verbose=verbose,
    )


def assert_allclose_dense_sparse(x, y, rtol=1e-07, atol=1e-9, err_msg=""):
    """Assert allclose for sparse and dense data.

    Both x and y need to be either sparse or dense, they
    can't be mixed.

    Parameters
    ----------
    x : {array-like, sparse matrix}
        First array to compare.

    y : {array-like, sparse matrix}
        Second array to compare.

    rtol : float, default=1e-07
        relative tolerance; see numpy.allclose.

    atol : float, default=1e-9
        absolute tolerance; see numpy.allclose. Note that the default here is
        more tolerant than the default for numpy.testing.assert_allclose, where
        atol=0.

    err_msg : str, default=''
        Error message to raise.
    """
    if sp.sparse.issparse(x) and sp.sparse.issparse(y):
        x = x.tocsr()
        y = y.tocsr()
        x.sum_duplicates()
        y.sum_duplicates()
        assert_array_equal(x.indices, y.indices, err_msg=err_msg)
        assert_array_equal(x.indptr, y.indptr, err_msg=err_msg)
        assert_allclose(x.data, y.data, rtol=rtol, atol=atol, err_msg=err_msg)
    elif not sp.sparse.issparse(x) and not sp.sparse.issparse(y):
        # both dense
        assert_allclose(x, y, rtol=rtol, atol=atol, err_msg=err_msg)
    else:
        raise ValueError(
            "Can only compare two sparse matrices, not a sparse matrix and an array."
        )


def set_random_state(estimator, random_state=0):
    """Set random state of an estimator if it has the `random_state` param.

    Parameters
    ----------
    estimator : object
        The estimator.
    random_state : int, RandomState instance or None, default=0
        Pseudo random number generator state.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.
    """
    if "random_state" in estimator.get_params():
        estimator.set_params(random_state=random_state)


def _is_numpydoc():
    try:
        import numpydoc  # noqa: F401
    except (ImportError, AssertionError):
        return False
    else:
        return True


try:
    _check_array_api_dispatch(True)
    ARRAY_API_COMPAT_FUNCTIONAL = True
except (ImportError, RuntimeError):
    ARRAY_API_COMPAT_FUNCTIONAL = False

try:
    import pytest

    skip_if_32bit = pytest.mark.skipif(_IS_32BIT, reason="skipped on 32bit platforms")
    fails_if_unstable_openblas = pytest.mark.xfail(
        _in_unstable_openblas_configuration(),
        reason="OpenBLAS is unstable for this configuration",
    )
    skip_if_no_parallel = pytest.mark.skipif(
        not joblib.parallel.mp, reason="joblib is in serial mode"
    )
    skip_if_array_api_compat_not_configured = pytest.mark.skipif(
        not ARRAY_API_COMPAT_FUNCTIONAL,
        reason="SCIPY_ARRAY_API not set, or versions of NumPy/SciPy too old.",
    )

    #  Decorator for tests involving both BLAS calls and multiprocessing.
    #
    #  Under POSIX (e.g. Linux or OSX), using multiprocessing in conjunction
    #  with some implementation of BLAS (or other libraries that manage an
    #  internal posix thread pool) can cause a crash or a freeze of the Python
    #  process.
    #
    #  In practice all known packaged distributions (from Linux distros or
    #  Anaconda) of BLAS under Linux seems to be safe. So we this problem seems
    #  to only impact OSX users.
    #
    #  This wrapper makes it possible to skip tests that can possibly cause
    #  this crash under OS X with.
    #
    #  Under Python 3.4+ it is possible to use the `forkserver` start method
    #  for multiprocessing to avoid this issue. However it can cause pickling
    #  errors on interactively defined functions. It therefore not enabled by
    #  default.

    if_safe_multiprocessing_with_blas = pytest.mark.skipif(
        sys.platform == "darwin", reason="Possible multi-process bug with some BLAS"
    )
    skip_if_no_numpydoc = pytest.mark.skipif(
        not _is_numpydoc(),
        reason="numpydoc is required to test the docstrings",
    )
except ImportError:
    pass


def check_skip_network():
    if int(os.environ.get("SKLEARN_SKIP_NETWORK_TESTS", 0)):
        raise SkipTest("Text tutorial requires large dataset download")


def _delete_folder(folder_path, warn=False):
    """Utility function to cleanup a temporary folder if still existing.

    Copy from joblib.pool (for independence).
    """
    try:
        if os.path.exists(folder_path):
            # This can fail under windows,
            #  but will succeed when called by atexit
            shutil.rmtree(folder_path)
    except OSError:
        if warn:
            warnings.warn("Could not delete temporary folder %s" % folder_path)


class TempMemmap:
    """
    Parameters
    ----------
    data
    mmap_mode : str, default='r'
    """

    def __init__(self, data, mmap_mode="r"):
        self.mmap_mode = mmap_mode
        self.data = data

    def __enter__(self):
        data_read_only, self.temp_folder = create_memmap_backed_data(
            self.data, mmap_mode=self.mmap_mode, return_folder=True
        )
        return data_read_only

    def __exit__(self, exc_type, exc_val, exc_tb):
        _delete_folder(self.temp_folder)


def create_memmap_backed_data(data, mmap_mode="r", return_folder=False):
    """
    Parameters
    ----------
    data
    mmap_mode : str, default='r'
    return_folder :  bool, default=False
    """
    temp_folder = tempfile.mkdtemp(prefix="sklearn_testing_")
    atexit.register(functools.partial(_delete_folder, temp_folder, warn=True))
    filename = op.join(temp_folder, "data.pkl")
    joblib.dump(data, filename)
    memmap_backed_data = joblib.load(filename, mmap_mode=mmap_mode)
    result = (
        memmap_backed_data if not return_folder else (memmap_backed_data, temp_folder)
    )
    return result


# Utils to test docstrings


def _get_args(function, varargs=False):
    """Helper to get function arguments."""

    try:
        params = signature(function).parameters
    except ValueError:
        # Error on builtin C function
        return []
    args = [
        key
        for key, param in params.items()
        if param.kind not in (param.VAR_POSITIONAL, param.VAR_KEYWORD)
    ]
    if varargs:
        varargs = [
            param.name
            for param in params.values()
            if param.kind == param.VAR_POSITIONAL
        ]
        if len(varargs) == 0:
            varargs = None
        return args, varargs
    else:
        return args


def _get_func_name(func):
    """Get function full name.

    Parameters
    ----------
    func : callable
        The function object.

    Returns
    -------
    name : str
        The function name.
    """
    parts = []
    module = inspect.getmodule(func)
    if module:
        parts.append(module.__name__)

    qualname = func.__qualname__
    if qualname != func.__name__:
        parts.append(qualname[: qualname.find(".")])

    parts.append(func.__name__)
    return ".".join(parts)


def check_docstring_parameters(func, doc=None, ignore=None):
    """Helper to check docstring.

    Parameters
    ----------
    func : callable
        The function object to test.
    doc : str, default=None
        Docstring if it is passed manually to the test.
    ignore : list, default=None
        Parameters to ignore.

    Returns
    -------
    incorrect : list
        A list of string describing the incorrect results.
    """
    from numpydoc import docscrape

    incorrect = []
    ignore = [] if ignore is None else ignore

    func_name = _get_func_name(func)
    if not func_name.startswith("sklearn.") or func_name.startswith(
        "sklearn.externals"
    ):
        return incorrect
    # Don't check docstring for property-functions
    if inspect.isdatadescriptor(func):
        return incorrect
    # Don't check docstring for setup / teardown pytest functions
    if func_name.split(".")[-1] in ("setup_module", "teardown_module"):
        return incorrect
    # Dont check estimator_checks module
    if func_name.split(".")[2] == "estimator_checks":
        return incorrect
    # Get the arguments from the function signature
    param_signature = list(filter(lambda x: x not in ignore, _get_args(func)))
    # drop self
    if len(param_signature) > 0 and param_signature[0] == "self":
        param_signature.remove("self")

    # Analyze function's docstring
    if doc is None:
        records = []
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error", UserWarning)
            try:
                doc = docscrape.FunctionDoc(func)
            except UserWarning as exp:
                if "potentially wrong underline length" in str(exp):
                    # Catch warning raised as of numpydoc 1.2 when
                    # the underline length for a section of a docstring
                    # is not consistent.
                    message = str(exp).split("\n")[:3]
                    incorrect += [f"In function: {func_name}"] + message
                    return incorrect
                records.append(str(exp))
            except Exception as exp:
                incorrect += [func_name + " parsing error: " + str(exp)]
                return incorrect
        if len(records):
            raise RuntimeError("Error for %s:\n%s" % (func_name, records[0]))

    param_docs = []
    for name, type_definition, param_doc in doc["Parameters"]:
        # Type hints are empty only if parameter name ended with :
        if not type_definition.strip():
            if ":" in name and name[: name.index(":")][-1:].strip():
                incorrect += [
                    func_name
                    + " There was no space between the param name and colon (%r)" % name
                ]
            elif name.rstrip().endswith(":"):
                incorrect += [
                    func_name
                    + " Parameter %r has an empty type spec. Remove the colon"
                    % (name.lstrip())
                ]

        # Create a list of parameters to compare with the parameters gotten
        # from the func signature
        if "*" not in name:
            param_docs.append(name.split(":")[0].strip("` "))

    # If one of the docstring's parameters had an error then return that
    # incorrect message
    if len(incorrect) > 0:
        return incorrect

    # Remove the parameters that should be ignored from list
    param_docs = list(filter(lambda x: x not in ignore, param_docs))

    # The following is derived from pytest, Copyright (c) 2004-2017 Holger
    # Krekel and others, Licensed under MIT License. See
    # https://github.com/pytest-dev/pytest

    message = []
    for i in range(min(len(param_docs), len(param_signature))):
        if param_signature[i] != param_docs[i]:
            message += [
                "There's a parameter name mismatch in function"
                " docstring w.r.t. function signature, at index %s"
                " diff: %r != %r" % (i, param_signature[i], param_docs[i])
            ]
            break
    if len(param_signature) > len(param_docs):
        message += [
            "Parameters in function docstring have less items w.r.t."
            " function signature, first missing item: %s"
            % param_signature[len(param_docs)]
        ]

    elif len(param_signature) < len(param_docs):
        message += [
            "Parameters in function docstring have more items w.r.t."
            " function signature, first extra item: %s"
            % param_docs[len(param_signature)]
        ]

    # If there wasn't any difference in the parameters themselves between
    # docstring and signature including having the same length then return
    # empty list
    if len(message) == 0:
        return []

    import difflib
    import pprint

    param_docs_formatted = pprint.pformat(param_docs).splitlines()
    param_signature_formatted = pprint.pformat(param_signature).splitlines()

    message += ["Full diff:"]

    message.extend(
        line.strip()
        for line in difflib.ndiff(param_signature_formatted, param_docs_formatted)
    )

    incorrect.extend(message)

    # Prepend function name
    incorrect = ["In function: " + func_name] + incorrect

    return incorrect


def _check_item_included(item_name, args):
    """Helper to check if item should be included in checking."""
    if args.include is not True and item_name not in args.include:
        return False
    if args.exclude is not None and item_name in args.exclude:
        return False
    return True


def _diff_key(line):
    """Key for grouping output from `context_diff`."""
    if line.startswith("  "):
        return "  "
    elif line.startswith("- "):
        return "- "
    elif line.startswith("+ "):
        return "+ "
    elif line.startswith("! "):
        return "! "
    return None


def _get_diff_msg(docstrings_grouped):
    """Get message showing the difference between type/desc docstrings of all objects.

    `docstrings_grouped` keys should be the type/desc docstrings and values are a list
    of objects with that docstring. Objects with the same type/desc docstring are
    thus grouped together.
    """
    msg_diff = ""
    ref_str = ""
    ref_group = []
    for docstring, group in docstrings_grouped.items():
        if not ref_str and not ref_group:
            ref_str += docstring
            ref_group.extend(group)
        diff = list(
            context_diff(
                ref_str.split(),
                docstring.split(),
                fromfile=str(ref_group),
                tofile=str(group),
                n=8,
            )
        )
        # Add header
        msg_diff += "".join((diff[:3]))
        # Group consecutive 'diff' words to shorten error message
        for start, group in groupby(diff[3:], key=_diff_key):
            if start is None:
                msg_diff += "\n" + "\n".join(group)
            else:
                msg_diff += "\n" + start + " ".join(word[2:] for word in group)
        # Add new lines at end of diff, to separate comparisons
        msg_diff += "\n\n"
    return msg_diff


def _check_consistency_items(
    items_docs,
    type_or_desc,
    section,
    n_objects,
    descr_regex_pattern="",
    ignore_types=tuple(),
):
    """Helper to check docstring consistency of all `items_docs`.

    If item is not present in all objects, checking is skipped and warning raised.
    If `regex` provided, match descriptions to all descriptions.

    Parameters
    ----------
    items_doc : dict of dict of str
        Dictionary where the key is the string type or description, value is
        a dictionary where the key is "type description" or "description"
        and the value is a list of object names with the same string type or
        description.

    type_or_desc : {"type description", "description"}
        Whether to check type description or description between objects.

    section : {"Parameters", "Attributes", "Returns"}
        Name of the section type.

    n_objects : int
        Total number of objects.

    descr_regex_pattern : str, default=""
        Regex pattern to match for description of all objects.
        Ignored when `type_or_desc="type description".

    ignore_types : tuple of str, default=()
        Tuple of parameter/attribute/return names for which type description
        matching is ignored. Ignored when `type_or_desc="description".
    """
    skipped = []
    for item_name, docstrings_grouped in items_docs.items():
        # If item not found in all objects, skip
        if sum([len(objs) for objs in docstrings_grouped.values()]) < n_objects:
            skipped.append(item_name)
        # If regex provided, match to all descriptions
        elif type_or_desc == "description" and descr_regex_pattern:
            not_matched = []
            for docstring, group in docstrings_grouped.items():
                if not re.search(descr_regex_pattern, docstring):
                    not_matched.extend(group)
            if not_matched:
                msg = textwrap.fill(
                    f"The description of {section[:-1]} '{item_name}' in {not_matched}"
                    f" does not match 'descr_regex_pattern': {descr_regex_pattern} "
                )
                raise AssertionError(msg)
        # Skip type checking for items in `ignore_types`
        elif type_or_desc == "type specification" and item_name in ignore_types:
            continue
        # Otherwise, if more than one key, docstrings not consistent between objects
        elif len(docstrings_grouped.keys()) > 1:
            msg_diff = _get_diff_msg(docstrings_grouped)
            obj_groups = " and ".join(
                str(group) for group in docstrings_grouped.values()
            )
            msg = textwrap.fill(
                f"The {type_or_desc} of {section[:-1]} '{item_name}' is inconsistent "
                f"between {obj_groups}:"
            )
            msg += msg_diff
            raise AssertionError(msg)
    if skipped:
        warnings.warn(
            f"Checking was skipped for {section}: {skipped} as they were "
            "not found in all objects."
        )


def assert_docstring_consistency(
    objects,
    include_params=False,
    exclude_params=None,
    include_attrs=False,
    exclude_attrs=None,
    include_returns=False,
    exclude_returns=None,
    descr_regex_pattern=None,
    ignore_types=tuple(),
):
    r"""Check consistency between docstring parameters/attributes/returns of objects.

    Checks if parameters/attributes/returns have the same type specification and
    description (ignoring whitespace) across `objects`. Intended to be used for
    related classes/functions/data descriptors.

    Entries that do not appear across all `objects` are ignored.

    Parameters
    ----------
    objects : list of {classes, functions, data descriptors}
        Objects to check.
        Objects may be classes, functions or data descriptors with docstrings that
        can be parsed by numpydoc.

    include_params : list of str or bool, default=False
        List of parameters to be included. If True, all parameters are included,
        if False, checking is skipped for parameters.
        Can only be set if `exclude_params` is None.

    exclude_params : list of str or None, default=None
        List of parameters to be excluded. If None, no parameters are excluded.
        Can only be set if `include_params` is True.

    include_attrs : list of str or bool, default=False
        List of attributes to be included. If True, all attributes are included,
        if False, checking is skipped for attributes.
        Can only be set if `exclude_attrs` is None.

    exclude_attrs : list of str or None, default=None
        List of attributes to be excluded. If None, no attributes are excluded.
        Can only be set if `include_attrs` is True.

    include_returns : list of str or bool, default=False
        List of returns to be included. If True, all returns are included,
        if False, checking is skipped for returns.
        Can only be set if `exclude_returns` is None.

    exclude_returns : list of str or None, default=None
        List of returns to be excluded. If None, no returns are excluded.
        Can only be set if `include_returns` is True.

    descr_regex_pattern : str, default=None
        Regular expression to match to all descriptions of included
        parameters/attributes/returns. If None, will revert to default behavior
        of comparing descriptions between objects.

    ignore_types : tuple of str, default=tuple()
        Tuple of parameter/attribute/return names to exclude from type description
        matching between objects.

    Examples
    --------
    >>> from sklearn.metrics import (accuracy_score, classification_report,
    ... mean_absolute_error, mean_squared_error, median_absolute_error)
    >>> from sklearn.utils._testing import assert_docstring_consistency
    ... # doctest: +SKIP
    >>> assert_docstring_consistency([mean_absolute_error, mean_squared_error],
    ... include_params=['y_true', 'y_pred', 'sample_weight'])  # doctest: +SKIP
    >>> assert_docstring_consistency([median_absolute_error, mean_squared_error],
    ... include_params=True)  # doctest: +SKIP
    >>> assert_docstring_consistency([accuracy_score, classification_report],
    ... include_params=["y_true"],
    ... descr_regex_pattern=r"Ground truth \(correct\) (labels|target values)")
    ... # doctest: +SKIP
    """
    from numpydoc.docscrape import NumpyDocString

    Args = namedtuple("args", ["include", "exclude", "arg_name"])

    def _create_args(include, exclude, arg_name, section_name):
        if exclude and include is not True:
            raise TypeError(
                f"The 'exclude_{arg_name}' argument can be set only when the "
                f"'include_{arg_name}' argument is True."
            )
        if include is False:
            return {}
        return {section_name: Args(include, exclude, arg_name)}

    section_args = {
        **_create_args(include_params, exclude_params, "params", "Parameters"),
        **_create_args(include_attrs, exclude_attrs, "attrs", "Attributes"),
        **_create_args(include_returns, exclude_returns, "returns", "Returns"),
    }

    objects_doc = dict()
    for obj in objects:
        if (
            inspect.isdatadescriptor(obj)
            or inspect.isfunction(obj)
            or inspect.isclass(obj)
        ):
            objects_doc[obj.__name__] = NumpyDocString(inspect.getdoc(obj))
        else:
            raise TypeError(
                "All 'objects' must be one of: function, class or descriptor, "
                f"got a: {type(obj)}."
            )

    n_objects = len(objects)
    for section, args in section_args.items():
        type_items = defaultdict(lambda: defaultdict(list))
        desc_items = defaultdict(lambda: defaultdict(list))
        for obj_name, obj_doc in objects_doc.items():
            for item_name, type_def, desc in obj_doc[section]:
                if _check_item_included(item_name, args):
                    # Normalize white space
                    type_def = " ".join(type_def.strip().split())
                    desc = " ".join(chain.from_iterable(line.split() for line in desc))
                    # Use string type/desc as key, to group consistent objs together
                    type_items[item_name][type_def].append(obj_name)
                    desc_items[item_name][desc].append(obj_name)

        _check_consistency_items(
            type_items,
            "type specification",
            section,
            n_objects,
            ignore_types=ignore_types,
        )
        _check_consistency_items(
            desc_items,
            "description",
            section,
            n_objects,
            descr_regex_pattern=descr_regex_pattern,
        )


def assert_run_python_script_without_output(source_code, pattern=".+", timeout=60):
    """Utility to check assertions in an independent Python subprocess.

    The script provided in the source code should return 0 and the stdtout +
    stderr should not match the pattern `pattern`.

    This is a port from cloudpickle https://github.com/cloudpipe/cloudpickle

    Parameters
    ----------
    source_code : str
        The Python source code to execute.
    pattern : str
        Pattern that the stdout + stderr should not match. By default, unless
        stdout + stderr are both empty, an error will be raised.
    timeout : int, default=60
        Time in seconds before timeout.
    """
    fd, source_file = tempfile.mkstemp(suffix="_src_test_sklearn.py")
    os.close(fd)
    try:
        with open(source_file, "wb") as f:
            f.write(source_code.encode("utf-8"))
        cmd = [sys.executable, source_file]
        cwd = op.normpath(op.join(op.dirname(sklearn_path), ".."))
        env = os.environ.copy()
        try:
            env["PYTHONPATH"] = os.pathsep.join([cwd, env["PYTHONPATH"]])
        except KeyError:
            env["PYTHONPATH"] = cwd
        kwargs = {"cwd": cwd, "stderr": STDOUT, "env": env}
        # If coverage is running, pass the config file to the subprocess
        coverage_rc = os.environ.get("COVERAGE_PROCESS_START")
        if coverage_rc:
            kwargs["env"]["COVERAGE_PROCESS_START"] = coverage_rc

        kwargs["timeout"] = timeout
        try:
            try:
                out = check_output(cmd, **kwargs)
            except CalledProcessError as e:
                raise RuntimeError(
                    "script errored with output:\n%s" % e.output.decode("utf-8")
                )

            out = out.decode("utf-8")
            if re.search(pattern, out):
                if pattern == ".+":
                    expectation = "Expected no output"
                else:
                    expectation = f"The output was not supposed to match {pattern!r}"

                message = f"{expectation}, got the following output instead: {out!r}"
                raise AssertionError(message)
        except TimeoutExpired as e:
            raise RuntimeError(
                "script timeout, output so far:\n%s" % e.output.decode("utf-8")
            )
    finally:
        os.unlink(source_file)


def _convert_container(
    container,
    constructor_name,
    columns_name=None,
    dtype=None,
    minversion=None,
    categorical_feature_names=None,
):
    """Convert a given container to a specific array-like with a dtype.

    Parameters
    ----------
    container : array-like
        The container to convert.
    constructor_name : {"list", "tuple", "array", "sparse", "dataframe", \
            "pandas", "series", "index", "slice", "sparse_csr", "sparse_csc", \
            "sparse_csr_array", "sparse_csc_array", "pyarrow", "polars", \
            "polars_series"}
        The type of the returned container.
    columns_name : index or array-like, default=None
        For pandas/polars container supporting `columns_names`, it will affect
        specific names.
    dtype : dtype, default=None
        Force the dtype of the container. Does not apply to `"slice"`
        container.
    minversion : str, default=None
        Minimum version for package to install.
    categorical_feature_names : list of str, default=None
        List of column names to cast to categorical dtype.

    Returns
    -------
    converted_container
    """
    if constructor_name == "list":
        if dtype is None:
            return list(container)
        else:
            return np.asarray(container, dtype=dtype).tolist()
    elif constructor_name == "tuple":
        if dtype is None:
            return tuple(container)
        else:
            return tuple(np.asarray(container, dtype=dtype).tolist())
    elif constructor_name == "array":
        return np.asarray(container, dtype=dtype)
    elif constructor_name in ("pandas", "dataframe"):
        pd = pytest.importorskip("pandas", minversion=minversion)
        result = pd.DataFrame(container, columns=columns_name, dtype=dtype, copy=False)
        if categorical_feature_names is not None:
            for col_name in categorical_feature_names:
                result[col_name] = result[col_name].astype("category")
        return result
    elif constructor_name == "pyarrow":
        pa = pytest.importorskip("pyarrow", minversion=minversion)
        array = np.asarray(container)
        array = array[:, None] if array.ndim == 1 else array
        if columns_name is None:
            columns_name = [f"col{i}" for i in range(array.shape[1])]
        data = {name: array[:, i] for i, name in enumerate(columns_name)}
        result = pa.Table.from_pydict(data)
        if categorical_feature_names is not None:
            for col_idx, col_name in enumerate(result.column_names):
                if col_name in categorical_feature_names:
                    result = result.set_column(
                        col_idx, col_name, result.column(col_name).dictionary_encode()
                    )
        return result
    elif constructor_name == "polars":
        pl = pytest.importorskip("polars", minversion=minversion)
        result = pl.DataFrame(container, schema=columns_name, orient="row")
        if categorical_feature_names is not None:
            for col_name in categorical_feature_names:
                result = result.with_columns(pl.col(col_name).cast(pl.Categorical))
        return result
    elif constructor_name == "series":
        pd = pytest.importorskip("pandas", minversion=minversion)
        return pd.Series(container, dtype=dtype)
    elif constructor_name == "pyarrow_array":
        pa = pytest.importorskip("pyarrow", minversion=minversion)
        return pa.array(container)
    elif constructor_name == "polars_series":
        pl = pytest.importorskip("polars", minversion=minversion)
        return pl.Series(values=container)
    elif constructor_name == "index":
        pd = pytest.importorskip("pandas", minversion=minversion)
        return pd.Index(container, dtype=dtype)
    elif constructor_name == "slice":
        return slice(container[0], container[1])
    elif "sparse" in constructor_name:
        if not sp.sparse.issparse(container):
            # For scipy >= 1.13, sparse array constructed from 1d array may be
            # 1d or raise an exception. To avoid this, we make sure that the
            # input container is 2d. For more details, see
            # https://github.com/scipy/scipy/pull/18530#issuecomment-1878005149
            container = np.atleast_2d(container)

        if constructor_name in ("sparse", "sparse_csr"):
            # sparse and sparse_csr are equivalent for legacy reasons
            return sp.sparse.csr_matrix(container, dtype=dtype)
        elif constructor_name == "sparse_csr_array":
            return sp.sparse.csr_array(container, dtype=dtype)
        elif constructor_name == "sparse_csc":
            return sp.sparse.csc_matrix(container, dtype=dtype)
        elif constructor_name == "sparse_csc_array":
            return sp.sparse.csc_array(container, dtype=dtype)


def raises(expected_exc_type, match=None, may_pass=False, err_msg=None):
    """Context manager to ensure exceptions are raised within a code block.

    This is similar to and inspired from pytest.raises, but supports a few
    other cases.

    This is only intended to be used in estimator_checks.py where we don't
    want to use pytest. In the rest of the code base, just use pytest.raises
    instead.

    Parameters
    ----------
    excepted_exc_type : Exception or list of Exception
        The exception that should be raised by the block. If a list, the block
        should raise one of the exceptions.
    match : str or list of str, default=None
        A regex that the exception message should match. If a list, one of
        the entries must match. If None, match isn't enforced.
    may_pass : bool, default=False
        If True, the block is allowed to not raise an exception. Useful in
        cases where some estimators may support a feature but others must
        fail with an appropriate error message. By default, the context
        manager will raise an exception if the block does not raise an
        exception.
    err_msg : str, default=None
        If the context manager fails (e.g. the block fails to raise the
        proper exception, or fails to match), then an AssertionError is
        raised with this message. By default, an AssertionError is raised
        with a default error message (depends on the kind of failure). Use
        this to indicate how users should fix their estimators to pass the
        checks.

    Attributes
    ----------
    raised_and_matched : bool
        True if an exception was raised and a match was found, False otherwise.
    """
    return _Raises(expected_exc_type, match, may_pass, err_msg)


class _Raises(contextlib.AbstractContextManager):
    # see raises() for parameters
    def __init__(self, expected_exc_type, match, may_pass, err_msg):
        self.expected_exc_types = (
            expected_exc_type
            if isinstance(expected_exc_type, Iterable)
            else [expected_exc_type]
        )
        self.matches = [match] if isinstance(match, str) else match
        self.may_pass = may_pass
        self.err_msg = err_msg
        self.raised_and_matched = False

    def __exit__(self, exc_type, exc_value, _):
        # see
        # https://docs.python.org/2.5/whatsnew/pep-343.html#SECTION000910000000000000000

        if exc_type is None:  # No exception was raised in the block
            if self.may_pass:
                return True  # CM is happy
            else:
                err_msg = self.err_msg or f"Did not raise: {self.expected_exc_types}"
                raise AssertionError(err_msg)

        if not any(
            issubclass(exc_type, expected_type)
            for expected_type in self.expected_exc_types
        ):
            if self.err_msg is not None:
                raise AssertionError(self.err_msg) from exc_value
            else:
                return False  # will re-raise the original exception

        if self.matches is not None:
            err_msg = self.err_msg or (
                "The error message should contain one of the following "
                "patterns:\n{}\nGot {}".format("\n".join(self.matches), str(exc_value))
            )
            if not any(re.search(match, str(exc_value)) for match in self.matches):
                raise AssertionError(err_msg) from exc_value
            self.raised_and_matched = True

        return True


class MinimalClassifier:
    """Minimal classifier implementation without inheriting from BaseEstimator.

    This estimator should be tested with:

    * `check_estimator` in `test_estimator_checks.py`;
    * within a `Pipeline` in `test_pipeline.py`;
    * within a `SearchCV` in `test_search.py`.
    """

    def __init__(self, param=None):
        self.param = param

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        check_classification_targets(y)
        self.classes_, counts = np.unique(y, return_counts=True)
        self._most_frequent_class_idx = counts.argmax()
        return self

    def predict_proba(self, X):
        check_is_fitted(self)
        X = check_array(X)
        proba_shape = (X.shape[0], self.classes_.size)
        y_proba = np.zeros(shape=proba_shape, dtype=np.float64)
        y_proba[:, self._most_frequent_class_idx] = 1.0
        return y_proba

    def predict(self, X):
        y_proba = self.predict_proba(X)
        y_pred = y_proba.argmax(axis=1)
        return self.classes_[y_pred]

    def score(self, X, y):
        from sklearn.metrics import accuracy_score

        return accuracy_score(y, self.predict(X))

    def __sklearn_tags__(self):
        return Tags(
            estimator_type="classifier",
            classifier_tags=ClassifierTags(),
            regressor_tags=None,
            transformer_tags=None,
            target_tags=TargetTags(required=True),
        )


class MinimalRegressor:
    """Minimal regressor implementation without inheriting from BaseEstimator.

    This estimator should be tested with:

    * `check_estimator` in `test_estimator_checks.py`;
    * within a `Pipeline` in `test_pipeline.py`;
    * within a `SearchCV` in `test_search.py`.
    """

    def __init__(self, param=None):
        self.param = param

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        self.is_fitted_ = True
        self._mean = np.mean(y)
        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)
        return np.ones(shape=(X.shape[0],)) * self._mean

    def score(self, X, y):
        from sklearn.metrics import r2_score

        return r2_score(y, self.predict(X))

    def __sklearn_tags__(self):
        return Tags(
            estimator_type="regressor",
            classifier_tags=None,
            regressor_tags=RegressorTags(),
            transformer_tags=None,
            target_tags=TargetTags(required=True),
        )


class MinimalTransformer:
    """Minimal transformer implementation without inheriting from
    BaseEstimator.

    This estimator should be tested with:

    * `check_estimator` in `test_estimator_checks.py`;
    * within a `Pipeline` in `test_pipeline.py`;
    * within a `SearchCV` in `test_search.py`.
    """

    def __init__(self, param=None):
        self.param = param

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def fit(self, X, y=None):
        check_array(X)
        self.is_fitted_ = True
        return self

    def transform(self, X, y=None):
        check_is_fitted(self)
        X = check_array(X)
        return X

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X, y)

    def __sklearn_tags__(self):
        return Tags(
            estimator_type="transformer",
            classifier_tags=None,
            regressor_tags=None,
            transformer_tags=TransformerTags(),
            target_tags=TargetTags(required=False),
        )


def _array_api_for_tests(array_namespace, device):
    try:
        array_mod = importlib.import_module(array_namespace)
    except (ModuleNotFoundError, ImportError):
        raise SkipTest(
            f"{array_namespace} is not installed: not checking array_api input"
        )

    if os.environ.get("SCIPY_ARRAY_API") is None:
        raise SkipTest("SCIPY_ARRAY_API is not set: not checking array_api input")

    from sklearn.externals.array_api_compat import get_namespace

    # First create an array using the chosen array module and then get the
    # corresponding (compatibility wrapped) array namespace based on it.
    # This is because `cupy` is not the same as the compatibility wrapped
    # namespace of a CuPy array.
    xp = get_namespace(array_mod.asarray(1))
    if (
        array_namespace == "torch"
        and device == "cuda"
        and not xp.backends.cuda.is_built()
    ):
        raise SkipTest("PyTorch test requires cuda, which is not available")
    elif array_namespace == "torch" and device == "mps":
        if os.getenv("PYTORCH_ENABLE_MPS_FALLBACK") != "1":
            # For now we need PYTORCH_ENABLE_MPS_FALLBACK=1 for all estimators to work
            # when using the MPS device.
            raise SkipTest(
                "Skipping MPS device test because PYTORCH_ENABLE_MPS_FALLBACK is not "
                "set."
            )
        if not xp.backends.mps.is_built():
            raise SkipTest(
                "MPS is not available because the current PyTorch install was not "
                "built with MPS enabled."
            )
    elif array_namespace == "torch" and device == "xpu":  # pragma: nocover
        if not hasattr(xp, "xpu"):
            # skip xpu testing for PyTorch <2.4
            raise SkipTest(
                "XPU is not available because the current PyTorch install was not "
                "built with XPU support."
            )
        if not xp.xpu.is_available():
            raise SkipTest(
                "Skipping XPU device test because no XPU device is available"
            )
    elif array_namespace == "cupy":  # pragma: nocover
        import cupy

        if cupy.cuda.runtime.getDeviceCount() == 0:
            raise SkipTest("CuPy test requires cuda, which is not available")
    return xp


def _get_warnings_filters_info_list():
    @dataclass
    class WarningInfo:
        action: "warnings._ActionKind"  # type: ignore[annotation-unchecked]
        message: str = ""  # type: ignore[annotation-unchecked]
        category: type[Warning] = Warning  # type: ignore[annotation-unchecked]

        def to_filterwarning_str(self):
            if self.category.__module__ == "builtins":
                category = self.category.__name__
            else:
                category = f"{self.category.__module__}.{self.category.__name__}"

            return f"{self.action}:{self.message}:{category}"

    return [
        WarningInfo("error", category=DeprecationWarning),
        WarningInfo("error", category=FutureWarning),
        WarningInfo("error", category=VisibleDeprecationWarning),
        # TODO: remove when pyamg > 5.0.1
        # Avoid a deprecation warning due pkg_resources usage in pyamg.
        WarningInfo(
            "ignore",
            message="pkg_resources is deprecated as an API",
            category=DeprecationWarning,
        ),
        WarningInfo(
            "ignore",
            message="Deprecated call to `pkg_resources",
            category=DeprecationWarning,
        ),
        # pytest-cov issue https://github.com/pytest-dev/pytest-cov/issues/557 not
        # fixed although it has been closed. https://github.com/pytest-dev/pytest-cov/pull/623
        # would probably fix it.
        WarningInfo(
            "ignore",
            message=(
                "The --rsyncdir command line argument and rsyncdirs config variable are"
                " deprecated"
            ),
            category=DeprecationWarning,
        ),
        # XXX: Easiest way to ignore pandas Pyarrow DeprecationWarning in the
        # short-term. See https://github.com/pandas-dev/pandas/issues/54466 for
        # more details.
        WarningInfo(
            "ignore",
            message=r"\s*Pyarrow will become a required dependency",
            category=DeprecationWarning,
        ),
        # warnings has been fixed from dateutil main but not released yet, see
        # https://github.com/dateutil/dateutil/issues/1314
        WarningInfo(
            "ignore",
            message="datetime.datetime.utcfromtimestamp",
            category=DeprecationWarning,
        ),
        # Python 3.12 warnings from joblib fixed in master but not released yet,
        # see https://github.com/joblib/joblib/pull/1518
        WarningInfo(
            "ignore", message="ast.Num is deprecated", category=DeprecationWarning
        ),
        WarningInfo(
            "ignore", message="Attribute n is deprecated", category=DeprecationWarning
        ),
        # Python 3.12 warnings from sphinx-gallery fixed in master but not
        # released yet, see
        # https://github.com/sphinx-gallery/sphinx-gallery/pull/1242
        WarningInfo(
            "ignore", message="ast.Str is deprecated", category=DeprecationWarning
        ),
        WarningInfo(
            "ignore", message="Attribute s is deprecated", category=DeprecationWarning
        ),
        # Plotly deprecated something which we're not using, but internally it's used
        # and needs to be fixed on their side.
        # https://github.com/plotly/plotly.py/issues/4997
        WarningInfo(
            "ignore",
            message=".+scattermapbox.+deprecated.+scattermap.+instead",
            category=DeprecationWarning,
        ),
        # TODO(1.10): remove PassiveAgressive
        WarningInfo(
            "ignore",
            message="Class PassiveAggressive.+is deprecated",
            category=FutureWarning,
        ),
    ]


def get_pytest_filterwarning_lines():
    warning_filters_info_list = _get_warnings_filters_info_list()
    return [
        warning_info.to_filterwarning_str()
        for warning_info in warning_filters_info_list
    ]


def turn_warnings_into_errors():
    warnings_filters_info_list = _get_warnings_filters_info_list()
    for warning_info in warnings_filters_info_list:
        warnings.filterwarnings(
            warning_info.action,
            message=warning_info.message,
            category=warning_info.category,
        )
