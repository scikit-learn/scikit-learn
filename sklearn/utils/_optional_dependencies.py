# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


def check_matplotlib_support(caller_name):
    """Raise ImportError with detailed error message if mpl is not installed.

    Plot utilities like any of the Display's plotting functions should lazily import
    matplotlib and call this helper before any computation.

    Parameters
    ----------
    caller_name : str
        The name of the caller that requires matplotlib.
    """
    try:
        import matplotlib  # noqa
    except ImportError as e:
        raise ImportError(
            "{} requires matplotlib. You can install matplotlib with "
            "`pip install matplotlib`".format(caller_name)
        ) from e


def check_pandas_support(caller_name):
    """Raise ImportError with detailed error message if pandas is not installed.

    Plot utilities like :func:`fetch_openml` should lazily import
    pandas and call this helper before any computation.

    Parameters
    ----------
    caller_name : str
        The name of the caller that requires pandas.

    Returns
    -------
    pandas
        The pandas package.
    """
    try:
        import pandas  # noqa

        return pandas
    except ImportError as e:
        raise ImportError("{} requires pandas.".format(caller_name)) from e
