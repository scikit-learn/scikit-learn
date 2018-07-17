.. _impute:

============================
Imputation of missing values
============================

.. currentmodule:: sklearn.impute

For various reasons, many real world datasets contain missing values, often
encoded as blanks, NaNs or other placeholders. Such datasets however are
incompatible with scikit-learn estimators which assume that all values in an
array are numerical, and that all have and hold meaning. A basic strategy to use
incomplete datasets is to discard entire rows and/or columns containing missing
values. However, this comes at the price of losing data which may be valuable
(even though incomplete). A better strategy is to impute the missing values,
i.e., to infer them from the known part of the data. See the :ref:`glossary`
entry on imputation.

The :class:`SimpleImputer` class provides basic strategies for imputing missing
values. Missing values can be imputed with a provided constant value, or using
the statistics (mean, median or most frequent) of each column in which the
missing values are located. This class also allows for different missing values
encodings.

The following snippet demonstrates how to replace missing values,
encoded as ``np.nan``, using the mean value of the columns (axis 0)
that contain the missing values::

    >>> import numpy as np
    >>> from sklearn.impute import SimpleImputer
    >>> imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    >>> imp.fit([[1, 2], [np.nan, 3], [7, 6]])       # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, fill_value=None, missing_values=nan, strategy='mean', verbose=0)
    >>> X = [[np.nan, 2], [6, np.nan], [7, 6]]
    >>> print(imp.transform(X))           # doctest: +NORMALIZE_WHITESPACE  +ELLIPSIS
    [[4.          2.        ]
     [6.          3.666...]
     [7.          6.        ]]

The :class:`SimpleImputer` class also supports sparse matrices::

    >>> import scipy.sparse as sp
    >>> X = sp.csc_matrix([[1, 2], [0, -1], [8, 4]])
    >>> imp = SimpleImputer(missing_values=-1, strategy='mean')
    >>> imp.fit(X)                  # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, fill_value=None, missing_values=-1, strategy='mean', verbose=0)
    >>> X_test = sp.csc_matrix([[-1, 2], [6, -1], [7, 6]])
    >>> print(imp.transform(X_test).toarray())      # doctest: +NORMALIZE_WHITESPACE
    [[3. 2.]
     [6. 3.]
     [7. 6.]]

Note that this format is not meant to be used to implicitly store missing values
in the matrix because it would densify it at transform time. Missing values encoded
by 0 must be used with dense input.

The :class:`SimpleImputer` class also supports categorical data represented as
string values or pandas categoricals when using the ``'most_frequent'`` or
``'constant'`` strategy::

    >>> import pandas as pd
    >>> df = pd.DataFrame([["a", "x"],
    ...                    [np.nan, "y"],
    ...                    ["a", np.nan],
    ...                    ["b", "y"]], dtype="category")
    ...
    >>> imp = SimpleImputer(strategy="most_frequent")
    >>> print(imp.fit_transform(df))      # doctest: +NORMALIZE_WHITESPACE
    [['a' 'x']
     ['a' 'y']
     ['a' 'y']
     ['b' 'y']]


:class:`SimpleImputer` can be used in a Pipeline as a way to build a composite
estimator that supports imputation. See :ref:`sphx_glr_auto_examples_plot_missing_values.py`.

.. _missing_indicator:

Marking imputed values
======================

The :class:`MissingIndicator` transformer is useful to transform a dataset into
corresponding binary matrix indicating the presence of missing values in the
dataset. This transformation is useful in conjunction with imputation. When
using imputation, preserving the information about which values had been
missing can be informative.

``NaN`` is usually used as the placeholder for missing values. However, it
enforces the data type to be float. The parameter ``missing_values`` allows to
specify other placeholder such as integer. In the following example, we will
use ``-1`` as missing values::

  >>> from sklearn.impute import MissingIndicator
  >>> X = np.array([[-1, -1, 1, 3],
  ...               [4, -1, 0, -1],
  ...               [8, -1, 1, 0]])
  >>> indicator = MissingIndicator(missing_values=-1)
  >>> mask_missing_values_only = indicator.fit_transform(X)
  >>> mask_missing_values_only
  array([[ True,  True, False],
         [False,  True,  True],
         [False,  True, False]])

The ``features`` parameter is used to choose the features for which the mask is
constructed. By default, it is ``'missing-only'`` which returns the imputer
mask of the features containing missing values at ``fit`` time::

  >>> indicator.features_
  array([0, 1, 3])

The ``features`` parameter can be set to ``'all'`` to returned all features
whether or not they contain missing values::
    
  >>> indicator = MissingIndicator(missing_values=-1, features="all")
  >>> mask_all = indicator.fit_transform(X)
  >>> mask_all
  array([[ True,  True, False, False],
         [False,  True, False,  True],
         [False,  True, False, False]])
  >>> indicator.features_
  array([0, 1, 2, 3])
