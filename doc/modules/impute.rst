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
i.e., to infer them from the known part of the data.

The :class:`SimpleImputer` class provides basic strategies for imputing missing
values, either using the mean, the median or the most frequent value of
the row or column in which the missing values are located. This class
also allows for different missing values encodings.

The following snippet demonstrates how to replace missing values,
encoded as ``np.nan``, using the mean value of the columns (axis 0)
that contain the missing values::

    >>> import numpy as np
    >>> from sklearn.impute import SimpleImputer
    >>> imp = SimpleImputer(missing_values='NaN', strategy='mean')
    >>> imp.fit([[1, 2], [np.nan, 3], [7, 6]])       # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, missing_values='NaN', strategy='mean', verbose=0)
    >>> X = [[np.nan, 2], [6, np.nan], [7, 6]]
    >>> print(imp.transform(X))           # doctest: +NORMALIZE_WHITESPACE  +ELLIPSIS
    [[4.          2.        ]
     [6.          3.666...]
     [7.          6.        ]]

The :class:`SimpleImputer` class also supports sparse matrices::

    >>> import scipy.sparse as sp
    >>> X = sp.csc_matrix([[1, 2], [0, 3], [7, 6]])
    >>> imp = SimpleImputer(missing_values=0, strategy='mean')
    >>> imp.fit(X)                  # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, missing_values=0, strategy='mean', verbose=0)
    >>> X_test = sp.csc_matrix([[0, 2], [6, 0], [7, 6]])
    >>> print(imp.transform(X_test))      # doctest: +NORMALIZE_WHITESPACE  +ELLIPSIS
    [[4.          2.        ]
     [6.          3.666...]
     [7.          6.        ]]

Note that, here, missing values are encoded by 0 and are thus implicitly stored
in the matrix. This format is thus suitable when there are many more missing
values than observed values.

A more sophisticated approach is to use the :class:`MICEImputer` class, which
implements the multivariate imputation by chained equations technique. While
the :class:`SimpleImputer` class fills in each missing value in a feature
column by a constant that is a function of that feature only, MICE fills in
each feature's missing values based on the values of the other features as
well. It does so in a round-robin fashion: at each step, a feature column is
designated as output and the other feature columns are treated as inputs. A
regressor is fit on the non-missing values of the output feature column and the
corresponding values in the input columns. Then, the regressor is applied so as
to predict the unknown values of the output feature. This is repeated for each
feature, and then is done for a number of imputation rounds. Here is an example
snippet::

    >>> import numpy as np
    >>> from sklearn.impute import MICEImputer
    >>> imp = MICEImputer(n_imputations=10, random_state=0)
    >>> imp.fit([[1, 2], [np.nan, 3], [7, np.nan]])
    MICEImputer(imputation_order='ascending', initial_strategy='mean',
          max_value=None, min_value=None, missing_values='NaN', n_burn_in=10,
          n_imputations=10, n_nearest_features=None, predictor=None,
          random_state=0, verbose=False)
    >>> X_test = [[np.nan, 2], [6, np.nan], [np.nan, 6]]
    >>> print(np.round(imp.transform(X_test)))
    [[ 1.  2.]
     [ 6.  4.]
     [13.  6.]]

Note that MICE uses sampling as part of its imputation process, so the exact
values of the imputed data will be slightly different each time you run it.

Both :class:`SimpleImputer` and :class:`MICEImputer` can be used in a Pipeline
as a way to build a composite estimator that supports imputation.
See :ref:`sphx_glr_auto_examples_plot_missing_values.py`.
