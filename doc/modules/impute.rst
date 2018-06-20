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


Univariate feature imputation
=============================

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

.. _mice:

Multivariate feature imputation
===============================

A more sophisticated approach is to use the :class:`ChainedImputer` class, which
implements MICE: Multivariate Imputation by Chained Equations. MICE is usually used
to generate multiple imputations, but :class:`ChainedImputer` generates a single
(averaged) imputation, which is why it is not named `MICEImputer`. MICE models each
feature with missing values as a function of other features, and uses that estimate
for imputation. It does so in a round-robin fashion: at each step, a feature column
is designated as output `y` and the other feature columns are treated as inputs `X`.
A regressor is fit on `(X, y)` for known `y`. Then, the regressor is used to predict
the unknown values of `y`. This is repeated for each feature, and then is done for a
number of imputation rounds, in a chained fashion. Here is an example snippet::

    >>> import numpy as np
    >>> from sklearn.impute import ChainedImputer
    >>> imp = ChainedImputer(n_imputations=10, random_state=0)
    >>> imp.fit([[1, 2], [np.nan, 3], [7, np.nan]])
    ChainedImputer(imputation_order='ascending', initial_strategy='mean',
            max_value=None, min_value=None, missing_values='NaN', n_burn_in=10,
            n_imputations=10, n_nearest_features=None, predictor=None,
            random_state=0, verbose=False)
    >>> X_test = [[np.nan, 2], [6, np.nan], [np.nan, 6]]
    >>> print(np.round(imp.transform(X_test)))
    [[ 1.  2.]
     [ 6.  4.]
     [13.  6.]]

Note that MICE has traditionally been used in the statistics community to obtain
multiple imputations, to which the entire subsequent analysis pipeline is applied.
This allows one to estimate the variance that is inherent in having missing values
in the first place. In the above example, however, we are using :class:`ChainedImputer`
in way that is common for machine learning applications: to fill-in missing values
as accurately as possible. We thus we average the last 10 imputations to obtain a
single one.

Both :class:`SimpleImputer` and :class:`ChainedImputer` can be used in a Pipeline
as a way to build a composite estimator that supports imputation.
See :ref:`sphx_glr_auto_examples_plot_missing_values.py`.
