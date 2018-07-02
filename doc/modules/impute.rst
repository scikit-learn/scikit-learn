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


Univariate vs. Multivariate Imputation
======================================

One type of imputation algorithm is univariate, which imputes values in the i-th
feature dimension using only non-missing values in that feature dimension
(e.g. :class:`impute.SimpleImputer`). By contrast, multivariate imputation
algorithms use the entire set of available feature dimensions to estimate the
missing values (e.g. :class:`impute.IterativeImputer`).


.. _single_imputer:

Univariate feature imputation
=============================

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
    >>> X = sp.csc_matrix([[1, 2], [0, 3], [7, 6]])
    >>> imp = SimpleImputer(missing_values=0, strategy='mean')
    >>> imp.fit(X)                  # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, fill_value=None, missing_values=0, strategy='mean', verbose=0)
    >>> X_test = sp.csc_matrix([[0, 2], [6, 0], [7, 6]])
    >>> print(imp.transform(X_test))      # doctest: +NORMALIZE_WHITESPACE  +ELLIPSIS
    [[4.          2.        ]
     [6.          3.666...]
     [7.          6.        ]]

Note that, here, missing values are encoded by 0 and are thus implicitly stored
in the matrix. This format is thus suitable when there are many more missing
values than observed values.

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

.. _iterative_imputer:


Multivariate feature imputation
===============================

A more sophisticated approach is to use the :class:`IterativeImputer` class,
which models each feature with missing values as a function of other features,
and uses that estimate for imputation. It does so in an iterated round-robin
fashion: at each step, a feature column is designated as output `y` and the
other feature columns are treated as inputs `X`. A regressor is fit on `(X, y)`
for known `y`. Then, the regressor is used to predict the missing values of `y`.
This is done for each feature in an iterative fashion, and then is repeated for
`n_iter` imputation rounds. The results of the final imputation round are
returned.

    >>> import numpy as np
    >>> from sklearn.impute import IterativeImputer
    >>> imp = IterativeImputer(n_iter=10, random_state=0)
    >>> imp.fit([[1, 2], [np.nan, 3], [7, np.nan]])
    IterativeImputer(imputation_order='ascending', initial_strategy='mean',
             max_value=None, min_value=None, missing_values=nan, n_iter=10,
             n_nearest_features=None, predictor=None, random_state=0,
             sample_posterior=False, verbose=False)
    >>> X_test = [[np.nan, 2], [6, np.nan], [np.nan, 6]]
    >>> print(np.round(imp.transform(X_test)))
    [[ 1.  2.]
     [ 6.  3.]
     [24.  6.]]

Both :class:`SimpleImputer` and :class:`IterativeImputer` can be used in a Pipeline
as a way to build a composite estimator that supports imputation.
See :ref:`sphx_glr_auto_examples_plot_missing_values.py`.


.. _multiple_imputation:

Multiple vs. Single Imputation
==============================

In the statistics community, it is common practice to perform multiple imputations,
generating, for example, `m` separate imputations for a single feature matrix.
Each of these `m` imputations is then put through the subsequent analysis pipeline
(e.g. feature engineering, clustering, regression, classification). The `m` final
analysis results (e.g. held-out validation errors) allow the data scientist
to obtain understanding of how analytic results may differ as a consequence
of the inherent uncertainty caused by the missing values. The above practice
is called multiple imputation.

Our implementation of :class:`IterativeImputer` was inspired by the R MICE
package (Multivariate Imputation by Chained Equations) [1], but differs from
it in setting single imputation to default instead of multiple imputation.
However, :class:`IterativeImputer` can also be used for multiple imputations
by applying it repeatedly to the same dataset with different random seeds when
``sample_posterior=True``.

See Chapter 4 of "Statistical Analysis with Missing Data" [2] for more discussion
on multiple vs. single imputations.

It is still an open problem as to how useful single vs. multiple imputation is in
the context of prediction and classification when the user is not interested in
measuring uncertainty due to missing values.

Note that a call to the ``transform`` method of :class:`IterativeImputer` is not
allowed to change the number of samples. Therefore multiple imputations cannot be
achieved by a single call to ``transform``.


References
==========
[1] Stef van Buuren, Karin Groothuis-Oudshoorn (2011). "mice: Multivariate Imputation
    by Chained Equations in R". Journal of Statistical Software 45: 1-67.

[2] Roderick J A Little and Donald B Rubin (1986). "Statistical Analysis with Missing
    Data". John Wiley & Sons, Inc., New York, NY, USA.
