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

Imputer transformers can be used in a Pipeline as a way to build a composite
estimator that supports imputation.  See
:ref:`sphx_glr_auto_examples_plot_missing_values.py`.


Simple univariate imputation
============================

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

.. _knnimpute:

Nearest neighbors imputation
===============================

The :class:`KNNImputer` class provides imputation for completing missing
values using the k-Nearest Neighbors approach. Each sample's missing values
are imputed using values from ``n_neighbors`` nearest neighbors found in the
training set. Note that if a sample has more than one feature missing, then
the sample can potentially have multiple sets of ``n_neighbors``
donors depending on the particular feature being imputed.

Each missing feature is then imputed as the average, either weighted or
unweighted, of these neighbors. Where the number of donor neighbors is less
than ``n_neighbors``, the training set average for that feature is used
for imputation. The total number of samples in the training set is, of course,
always greater than or equal to the number of nearest neighbors available for
imputation, depending on both the overall sample size as well as the number of
samples excluded from nearest neighbor calculation because of too many missing
features (as controlled by ``row_max_missing``).
For more information on the methodology, see ref. [#]_.

The following snippet demonstrates how to replace missing values,
encoded as ``np.nan``, using the mean feature value of the two nearest
neighbors of the rows that contain the missing values::

    >>> import numpy as np
    >>> from sklearn.impute import KNNImputer
    >>> nan = np.nan
    >>> X = [[1, 2, nan], [3, 4, 3], [nan, 6, 5], [8, 8, 7]]
    >>> imputer = KNNImputer(n_neighbors=2, weights="uniform")
    >>> imputer.fit_transform(X)
    array([[1. , 2. , 4. ],
           [3. , 4. , 3. ],
           [5.5, 6. , 5. ],
           [8. , 8. , 7. ]])

.. [#] Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor
    Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing value
    estimation methods for DNA microarrays, BIOINFORMATICS Vol. 17 no. 6, 2001
    Pages 520-525.

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

When using the :class:`MissingIndicator` in a :class:`Pipeline`, be sure to use
the :class:`FeatureUnion` or :class:`ColumnTransformer` to add the indicator
features to the regular features. First we obtain the `iris` dataset, and add
some missing values to it.

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.impute import SimpleImputer, MissingIndicator
  >>> from sklearn.model_selection import train_test_split
  >>> from sklearn.pipeline import FeatureUnion, make_pipeline
  >>> from sklearn.tree import DecisionTreeClassifier
  >>> X, y = load_iris(return_X_y=True)
  >>> mask = np.random.randint(0, 2, size=X.shape).astype(np.bool)
  >>> X[mask] = np.nan
  >>> X_train, X_test, y_train, _ = train_test_split(X, y, test_size=100,
  ...                                                random_state=0)

Now we create a :class:`FeatureUnion`. All features will be imputed using
:class:`SimpleImputer`, in order to enable classifiers to work with this data.
Additionally, it adds the the indicator variables from
:class:`MissingIndicator`.

  >>> transformer = FeatureUnion(
  ...     transformer_list=[
  ...         ('features', SimpleImputer(strategy='mean')),
  ...         ('indicators', MissingIndicator())])
  >>> transformer = transformer.fit(X_train, y_train)
  >>> results = transformer.transform(X_test)
  >>> results.shape
  (100, 8)

Of course, we cannot use the transformer to make any predictions. We should
wrap this in a :class:`Pipeline` with a classifier (e.g., a
:class:`DecisionTreeClassifier`) to be able to make predictions.

  >>> clf = make_pipeline(transformer, DecisionTreeClassifier())
  >>> clf = clf.fit(X_train, y_train)
  >>> results = clf.predict(X_test)
  >>> results.shape
  (100,)

