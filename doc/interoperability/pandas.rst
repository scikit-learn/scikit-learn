
.. _pandas:

=======================
Pandas Interoperability
=======================

This part of the User Guide aims to highlight some of the particularities that
arise when using pandas datastructures as input for scikit-learn.

The basics of using pandas and Scikit-learn
==================================================================

Scikit-learn supports the use of
`pandas DataFrames <http://pandas.pydata.org/pandas-docs/stable/>`__
implicitly. However, this implicit support has its conditions and potential
pitfalls and some (not all) of these will be outlined below.

Reasons for these conditions and the difficulty of supporting pandas, stems
from the fact that pandas data structures such as Series and DataFrames can
hold heterogenous datatypes (different columns can contain different
datatypes). Therefore, Series and DataFrames can be thought of as containers
for arrays that hold the actual data in an homogenous format. Some of these
underlying data structures might not necessarily be representable as NumPy
arrays (e.g. Categorical data). Pandas supports these alternate data types using
pandas specific or 3rd party libraries that extend the NumPy type system. The
difficulty of supporting pandas data structure thus comes from pandas
capability of supporting heterogenous data in one data structure.

Every Scikit-learn estimator/transformer/pipeline
(for the rest of this section we shall call these primitives)
supports the use of DataFrames as inputs. Most do this by obtaining a
`NumPy array <https://docs.scipy.org/doc/numpy/user/>`__ using
the :meth:`~numpy.asarray` on a DataFrame object. The only exception where a
a DataFrame is being used explicitly is the
:class:`~sklearn.compose.ColumnTransformer` which is briefly
discussed below `Dealing with heterogenous data`_.

.. note::
  Starting with pandas version 0.24.0, it is encouraged to obtain
  NumPy arrays from DataFrames or Series via :meth:`.to_numpy()` instead of
  using :meth:`.values`. More details on this can be found in the
  `release notes <http://pandas-docs.github.io/pandas-docs-travis/whatsnew/v0.24.0.html#accessing-the-values-in-a-series-or-index>`__
  and the documentation `here <http://pandas.pydata.org/pandas-docs/stable/getting_started/basics.html#basics-dtypes>`__
  and `here <http://pandas.pydata.org/pandas-docs/stable/getting_started/basics.html#attributes-and-underlying-data>`__.

There are several conditions on using a DataFrame as an input to
most Scikit-learn estimators, one of which is that the data in the
DataFrame columns used by the estimator are of numerical type. Other conditions
and pitfalls are described in subsequent sections. The numerical condition can
be checked e.g. using something like::

  >>> import pandas as pd
  >>> from pandas.api.types import is_numeric_dtype
  >>> df = pd.DataFrame({'foo': [1, 2, 3],
  ...                    'bar': [1.2, 2.3, 3.]})
  >>> all([is_numeric_dtype(df[x]) for x in df.columns])
  True
  >>> df = pd.DataFrame({'foo': [1,2,3],
  ...                    'bar': [1.2, 2.3, 3.],
  ...                    'baz': ['foo', 'bar', 'baz']})
  >>> all([is_numeric_dtype(df[x]) for x in df.columns])
  False

There are also a variety of classes such as pipelines and model selection
processes that will pass a DataFrame along "as is" to the nested estimators
(see the next section for an example). However, this is not guaranteed and some
of the issues arising as a result are outlined below and sources
(where available) of
discussions and work-arounds (the latter is provided without guarantee that the
workaround will still work) are provided. It should also be mentioned that if
the DataFrame contains heterogenous data, the :meth:`~numpy.asarray()` function will
create an in-memory copy of the DataFrame, thus using a NumPy array in the
first place can be more memory efficient as well as avoiding some of the
potential pitfalls when using DataFrames. If the DataFrame contains only
homogenous data in the first place, no in-memory copy will be created using
:meth:`~numpy.asarray()`.

Pandas in **not** Pandas out
============================

Some primitives in Scikit-learn support the pandas in pandas out, however it
should generally be assumed that you get a Numpy array as an output when
providing a DataFrame as an input.

Example for Pandas working with Scikit-learn primitive::

  >>> import numpy as np
  >>> from sklearn.datasets import load_iris
  >>> import pandas as pd
  >>>
  >>> # make a dataframe
  >>> iris = load_iris()
  >>> X, y = iris.data[:-1,:], iris.target[:-1]
  >>> iris_pd = pd.DataFrame(X)
  >>> iris_pd.columns = iris.feature_names
  >>> iris_pd['target'] = y
  >>>
  >>> from sklearn.model_selection import train_test_split
  >>> train, test = train_test_split(iris_pd.copy(), test_size= 0.3)
  >>>
  >>> type(train)
  <class 'pandas.core.frame.DataFrame'>

However, using the :class:`~sklearn.preprocessing.StandardScaler` returns a
NumPy array instead of a DataFrame even though we use a DataFrame as input::

  >>> from sklearn.preprocessing import StandardScaler
  >>>
  >>> scaler = StandardScaler()
  >>> X = scaler.fit_transform(train)
  >>> type(X)
  <class 'numpy.ndarray'>

As this example shows, at the moment it is not guaranteed that Scikit-learn
primitivies with :meth:`.fit`, :meth:`.transform` (and :meth:`.predict`)
capability support pandas in pandas out. However, there are ways around this
such as an example given
`here <https://github.com/scikit-learn/scikit-learn/issues/5523#issuecomment-171674105>`__
show, where adding additional functionality to the StandardScaler class adds
the pandas in pandas out capability. Care should be taken as this does not
take care of the column ordering problem that is discussed in the next section.

The column ordering problem
===========================

Because Scikit-learn transforms DataFrames to NumPy arrays, it should be
assumed, that all information and benefits of column names is lost and that
from that point forward, only column order and not column labels stay relevant.
This can cause problems in general when predicting unseen data using a previously
trained estimator and applying it to the new data as it does not matter
that the unseen/new data has the same data columns and labels, they still
**must** be provided in the correct order too.
Scikit-learn does not check that the column order is consistent nor does
it do any automatic re-ordering of DataFrame columns!

An example of how this might impact your future prediction can be seen in the
example given below::

  >>> from sklearn.datasets import load_iris
  >>> import pandas as pd
  >>>
  >>> # make a dataframe
  >>> iris = load_iris()
  >>> X, y = iris.data[:-1,:], iris.target[:-1]
  >>> iris_pd = pd.DataFrame(X)
  >>> iris_pd.columns = iris.feature_names
  >>> iris_pd['target'] = y
  >>>
  >>> from sklearn.model_selection import train_test_split
  >>> train, test = train_test_split(iris_pd, test_size= 0.3, random_state=42)
  >>>
  >>> feature_columns_train = ['sepal length (cm)','sepal width (cm)',
  ...                          'petal length (cm)','petal width (cm)']
  >>> # last two correct order
  >>> feature_columns_test = ['sepal length (cm)','sepal width (cm)',
  ...                         'petal width (cm)','petal length (cm)']
  >>> # last two switched order
  >>>
  >>> from sklearn.linear_model import LogisticRegression
  >>> lg = LogisticRegression(n_jobs=4, random_state=123, verbose=0,
  ...                         penalty='l2', C=1.0,
  ...                         solver='lbfgs', multi_class='auto')
  >>> lg.fit(train[feature_columns_train], train['target'])
  LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
                     intercept_scaling=1, l1_ratio=None, max_iter=100,
                     multi_class='auto', n_jobs=4, penalty='l2', random_state=123,
                     solver='lbfgs', tol=0.0001, verbose=0, warm_start=False)
  >>>
  >>> res1 = lg.predict(test[feature_columns_train])
  >>> res1[:5]
  array([1, 0, 2, 1, 1])
  >>> # result is actually
  >>> res2 = lg.predict(test[feature_columns_test])
  >>> res2[:5]
  array([0, 0, 2, 0, 0])


At the time of writing, it is the user's responsibility to ensure that the
column ordering in the data used for training the estimator is the same as the
ordering of the data used for prediction. There is an ongoing discussion
whether or not this will change in the future and this
`issue <https://github.com/scikit-learn/scikit-learn/issues/7242>`__ should be
watched and used to update this paragraph in the future. A simple and straight-
forward way of ensuring that column ordering and column labels are the same is
using something like `df[list of column names]` to enforce the
correct ordering.

Handling Categorical data
=========================

For a general guide on how to get started with categorical features please refer
to :term:`categocrical feature` and :ref:`preprocessing_categorical_features`.
It is worth noting that as of :ref:`changes_0_20_3`, both
:class:`~sklearn.preprocessing.OneHotEncoder` and
:class:`~sklearn.preprocessing.OrdinalEncoder`
support string or Categorical columns coming straight from pandas DataFrames.


Dealing with heterogenous data
==============================

Many modern datasets used with Scikit-learn contain heterogenous data. For the
purpose of adding bespoke preprocessing steps for separate columns, Scikit-
learn provides an experimental :class:`~sklearn.compose.ColumnTransformer` API
(:ref:`column_transformer`).
This API (which might change in the future) allows the definition of different
transformation steps to be applied to different columns in either arrays,
sparse matrices or pandas DataFrames.

Dealing with missing values
===========================

As per the glossary, most Scikit-learn primitives do not work with missing
values. If they do, NaN is the preferred representation of missing values. For
more details, see :term:`missing values`. Non-numeric data is now also supported
via the ``'most_frequent'`` or ``'constant'`` of the
:class:`~sklearn.impute.SimpleImputer` class. For details see :ref:`impute`.


Sparse DataFrames Handling
=============================

.. note::
  **Issue:**
  ``Sparse DataFrames`` are not automatically converted to ``scipy.sparse``
  matrices.

In general, Sparse data structures (i.e. DataFrames, Series, Arrays) are memory
optimised structures of their standard counterparts. They work on the principle
that they contain a lot of NaN, 0, or another repeating value (this can be
specified), and as such a lot of memory can be saved, which means one can
potentially work with datasets that would otherwise be too large to fit into
available memory. However one has to be careful they don't get converted into
the dense format by mistake.

In Pandas, the main sparse data structures is: :class:`~pandas.SparseArray`.
However, Scikit-learn does not support sparse Pandas structures and by default
they will be converted to dense numpy arrays. The best way to use sparse
arrays in Scikit-learn is to convert them manually to sparce Scipy matrices.
The methods: :meth:`~pandas.DataFrame.to_sparse(fill_value=0)` and
:meth:`~pandas.SparseDataFrame.to_dense()` can be
used to convert between normal and sparse data structures.
The :meth:`~pandas.SparseDataFrame.density` property can be called on the
sparse structures to report sparseness.

In scipy.sparse we have a number of various sparse matrix classes, Scikit-learn
mostly uses CSR and CSC formats.

Example Usage
-------------

  >>> import numpy as np
  >>> import pandas as pd
  >>> from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, issparse
  >>>
  >>> arr = np.random.random(size=(1000, 1000))
  >>> arr[arr < .9] = 0
  >>>
  >>> sparse_df = pd.SparseDataFrame(arr, default_fill_value=0)
  >>> coo = sparse_df.to_coo()
  >>> #or
  >>> coo = coo_matrix(sparse_df)
  >>>
  >>> csr = coo.tocsr()
  >>> csc = coo.tocsc()
  >>>
  >>> print('Confirm both are sparse:',
  ...       issparse(coo) == issparse(csr) == issparse(csc) == True)
  Confirm both are sparse: True
  >>> print('Confirm same amount of non-empty values:',
  ...       coo.nnz == csr.nnz == csc.nnz)
  Confirm same amount of non-empty values: True


The code above highlights the following three elements:

1) If your sparse value is not NaN then it is important to specify
*default_fill_value* property when creating your pandas DataFrame, otherwise no
space saving will occur. Check this using the
:attr:`~pandas.SparseDataFrame.density()` property, which
should be less than 100% if successful. When creating the scipy sparse matrix,
this *default_fill_value* will be used for use as the sparse value (nnz).

2) Either the :meth:`~pandas.SparseDataFrame.to_coo()` method on the pandas
SparseDataFrame, or :class:`~scipy.sparse.coo_matrix` constructor are
alternative ways you can convert to a scipy sparse datastructure.

3) It is generally better to convert from your pandas Dataframe first to a
:class:`~scipy.sparse.coo_matrix`, as this is far quicker to construct,
and from this to then convert to a Compressed Row
:class:`~scipy.sparse.csr_matrix`, or Compressed Column
:class:`~scipy.sparse.csc_matrix` sparse matrix using the
:meth:`~scipy.sparse.csc_matrix.tocsr()` or
:meth:`~scipy.sparse.csr_matrix.tocsc()` methods respectively.
