
.. _pandas:

=======================
Pandas Interoperability
=======================

The basics of using pandas and Scikit-learn
==================================================================

In principle, Scikit-learn supports the use of
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

Every Scikit-learn estimator supports the use of DataFrames which is achieved
by obtaining a `NumPy array <https://docs.scipy.org/doc/numpy/user/>`__ using
the :meth:`.values` property of the DataFrame class.

.. note::
  Starting with pandas version 0.24.0, it is encouraged to obtain
  NumPy arrays from DataFrames or Series via :meth:`.to_numpy()` instead of
  using :meth:`.values`. More details on this can be found in the
  `release notes <http://pandas-docs.github.io/pandas-docs-travis/whatsnew/v0.24.0.html#accessing-the-values-in-a-series-or-index>`__
  and the documentation `here <http://pandas.pydata.org/pandas-docs/stable/basics.html#basics-dtypes>`__
  and `here <http://pandas.pydata.org/pandas-docs/stable/basics.html#attributes-and-underlying-data>`__.

There are several
conditions on such an input DataFrame one of which is that the data in the
DataFrame columns used by the estimator are of numerical type. Other conditions
and pitfalls are described in subsequent sections. The numerical condition can
be checked e.g. using something like:

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
of the issues with this are outlined below and sources (where available) of
discussions and work-arounds (the latter is provided without guarantee that the
workaround will still work) are provided. It should also be mentioned that if
the DataFrame contains heterogenous data, the :meth:`.values` property will
create an in-memory copy of the DataFrame, thus using a NumPy array in the
first place can be more memory efficient as well as avoiding some of the
potential pitfalls when using DataFrames. If the DataFrame contains only
homogenous data in the first place, no in-memory copy will be created using
:meth:`.values`.

Pandas in Pandas out
====================

It is currently not the case, that all estimators in scikit-learn return a
DataFrame as an output when provided with a DataFrame as an input. There is an
ongoing discussion and drive to improve support on this matter although there
are also voices that suggest that in general, scikit-learn should only be
expected to work smoothly with numpy arrays and have some basic support for
DataFrames.

Example for pandas working with scikit-learn estimator/transformer:

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

However, removing some random values from the dataset and using the
:class:`~sklearn.impute.SimpleImputer` to replace the NaNs returns a numpy
array instead of a DataFrame even though we use a DataFrame as input.

>>> rng = np.random.RandomState(42)
>>> # selecting some random indices to replace
>>> idx = train.index[rng.binomial(1, 0.2, train.shape[0]).astype(bool)]
>>> train.loc[idx, 'sepal length (cm)'] = np.nan
>>>
>>> from sklearn.impute import SimpleImputer
>>>
>>> imputer = SimpleImputer()
>>> X = imputer.fit_transform(train)
>>> type(X)
<class 'numpy.ndarray'>

Independent of this, at the moment it is not guaranteed that scikit-learn
operators with :meth:`.fit`, :meth:`.transform` (and :meth:`.predict`)
capability support pandas in pandas out. However, there are ways around this
such as an example given
`here <https://github.com/scikit-learn/scikit-learn/issues/5523#issuecomment-171674105>`__
show, where adding additional functionality to the StandardScaler class adds
the pandas in pandas out capability. Care should be taken that this does not
take care of the column ordering problem that is discussed in the next section.

The column ordering problem
===========================

Because Scikit-learn transforms DataFrames to numpy arrays, it should be
assumed, that all information and benefits of column names is lost and that
from that point forward, only column order and not column labels stay relevant.
This can cause problems when e.g. pickling a trained estimator and later
applying it to a new DataFrame that, while having the same data columns and
labels, has those in a different order compared to the original DataFrame.
Intuitively it might be assumed that because Scikit-learn handles the use of
DataFrames so smoothly in most cases, the same goes for re-ordering labeled
DataFrames but this is **not** the case.

An example of how this might impact your future prediction can be seen in the
example given below (original with slight modifications adjusting for current
API, thanks to `SauceCat <https://github.com/scikit-learn/scikit-learn/issues/7242#issue-173131995>`__).

>>> # for simplification, consider a very simple case
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
          intercept_scaling=1, max_iter=100, multi_class='auto', n_jobs=4,
          penalty='l2', random_state=123, solver='lbfgs', tol=0.0001,
          verbose=0, warm_start=False)
>>>
>>> prob1 = lg.predict_proba(test[feature_columns_train])
>>> prob1[:5]
array([[4.09709461e-03, 8.21100411e-01, 1.74802495e-01],
       [9.42618164e-01, 5.73813720e-02, 4.64354721e-07],
       [2.72655051e-07, 5.28875447e-03, 9.94710973e-01],
       [6.86315850e-03, 7.80379358e-01, 2.12757484e-01],
       [1.64263139e-03, 7.43621534e-01, 2.54735834e-01]])
>>> # result is actually
>>> prob2 = lg.predict_proba(test[feature_columns_test])
>>> prob2[:5]
array([[7.92829716e-01, 1.79085973e-01, 2.80843105e-02],
       [9.95986933e-01, 4.01303839e-03, 2.87377384e-08],
       [2.47995509e-03, 7.79557758e-03, 9.89724467e-01],
       [7.09780229e-01, 2.39891794e-01, 5.03279763e-02],
       [5.62705633e-01, 3.48565301e-01, 8.87290655e-02]])


At the time of writing, it is the users responsibility to ensure that the
column ordering in the data used for training the estimator is the same as the
ordering of the data used for prediction. There is an ongoing discussion
whether or not this will change in the future and this
`issue <https://github.com/scikit-learn/scikit-learn/issues/7242>`__ should be
watched and used to update this paragraph in the future. A simple and straight-
forward way of ensuring that column ordering and column labels are the same is
using something like :meth:`df.loc[:, list of column names]` to enforce the
correct ordering.

Handling Categorical data
=========================

Section to be extended.

See the following references to get started:

- https://scikit-learn.org/stable/glossary.html#term-categorical-feature
- https://scikit-learn.org/stable/modules/preprocessing.html#preprocessing-categorical-features
- https://github.com/scikit-learn-contrib/sklearn-pandas


Dealing with heterogenous data
==============================

Many modern datasets used with Scikit-learn contain heterogenous data. For the
purpose of adding bespoke preprocessing steps for separate columns, Scikit-
learn provides an experimental :class:`~sklearn.compose.ColumnTransformer` API.
This API (which might change in the future) allows the definition of different
transformation steps to be applied to different columns in either arrays,
sparse matrices or pandas DataFrames.

Dealing with missing values
===========================

As per the glosary, most scikit-learn estimators do not work with missing
values. If they do, NaN is the preferred representation of missing values. For
more details, see https://scikit-learn.org/stable/glossary.html#term-missing-values.


Sparse DataFrames Handling
=============================

**Issue:**
``Sparse DataFrames`` are not automatically converted to ``scipy.sparse``
matrices.

This is an issue which has vastly improved from pandas version 0.21.1 onwards.
The conversation from dataframes has been largely optimized and are much faster
to convert.

In general, Sparse datastructures (i.e. DataFrames, Series, Arrays) are memory
optimised structures of their standard counterparts. They work on the principle
that they contain a lot of NaN, 0, or another repeating value (this can be
specified), and as such a lot of memory can be saved, which means one can
potentially work with datasets that would otherwise be too large to fit into
available memory. However one has to be careful they don't get converted into
the dense format by mistake.

In Pandas, the sparse datastructrures are: :class:`~pandas.SparseDataFrame`,
:class:`~pandas.SparseSeries` and :class:`~pandas.SparseArray`.
The methods: :meth:`.to_sparse(fill_value=0)` and :meth:`.to_dense()` can be
used to convert between normal and sparse data structures.
The `.density` property can be called on the sparse structures to report
sparseness.

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
  >>> print('Density: {:.2%}'.format(sparse_df.density))
  Density: 10.00%
  >>>
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
space saving will occur. Check this using the :meth:`.density` property, which
should be less than 100% if successful. When creating the scipy sparse matrix,
this *default_fill_value* will be used for use as the sparse value (nnz).

2) Either the :meth:`.to_coo()` method on the pandas dataframe, or
:meth:`coo_matrix()` constructor are alternative ways you can convert to a
scipy sparse datastructure.

3) It is generally better to convert from your pandas Dataframe first to a
:class:`coo_matrix`, as this is far quicker to construct, and from this to then
convert to a Compressed Row :class:`csr_matrix`, or Compressed Column
:class:`csc_matrix` sparse matrix using the :meth:`.tocsr()` or
:meth:`.tocsc()` methods respectively.
