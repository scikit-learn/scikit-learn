
.. _pandas:

==================================================================
The basics of using pandas and scikit-learn
==================================================================

In general, scikit-learn supports the use of pandas DataFrames somewhat implicitly. However, this implicit support has its conditions and potential pitfalls and some (not all) of these will be outlined below.

Every scikit-learn estimator supports the use of DataFrames which is achieved by obtaining a numpy array using the `.values` property of the DataFrame class. As mentioned above, there are several conditions on such an input DataFrame one of which is that the data in the DataFrame columns used with the estimator are of type numeric. This can be checked e.g. using something like:

  >>> import pandas as pd
  >>> from pandas.api.types import is_numeric_dtype
  >>> df = pd.DataFrame({'foo': [1,2,3],
                         'bar': [1.2, 2.3, 3.]})
  >>> all([is_numeric_dtype(df[x]) for x in df.columns]) # should return True
  >>> df = pd.DataFrame({'foo': [1,2,3],
                         'bar': [1.2, 2.3, 3.],
                         'baz': ['foo', 'bar', 'baz']})
  >>> all([is_numeric_dtype(df[x]) for x in df.columns]) # should return False

There are also a variety of classes such as pipelines or model selection that will pass a DataFrame along "as is" to the nested estimators. However, this is not guaranteed and some of the issues with this are outlined below and sources (where available) of discussions and work-arounds (the latter is provided without guarantee that it will work and care should be taken) are provided. It should also be mentioned that the `.values` property will create an in memory copy of the DataFrame thus using a numpy array in the first place can be more memory efficient as well as avoiding some of the potential pitfalls when using DataFrames.

Pandas in Pandas out
====================

It is currently not the case, that all estimators in scikit-learn return a DataFrame as output when provided with a DataFrame as an input. There is an ongoing discussion and drive to improve support on this matter although there are also voices that suggest that in general, scikit-learn should only be expected to work smoothly with numpy arrays and have some basic support for DataFrames.

Independent of that discussion, at the moment it is not guaranteed that scikit-learn operators with `.fit`, `.transform` (and `.predict`) capability support pandas in pandas out. However, there are ways around this such as an example given here(https://github.com/scikit-learn/scikit-learn/issues/5523#issuecomment-171674105) show, where adding additional functionality to the StandardScaler class adds the pandas in pandas out capability. Care should be taken that this does not take care of the column ordering problem that is discussed in the next section.

The column ordering problem
===========================

Because scikit-learn transforms DataFrames to numpy arrays, it should be assumed, that all information and benefits of column names is lost and that from that point forward, only column order and not column labels stay relevant. This can cause problems when e.g. pickling a trained estimator and later applying it to a new DataFrame that while having the same data columns and labels has those in a different order compared to the original DataFrame. Intuitively it might be assumed that because scikit-learn handles the use of DataFrames so smoothly, the same goes for re-ordering labeled DataFrames but this is **not** the case.

An example of how this might impact your future prediction can be seen in the example given here(https://github.com/scikit-learn/scikit-learn/issues/7242#issue-173131995).

At this stage, it is the users responsibility, to ensure that the column ordering in the data used for training the estimator is the same as the ordering of the data used for prediction. There is an ongoing discussion whether or not this will change in the future and this issue (https://github.com/scikit-learn/scikit-learn/issues/7242) should be watched and used to update this paragraph in the future. A simple and straight-forward way of ensuring that column ordering and column labels are the same is using something like `df.loc[:, list of column names]` to enforce the correct ordering.

Handling Categorical data
=========================

Section to be extended.

See the following references to get started:

- https://scikit-learn.org/stable/glossary.html#term-categorical-feature
- https://scikit-learn.org/stable/modules/preprocessing.html#preprocessing-categorical-features
- https://github.com/scikit-learn-contrib/sklearn-pandas


Dealing with heterogenous data
==============================

Many modern datasets used with Scikit-learn contain heterogenous data. For the purpose of adding bespoke preprocessing steps for separate columns, Scikit-learn provides an experimental ColumnTransformer API (https://scikit-learn.org/stable/modules/compose.html#column-transformer). This API (which might change in the future) allows the definition of different transformation steps to be applied to different columns in either arrays, sparse matrices or pandas DataFrames.

Dealing with missing values
===========================

As per the glosary, most scikit-learn estimators do not work with missing values. If they do, NaN is the preferred representation of missing values. For more details, see https://scikit-learn.org/stable/glossary.html#term-missing-values.
