.. currentmodule:: sklearn.preprocessing

.. _preprocessing_targets:

==========================================
Transforming the prediction target (``y``)
==========================================

Transforming target in regression
---------------------------------

:class:`TransformedTargetRegressor` transforms the targets ``y`` before fitting a
regression model. The predictions are mapped back to the original space via an
inverse transform. It takes as an argument the regressor that will be used for
prediction, and the transformer that will be applied to the target variable::

  >>> import numpy as np
  >>> from sklearn.datasets import load_boston
  >>> from sklearn.preprocessing import (TransformedTargetRegressor,
  ...                                    QuantileTransformer)
  >>> from sklearn.linear_model import LinearRegression
  >>> from sklearn.model_selection import train_test_split
  >>> boston = load_boston()
  >>> X = boston.data
  >>> y = boston.target
  >>> transformer = QuantileTransformer(output_distribution='normal')
  >>> regressor = LinearRegression()
  >>> regr = TransformedTargetRegressor(regressor=regressor,
  ...                                   transformer=transformer)
  >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
  >>> regr.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score: {0:.2f}'.format(regr.score(X_test, y_test)))
  R2 score: 0.67
  >>> raw_target_regr = LinearRegression().fit(X_train, y_train)
  >>> print('R2 score: {0:.2f}'.format(raw_target_regr.score(X_test, y_test)))
  R2 score: 0.64

For simple transformations, instead of a Transformer object, a pair of
functions can be passed, defining the transformation and its inverse mapping::

  >>> from __future__ import division
  >>> def func(x):
  ...     return np.log(x)
  >>> def inverse_func(x):
  ...     return np.exp(x)

Subsequently, the object is created as::

  >>> regr = TransformedTargetRegressor(regressor=regressor,
  ...                                   func=func,
  ...                                   inverse_func=inverse_func)
  >>> regr.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score: {0:.2f}'.format(regr.score(X_test, y_test)))
  R2 score: 0.65

By default, the provided functions are checked at each fit to be the inverse of
each other. However, it is possible to bypass this checking by setting
``check_inverse`` to ``False``::

  >>> def inverse_func(x):
  ...     return x
  >>> regr = TransformedTargetRegressor(regressor=regressor,
  ...                                   func=func,
  ...                                   inverse_func=inverse_func,
  ...                                   check_inverse=False)
  >>> regr.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score: {0:.2f}'.format(regr.score(X_test, y_test)))
  R2 score: -4.50

.. note::

   The transformation can be triggered by setting either ``transformer`` or the
   pair of functions ``func`` and ``inverse_func``. However, setting both
   options will raise an error.

Label binarization
------------------

:class:`LabelBinarizer` is a utility class to help create a label indicator
matrix from a list of multi-class labels::

    >>> from sklearn import preprocessing
    >>> lb = preprocessing.LabelBinarizer()
    >>> lb.fit([1, 2, 6, 4, 2])
    LabelBinarizer(neg_label=0, pos_label=1, sparse_output=False)
    >>> lb.classes_
    array([1, 2, 4, 6])
    >>> lb.transform([1, 6])
    array([[1, 0, 0, 0],
           [0, 0, 0, 1]])

For multiple labels per instance, use :class:`MultiLabelBinarizer`::

    >>> lb = preprocessing.MultiLabelBinarizer()
    >>> lb.fit_transform([(1, 2), (3,)])
    array([[1, 1, 0],
           [0, 0, 1]])
    >>> lb.classes_
    array([1, 2, 3])

Label encoding
--------------

:class:`LabelEncoder` is a utility class to help normalize labels such that
they contain only values between 0 and n_classes-1. This is sometimes useful
for writing efficient Cython routines. :class:`LabelEncoder` can be used as
follows::

    >>> from sklearn import preprocessing
    >>> le = preprocessing.LabelEncoder()
    >>> le.fit([1, 2, 2, 6])
    LabelEncoder()
    >>> le.classes_
    array([1, 2, 6])
    >>> le.transform([1, 1, 2, 6])
    array([0, 0, 1, 2])
    >>> le.inverse_transform([0, 0, 1, 2])
    array([1, 1, 2, 6])

It can also be used to transform non-numerical labels (as long as they are
hashable and comparable) to numerical labels::

    >>> le = preprocessing.LabelEncoder()
    >>> le.fit(["paris", "paris", "tokyo", "amsterdam"])
    LabelEncoder()
    >>> list(le.classes_)
    ['amsterdam', 'paris', 'tokyo']
    >>> le.transform(["tokyo", "tokyo", "paris"])
    array([2, 2, 1])
    >>> list(le.inverse_transform([2, 2, 1]))
    ['tokyo', 'tokyo', 'paris']
