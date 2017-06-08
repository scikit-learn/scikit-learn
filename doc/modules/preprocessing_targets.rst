.. currentmodule:: sklearn.preprocessing

.. _preprocessing_targets:

==========================================
Transforming the prediction target (``y``)
==========================================

Transforming target in regression
---------------------------------

:class:`TransformedTargetRegressor` is a meta-estimator transforming the target
before fitting a regression model and inverting back the prediction to the
original space. The meta-estimator take a regressor and a transformer as
argument::

  >>> import numpy as np
  >>> from sklearn.datasets import load_diabetes
  >>> from sklearn import preprocessing
  >>> from sklearn.linear_model import LinearRegression
  >>> from sklearn.model_selection import train_test_split
  >>> diabetes = load_diabetes()
  >>> X = diabetes.data
  >>> y = diabetes.target
  >>> transformer = preprocessing.StandardScaler()
  >>> regressor = LinearRegression()
  >>> clf = preprocessing.TransformedTargetRegressor(regressor=regressor,
  ...                                                transformer=transformer)
  >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
  >>> clf.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score:', clf.score(X_test, y_test))
  R2 score: 0.359400909897

The transformer can also be replaced by a function and an inverse function when
creating the meta-estimator. We can define the following two functions::

  >>> from __future__ import division
  >>> def func(x):
  ...     return x * 2
  >>> def inverse_func(x):
  ...     return x / 2

Subsequently, the meta-classifier is created as::

  >>> clf = preprocessing.TransformedTargetRegressor(regressor=regressor,
  ...                                                func=func,
  ...                                                inverse_func=inverse_func)
  >>> clf.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score:', clf.score(X_test, y_test))
  R2 score: 0.359400909897

By default, the provided function are checked to be the inverse of each
other. However, it is possible to bypass this checking by setting
``check_inverse`` to ``False``::

  >>> def inverse_func(x):
  ...     return x + 2
  >>> clf = preprocessing.TransformedTargetRegressor(regressor=regressor,
  ...                                                func=func,
  ...                                                inverse_func=inverse_func,
  ...                                                check_inverse=False)
  >>> clf.fit(X_train, y_train) # doctest: +ELLIPSIS
  TransformedTargetRegressor(...)
  >>> print('R2 score:', clf.score(X_test, y_test))
  R2 score: -5.63480137539

.. note::

   Either ``transformer`` or ``func`` and ``inverse_func`` cannot be set at the
   same time. An error will be raised otherwise.

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
