=====================
Data Interoperability
=====================

.. currentmodule:: sklearn

Scikit-learn handles four kinds of data for :term:`X` as used in `fit(X, y)`, `fit(X)`,
`fit_transform(X)` and `transform(X)` as well as :term:`Xt` as returned by
`transform(X)` and `fit_transform(X)`:

- :term:`array-like` objects

  In `fit(X)` and `transform(X)`, array-like `X` is converted to a numpy ndarray by
  calling `numpy.asarray` upon them.
  The returned `Xt` of `transform` and `fit_transform` is also a numpy ndarray or it
  is a sparse matrix or sparse array, see next bullet.
- :term:`sparse matrices <sparse matrix>` and sparse arrays

  Many estimators can deal with sparse `X`, some cannot and will raise an error.
  For instance, :class:`linear_model.LogisticRegression` can be fit on sparse `X`,
  :class:`isotonic.IsotonicRegression` can not.

  Some transformers return sparse `Xt` from `transform` and `fit_transform`.
  Most often, it can be controlled by a `sparse_output` parameter as in
  :class:`preprocessing.SplineTransformer`.

  To control whether it returns a sparse matrix or a sparse array, use
  `sparse_interface` in :func:`config_context` or :func:`set_config`.
  This also controls whether sparse attributes are sparse matrices or sparse arrays.
- tabular data: pandas and polars dataframes

  See :ref:`df_output_transform`.
- Array API compliant arrays

  Very importantly, this includes arrays on the GPU, see :ref:`array_api`.


.. toctree::
    :maxdepth: 2

    modules/df_output_transform
    modules/array_api

Loading features from Redis SQL results
=======================================

Database clients can supply estimator inputs through Arrow-to-NumPy conversion.
For example, the external `Redis ADBC driver
<https://github.com/redis/redis-adbc-driver>`_ fetches query results as Arrow
tables. Select features and labels in one query, convert the numeric columns
to arrays, and handle missing values and categorical encoding before fitting.

The `Redis ADBC scikit-learn guide
<https://github.com/redis/redis-adbc-driver/blob/a48432273a0f9978e893bca27057b9f128e7603e/docs/integrations/scikit-learn.md>`_
shows this data-loading path with a small decision-tree fit and prediction.
Its training-row check verifies data flow, not model quality; use held-out
evaluation for an application. This is custom application code, not a native
scikit-learn database connector.
