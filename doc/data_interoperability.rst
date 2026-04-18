=====================
Data Interoperability
=====================

.. currentmodule:: sklearn

Scikit-learn handles four kinds of data for :term:`X` as in `fit(X)` and
`fit(X, y)` as well as :term:`Xt` as returned by `transform(X)`:

- :term:`array-like` objects

  They are converted to a numpy ndarray by calling `numpy.asarray` upon them.
  The returned `Xt` of `transform` and `fit_transform` is also a numpy ndarray.
- sparse matrices and sparse arrays

  To control whether `transform` returns a sparse matrix or a sparse array, use
  `sparse_interface` in :func:`config_context` or :func:`set_config`. This also
  controls whether sparse attributes are sparse matrices or sparse arrays.
  Note that it is estimator specific if a sparse object is returned. Most often, it
  can be controlled by `sparse_output` as in :class:`preprocessing.SplineTransformer`.
- tabular data: pandas and polars dataframes

  See :ref:`df_output_transform`.
- Array API compliant arrays

  Very importantly, this includes arrays on the GPU, see :ref:`array_api`.


.. toctree::
    :maxdepth: 2

    modules/df_output_transform
    modules/array_api
