.. Places parent toc into the sidebar

:parenttoc: True

.. include:: includes/big_toc_css.rst

.. _array_api:

================================
Array API support (experimental)
================================

The `Array API <https://data-apis.org/array-api/latest/>`_ specification defines
a standard API for all array manipulation libraries with a NumPy-like API.

Some scikit-learn estimators that primarily rely on NumPy (as opposed to using
Cython) to implement the algorithmic logic of their `fit`, `predict` or
`transform` methods can be configured to accept any Array API compatible input
datastructures and automatically dispatch operations to the underlying namespace
instead of relying on NumPy.

At this stage, this support is **considered experimental** and must be enabled
explicitly as explained in the following.

Example usage
=============

Here is an example code snippet to demonstrate how to use `CuPy
<https://cupy.dev/>`_ to run
:class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis` on a GPU::

    >>> from sklearn.datasets import make_classification
    >>> from sklearn import config_context
    >>> from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    >>> import cupy.array_api as xp

    >>> X, y = make_classification(random_state=0)
    >>> X = xp.asarray(X)
    >>> y = xp.asarray(y)
    >>> X.device
    <CUDA Device 0>

    >>> with config_context(array_api_dispatch=True):
    ...     lda = LinearDiscriminantAnalysis()
    ...     lda.fit(X, y)
    ...     X_trans = lda.transform(X)
    ...
    >>> X_trans.device
    <CUDA Device 0>


.. _array_api_estimators:

Estimators with support for `Array API`-compatible inputs
=========================================================

- :class:`sklearn.discriminant_analysis.LinearDiscriminantAnalysis` (with `solver="svd"`)

Coverage for more estimators is expected to grow over time. Please follow the
dedicated `meta-issue on GitHub
<https://github.com/scikit-learn/scikit-learn/issues/22352>`_ to track progress.
