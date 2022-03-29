.. _array_api:

================================
Array API support (experimental)
================================

.. currentmodule:: sklearn

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

    >>> X_np, y_np = make_classification(random_state=0)
    >>> X_cu = xp.asarray(X_np)
    >>> y_cu = xp.asarray(y_np)
    >>> X_cu.device
    <CUDA Device 0>

    >>> with config_context(array_api_dispatch=True):
    ...     lda = LinearDiscriminantAnalysis()
    ...     X_trans = lda.fit_transform(X_cu, y_cu)
    >>> X_trans.device
    <CUDA Device 0>

After the model is trained, fitted attributes that are arrays will also be
from the same Array API namespace as the training data. For example, if CuPy's
Array API namespace was used for training, then fitted attributes will be on the
GPU. Moving the fitted attributes to the CPU depends on the array library. For
example, if you use `cupy.array_api` then the transfer can be done as follows::

.. code-block:: python

    >>> for key, attribute in vars(lda).items():
    ...     if not hasattr(attribute, "__array_namespace__"):
    ...         # Does not implement the Array API specification
    ...         continue
    ...     # Convert `cupy.array_api` array into a NumPy ndarray
    ...     setattr(lda, key, attribute._array.get())
    >>> # transform on a CPU
    >>> X_trans = lda.transform(X_np)
    >>> type(X_trans)
    <class 'numpy.ndarray'>

.. _array_api_estimators:

Estimators with support for `Array API`-compatible inputs
=========================================================

- :class:`discriminant_analysis.LinearDiscriminantAnalysis` (with `solver="svd"`)

Coverage for more estimators is expected to grow over time. Please follow the
dedicated `meta-issue on GitHub
<https://github.com/scikit-learn/scikit-learn/issues/22352>`_ to track progress.
