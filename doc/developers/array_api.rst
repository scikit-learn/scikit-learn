.. _array_api_developers:

========================================
Adding array API support to an estimator
========================================

.. currentmodule:: sklearn

This section is intended for contributors who add or maintain array API support
in scikit-learn estimators. For how to use array API support, see the
:ref:`user guide <array_api>`.

The helpers described here are private (they live in
``sklearn/utils/_array_api.py``) and can change without deprecation.

Overview
========

An estimator that supports the array API:

- runs the computations of ``fit`` in the array library (namespace) and on the
  device of ``X``;
- stores its fitted attributes as NumPy arrays;
- runs ``predict``, ``transform`` and other inference methods in the namespace
  and on the device of the ``X`` they are called with, and returns arrays from
  that namespace on that device.

See :ref:`array_api_fitted_attributes` for this contract from the point of view
of users.

.. TODO(#34604): remove this note once `_NUMPY_FITTED_ATTRS` is removed.

.. note::
   Storing fitted attributes as NumPy arrays is being rolled out estimator by
   estimator. Only estimators whose class name is listed in
   ``_NUMPY_FITTED_ATTRS`` follow the contract described here. The other
   estimators store their fitted attributes in the namespace of ``X`` and call
   ``check_same_namespace`` in their inference methods.

Fitting
=======

Fitted attributes are converted to NumPy by the ``_fit_context`` decorator
(``sklearn/base.py``) after the decorated method returns. Methods decorated with
it, usually ``fit`` and ``partial_fit``, therefore do not need any extra code,
and every ``return`` statement of such a method is covered. A method that sets
fitted attributes without being decorated, or without calling a decorated
method, does not get its attributes converted.

The conversion only looks at the estimator's own attributes and the containers
(``dict``, ``list``, ``tuple``, ...) stored in them. Nested estimators are left
alone, because each of them converts its own attributes when its own decorated
fit method returns.

An attribute that holds arrays but is neither an array, a container nor an
estimator can define ``__sklearn_array_api_convert__(self, converter)``, which
returns a copy of the object with ``converter`` applied to each of its arrays.

Inference methods
=================

In inference methods, access fitted attributes through ``_fitted_attrs_like``,
which returns a read-only view of the estimator. Reading an attribute from the
view moves it to the namespace and device of ``X``::

    from sklearn.utils._array_api import _fitted_attrs_like

    def decision_function(self, X):
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)
        fitted = _fitted_attrs_like(self, X)
        return X @ fitted.coef_.T + fitted.intercept_

If the namespace and device of ``X`` are already known, pass them instead of
``X`` to avoid inspecting ``X`` again::

    xp, _, device = get_namespace_and_device(X)
    fitted = _fitted_attrs_like(self, xp=xp, device=device)

Keep the following in mind:

- Read attributes from ``self`` when only their metadata is needed (for example
  ``self.coef_.ndim`` or ``self.classes_.shape[0]``). Metadata does not depend
  on the namespace, and reading it from the view would copy the array for
  nothing.
- A moved attribute is cached in the view, so reading it a second time from the
  same view does not copy it again. The cache lives as long as the view. An
  inference method that calls another inference method (for example ``predict``
  calling ``decision_function``) copies the attributes that both of them use
  twice.
- The view cannot be assigned to. Assign to the estimator instead.

State that must follow ``X``
============================

Fitted state that is read back during fitting is a NumPy array, even if the
previous fit was done with a different namespace. With ``warm_start=True`` or in
``partial_fit``, move the stored state to the namespace and device used for the
computation before using it, for example with
``move_to(self.coef_, xp=xp, device=device)``.

``classes_`` is a NumPy array, which allows it to hold string class labels that
cannot be represented in every namespace. To turn scores into class labels,
compute the indices in the namespace of ``X``, move them to NumPy, take the
labels from ``classes_`` and move numerical labels back to the namespace and
device of ``X``. String labels are returned as a NumPy array. See
``LinearClassifierMixin.predict`` in ``sklearn/linear_model/_base.py``.

Code that combines ``classes_`` with arrays that follow ``X``, for example
predictions or ``labels=estimator.classes_`` passed to a metric, has to convert
one of them explicitly.

Meta-estimators and nested fits
===============================

.. TODO: update this section once the conversion of nested fits is settled,
   see the discussion in https://github.com/scikit-learn/scikit-learn/pull/35053

The attributes of a sub-estimator are converted to NumPy as soon as the
sub-estimator's own fit method returns. This also happens while the fit method
of the meta-estimator is still running. A meta-estimator that uses the fitted
attributes of its sub-estimators during ``fit`` therefore sees NumPy arrays,
even when ``X`` is from another namespace, and has to move them (or use
``_fitted_attrs_like``) before combining them with ``X``.

Testing
=======

Add the ``array_api_support`` tag to an estimator's set of tags to indicate that
it supports the array API. This enables dedicated checks as part of the common
tests. They verify that:

- the fitted attributes and the outputs of the inference methods are the same
  when using NumPy and array API inputs, including inputs where ``y`` and
  ``sample_weight`` use a different namespace or device than ``X``, and string
  ``y``;
- all fitted attributes that are arrays are NumPy arrays;
- the outputs of the inference methods are in the namespace and on the device of
  their input;
- an estimator fitted with one namespace can be used for inference with another,
  and the other way around (``check_array_api_cross_namespace_inference``).

To run these checks you need to install
`array-api-strict <https://data-apis.org/array-api-strict/>`__ in your
test environment. This allows you to run checks without having a
GPU. To run checks on real GPU devices you also need to install
`PyTorch <https://pytorch.org/>`__, `CuPy <https://cupy.dev/>`__, and/or
`dpnp <https://intelpython.github.io/dpnp/>`__, and have compatible GPU
hardware. Full GPU coverage is expected to be split across machines because
CUDA, MPS, and Intel GPU backends require different hardware. Checks that can
not be executed or have missing dependencies will be automatically skipped.
Therefore it's important to run the tests with the `-v` flag to see which
checks are skipped:

.. code-block:: console

    $ pip install array-api-strict  # and other libraries as needed
    $ pytest -k "array_api" -v

On macOS, set ``PYTORCH_ENABLE_MPS_FALLBACK=1`` to run the checks on the MPS
device (see :ref:`mps_support`).

Running the scikit-learn tests against `array-api-strict` should help reveal
most code problems related to handling multiple device inputs via the use of
simulated non-CPU devices. This allows for fast iterative development and debugging of
array API related code.

However, to ensure full handling of PyTorch, CuPy, or dpnp inputs allocated on
actual GPU devices, it is necessary to run the tests against those libraries and
hardware. This can either be achieved by using
`Google Colab <https://gist.github.com/EdAbati/ff3bdc06bafeb92452b3740686cc8d7c>`__
for CUDA or leveraging our CI infrastructure on pull requests. CUDA and Intel
GPU tests are manually triggered by maintainers. Intel GPU testing for PyTorch
XPU and dpnp is run on a dedicated self-hosted runner:
`probabl-ai/scikit-learn-intel-workflow
<https://github.com/probabl-ai/scikit-learn-intel-workflow>`__.
