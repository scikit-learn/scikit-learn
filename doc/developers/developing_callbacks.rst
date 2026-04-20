.. _developing_callbacks:

====================
Developing callbacks
====================

.. currentmodule:: sklearn.callback

The callback protocol
---------------------

To be compatible with scikit-learn estimators, callbacks must implement the
:class:`FitCallback` `protocol
<https://typing.python.org/en/latest/spec/protocol.html>`__:

.. code-block:: python

    class FitCallback(Protocol):

        def setup(self, estimator, context) -> None: ...

        def on_fit_task_begin(
            self,
            estimator,
            context,
            *,
            X=None,
            y=None,
            metadata=None,
            fitted_estimator=None
        ) -> None: ...

        def on_fit_task_end(
            self,
            estimator,
            context,
            *,
            X=None,
            y=None,
            metadata=None,
            fitted_estimator=None
        ) -> bool: ...

        def teardown(self, estimator, context) -> None: ...

The methods of the protocol, referred to as callback hooks, will be called at specific
steps during the fitting process of the estimator the callback is registered on:

- :meth:`~FitCallback.setup` and :meth:`~FitCallback.teardown`

  These hooks are only called once, respectively at the start and end of `fit`. They
  take care of setting up and tearing down the callback, like allocating and freeing
  resources for instance.

- :meth:`~FitCallback.on_fit_task_begin` and :meth:`~FitCallback.on_fit_task_end`

  These hooks are called at the beginning and end of each :ref:`task
  <callback_task_definition>` during `fit`.

  In concrete implementations of callbacks, only the optional keyword-only arguments
  actually used by the hook should be explicitly declared in the hook signature. The
  presence of an argument in the signature signals that the hook requires that argument,
  which allows the callback framework to avoid computing values that are not used by any
  registered callback.

  .. warning::

      These arguments must be defined as **keyword only**. If the kwargs are not keyword
      only, the values will not be provided to the hooks.

  Even if requested, the optional arguments might or might not be provided by the
  estimator, depending on its ability to produce them at this task. Thus the
  implementation of the hooks should not expect to always receive a value for each of
  them and adapt their behavior accordingly.

  .. dropdown:: Interrupting `fit`

      The :meth:`~FitCallback.on_fit_task_end` hook returns a boolean, which when set to
      `True`, requests the estimator to stop the `fit` process at this task. Note that
      estimators that don't aim to be interruptable will ignore this request and
      continue with the next task.

All the hooks receive as mandatory arguments the estimator instance calling the callback
and the :class:`~CallbackContext` object holding the contextual information that allows
to identify the task that is being processed as public attributes. See
:class:`~CallbackContext` for more details.

.. dropdown:: The `estimator` argument

    The estimator instance received by the hooks as a mandatory argument is in the same
    state as it was when calling the hook during `fit`. Therefore it is not expected to
    be fully fitted (except for the :meth:`~FitCallback.teardown` hook). Callbacks
    should not rely on it to predict, transform, etc ... but rather use the
    `fitted_estimator` when available.

Auto-propagated callbacks
-------------------------

Auto-propagated callbacks, i.e. callbacks that are expected to be propagated from
meta-estimators to their sub-estimators, must implement the
:class:`~AutoPropagatedCallback` protocol, an extension of the :class:`~FitCallback`
protocol:

.. code-block:: python

    class AutoPropagatedCallback(FitCallback, Protocol):

        @property
        def max_propagation_depth(self) -> int | None: ...

By contrast with regular callbacks that are only invoked at the tasks of the estimator
on which they are registered, auto-propagated callbacks are invoked at the tasks of
all the estimators in estimator compositions, up to the maximum propagation depth. If
set to 0, the callback is not propagated to sub-estimators and only invoked at the tasks
of the top-level estimator. If set to `None`, the callback is propagated to
sub-estimators at all nesting levels.

Auto-propagated callbacks should be registered on the top-level estimator. If the
top-level estimator does not support callbacks, they can be registered on
sub-estimators and are expected to work, possibly not at full capacity.

.. note::

    The :meth:`~FitCallback.setup` and :meth:`~FitCallback.teardown` hooks of an
    auto-propagated callback are also called only once, at the beginning and end of the
    top-level estimator's `fit` method. They are not called for any of its
    sub-estimators.

Callback shared state
---------------------

Since the estimator on which the callback is registered may be cloned and fitted
multiple times in a meta-estimator, callbacks should behave as if the same callback
instance were registered on multiple estimators. Therefore, `setup` / `teardown` should
not reset the state of the callback, and `on_fit_task_begin` / `on_fit_task_end` should
accumulate data across all fits. Indeed resetting state in `setup`/`teardown` would drop
information collected from previous or concurrent fits.

.. TODO: document the callback Manager.

Minimal example
---------------

Here is an example implementation of a simple custom callback that prints a message
every time it is invoked:

.. code-block:: python

    class MyCallback:

        def setup(self, estimator, context):
            print(f"Setup hook is being called in the {context.task_name} task.")

        def teardown(self, estimator, context):
            print(f"Teardown hook is being called in the {context.task_name} task.")

        def on_fit_task_begin(self, estimator, context, *, X=None):
            msg = f"{context.task_name} task is starting."
            if X is not None:
                msg += f" With training data of shape {X.shape}."
            print(msg)

        def on_fit_task_end(
            self, estimator, context, *, X=None, y=None, fitted_estimator=None
        ):
            msg = f"{context.task_name} task is ending."
            mean_squared_error = ((y - fitted_estimator.predict(X))**2).mean()
            msg += f" With a mean squared error of {mean_squared_error}."
            print(msg)


.. TODO: add a link to a more detailed gallery example of a custom callback.
