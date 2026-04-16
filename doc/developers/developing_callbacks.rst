.. _developing_callbacks:

====================
Developing callbacks
====================

.. currentmodule:: sklearn.callback

This section is intended for developers who wish to implement callbacks compatible
with the scikit-learn callback API.

The callback protocol
---------------------

To be compatible with scikit-learn estimators, callbacks must implement the
:class:`FitCallback` `protocol
<https://typing.python.org/en/latest/spec/protocol.html>`__:

.. code-block:: python

    class FitCallback(Protocol):

        def setup(self, estimator, context): -> None

        def on_fit_task_begin(
            self,
            estimator,
            context,
            *,
            X=None,
            y=None,
            metadata=None,
            fitted_estimator=None
        ): -> None

        def on_fit_task_end(
            self,
            estimator,
            context,
            *,
            X=None,
            y=None,
            metadata=None,
            fitted_estimator=None
        ): -> bool

        def teardown(self, estimator, context): -> None

These methods, referred to as callback hooks, will be called at different specific steps
during the fitting process of the estimator the callback is registered on.

All the hooks receive as mandatory arguments the estimator instance calling the callback
(`estimator`), and the :class:`~CallbackContext` object (`context`) holding the
contextual information of the task representing the whole `fit` function. For more
details about the attributes of the :class:`~CallbackContext` object, see :class:`its
documentation <CallbackContext>`.

.. dropdown:: The estimator argument

    The estimator instance received by the hooks as a mandatory argument is in the same
    state as it was when calling the hook during `fit`. Therefore it is not expected to
    be fully fitted (except for the :meth:`~Callback.teardown` hook). Thus this
    estimator instance is probably not able to call its `predict` or `transform`
    methods.


`setup` and `teardown`
^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~Callback.setup` and :meth:`~Callback.teardown` hooks are only called once,
respectively at the start and end of `fit`. They take care of setting up and tearing
down the callback.

`on_fit_task_begin` and `on_fit_task_end`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~FitCallback.on_fit_task_begin` and :meth:`~FitCallback.on_fit_task_end`
hooks are called respectively at the beginning and end of each task in `fit`.

These hooks can accept optional keyword only arguments in their signatures. The possible
key and values for these optional keyword arguments are:

- `X` : array-like
    The training data of the estimator.

- `y` : array-like
    The training targets of the estimator.

- `metadata` : dict
    A dictionary containing the training metadata.

- `fitted_estimator` : estimator instance
    An estimator instance ready to predict, transform, etc ... as if the fit had stopped
    at the end of the current task.

These optional arguments might or might not be provided by the estimator, depending on
its ability to generate the requested arguments in the current task. Thus the hooks
should not expect to always receive a value for each one of these arguments.

It is not mandatory to use these arguments in the hook signatures, and only the
arguments actually used by the hook should be added. The presence of the argument in the
signature signals whether the hook requires that argument or not. To minimise the
performance overhead of callbacks, only the arguments required by at least one
registered callback are provided by the estimator.

.. warning::

    These keyword argument must be defined as **keyword only**. If the kwargs are not
    keyword only, the values will not be provided to the hooks. For example, the
    following hook implementation::

        def on_fit_task_begin(self, estimator, context, X=None, y=None):
            ...

    will result in the hook never receiving any value for ``X`` or ``y``. To receive
    these arguments (when made available by the estimator), the implementation must be::

        def on_fit_task_begin(self, estimator, context, *, X=None, y=None):
            ...

The :meth:`~FitCallback.on_fit_task_end` hook returns a boolean, which, when set to
`True`, can request the estimator to interrupt its `fit` function, for example to
implement early stopping. This interruption request might or might not be taken into
account by the estimator in each task, depending on how the estimator's `fit` function
is implemented.

Auto-propagated callbacks
-------------------------

In addition to :class:`~FitCallback`, a second protocol,
:class:`~AutoPropagatedCallback`, inherits from the first one but adding a
:meth:`~AutoPropagatedCallback.max_propagation_depth` property. This protocol, and thus
this property, identifies callbacks that, when attached to a :term:`meta-estimator`,
should be propagated to the sub-estimators.

In a situation where a meta-estimator fits sub-estimators (which can potentially also
fit sub-sub-estimators, etc, ...), a non-propagated callback can be attached to an
estimator at any depth of this estimator nesting, and it will have its hooks called only
in this estimator's `fit` tasks (and in that estimator's clones `fit` task if the
estimator is cloned in a loop, such as in a grid search).

In contrary, an auto-propagated callback requires to be registered to the top level
meta-estimator of the estimator nesting, and its :meth:`~FitCallback.on_fit_task_begin`
and :meth:`~FitCallback.on_fit_task_end` hooks will be called in all `fit` tasks, up
to an estimator depth of :meth:`~AutoPropagatedCallback.max_propagation_depth`. An
auto-propagated callback can thus aggregate information across multiple levels of
meta-estimators and sub-estimators.

The :meth:`~Callback.setup` and :meth:`~Callback.teardown` hooks of an auto-propagated
callback will still be only called once, at the beginning and end of the top level
meta-estimator's `fit` function.

Minimal example
---------------

Here is a minimal example of a custom callback::


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

            return False

.. TODO: add a link to a more detailed gallery example of a custom callback.
