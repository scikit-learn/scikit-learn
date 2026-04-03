.. _custom_callbacks:

===========================
Developing custom Callbacks
===========================

.. currentmodule:: sklearn.callback

Scikit-learn offers a callback API to use built-in or custom callbacks with compatible
estimators. This section is intended for developers who wish to implement custom
callbacks.

A callback is an object that can be registered to an estimator. The callback holds
methods which will get called at different specific steps during the fitting process of
the estimator it is attached to. These methods are called hooks and the steps in ``fit``
that trigger the hooks are referred to as ``fit`` tasks.

The hooks can receive contextual information about the task in which they get called
through objects of the :class:`~CallbackContext` class and optional keyword arguments in
the hooks signature.

In scikit-learn, callbacks are classes following the :class:`~FitCallback` protocol.
This protocol defines four hook methods that each callback class must implement and
which are invoked at different steps of the estimator's fit process. To be a compatible
callback, a custom class must then simply implement these four hooks.

Callback hooks
--------------

The four hooks a callback needs to implement are:

.. method:: FitCallback.setup(estimator, context) -> None
    :noindex:

    Called only once at the start of the ``fit`` method, this hook is responsible for
    the callback's setup.

    The arguments are the estimator instance calling the callback, and the
    :class:`~CallbackContext` object holding the contextual information of the task
    representing the whole ``fit`` function.

.. method:: FitCallback.teardown(estimator, context) -> None
    :noindex:

    Called only once at the end of the ``fit`` method, this hook is responsible for the
    teardown of the callback. It is called even if the fit method exits on an error.

    The arguments are the estimator instance calling the callback, and the
    :class:`~CallbackContext` object holding the contextual information of the task
    representing the whole ``fit`` function.

.. method:: FitCallback.on_fit_task_begin(estimator, context, *, **kwargs) -> None
    :noindex:

    This hook is called at the beginning of each task in the ``fit`` function of the
    estimator the callback is registered to.

    The mandatory arguments are the estimator instance calling the callback, and the
    :class:`~CallbackContext` object holding the contextual information of the task.
    Additional keyword only arguments can be defined in the method signature, and they
    might or might not be provided by the estimator, depending on its ability to
    generate the requested arguments in the current task. The list of the possible
    keyword arguments that can be requested by this hook is described in the
    :ref:`related section <callback_hooks_kwargs>` below.

.. method:: FitCallback.on_fit_task_end(estimator, context, *,  **kwargs) -> bool
    :noindex:

    This hook is called at the end of each task in the ``fit`` function of the
    estimator the callback is registered to.

    The mandatory arguments are the estimator instance calling the callback, and the
    :class:`~CallbackContext` object holding the contextual information of the task.
    Additional keyword only arguments can be defined in the method signature, and they
    might or might not be provided by the estimator, depending on its ability to
    generate the requested arguments in the current task. The list of the possible
    keyword arguments that can be requested by this hook is described in the
    :ref:`related section <callback_hooks_kwargs>` below.

    This hook returns a boolean, which, when set to ``True``, can request the estimator
    to interrupt its ``fit`` function, for example to implement early stopping. This
    interruption request might or might not be taken into account by the estimator in
    each task, depending on how the estimator's ``fit`` function is implemented.

.. warning::

    The optional keyword argument that can be used in the signature of
    :meth:`~FitCallback.on_fit_task_begin` and :meth:`~FitCallback.on_fit_task_end` must
    be defined as **keyword only**. If the kawrgs are not keyword only, the values will
    not be provided to the hooks. For example, the following hook implementation::

        def on_fit_task_begin(self, estimator, context, X=None, y=None):
            ...

    will result in the hook never receiving any value for ``X`` or ``y``. To receive
    these arguments (when made available by the estimator), the implementation must be::

        def on_fit_task_begin(self, estimator, context, *, X=None, y=None):
            ...


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
in this estimator's ``fit`` tasks (and in that estimator's clones ``fit`` task if the
estimator is cloned in a loop, such as in a grid search).

In contrary, an auto-propagated callback requires to be registered to the root
meta-estimator of the estimator nesting, and its :meth:`~FitCallback.on_fit_task_begin`
and :meth:`~FitCallback.on_fit_task_end` hooks will be called in all ``fit`` tasks, up
to an estimator depth of :meth:`~AutoPropagatedCallback.max_propagation_depth`. An
auto-propagated callback can thus aggregate information across multiple levels of
meta-estimators and sub-estimators.

The :meth:`~Callback.setup` and :meth:`~Callback.teardown` hooks of an auto-propagated
callback will still be only called once, at the beginning and end of the root
meta-estimator's ``fit`` function.

The CallbackContext attributes
------------------------------

The callback hooks receive as a mandatory argument a :class:`~CallbackContext` object
which holds contextual information about the task in which the callback hook is being
called. This contextual information takes the form of public attributes of the
:class:`~CallbackContext` object. These public attributes are:

- ``task_name`` : str
        The name of the task this context is responsible for.

- ``task_id`` : int
    The identifier of the task this context is responsible for.

- ``max_subtasks`` : int or None
    The maximum number of children tasks for this task. 0 means it does not have any
    child, None means the maximum number of subtasks is not known in advance.

- ``estimator_name`` : str
    The name of the estimator that holds this context.

- ``parent`` : CallbackContext or None
    The parent context of this context. None if this context is the root.

- ``root_uuid`` : uuid.UUID instance
    The UUID of the root context. All contexts in the same task tree have the same
    root UUID that is used to identify the task tree itself.

- ``source_estimator_name`` : str or None
    If the estimator of the current task is a sub-estimator inside a meta-estimator, and
    the current task is the root task of that sub-estimator's ``fit``, this attribute
    contains the name of the meta-estimator. Otherwise it is None.

- ``source_task_name`` : str or None
    If the estimator of the current task is a sub-estimator inside a meta-estimator, and
    the current task is the root task of that sub-estimator's ``fit``, this attribute
    contains the name of the task in the meta-estimator's ``fit`` that called the
    sub-estimator's ``fit``. The current task and the source task correspond to the same
    task : the ``fit`` of the sub-estimator.

.. _callback_hooks_kwargs:

Optional keyword arguments passed to hooks
------------------------------------------

The :meth:`~FitCallback.on_fit_task_begin` and :meth:`~FitCallback.on_fit_task_end`
hooks can request extra contextual information about the ``fit`` task by having optional
keyword only arguments in their signature. The list of the possible keys and
corresponding values for these ``kwargs`` is as follows.

- `X` : array-like
    The training data of the estimator.

- `y` : array-like
    The training targets of the estimator.

- `metadata` : dict
    A dictionary containing the training and validation metadata.

- `fitted_estimator` : estimator instance
    An estimator instance ready to predict, transform, etc ... as if the fit had stopped
    at the end of the current task.

Minimal example
---------------

Here is a minimal example of a custom non-propagated callback::

    import numpy as np


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
            if (
                fitted_estimator is not None
                and X is not None
                and y is not None
                and hasattr(fitted_estimator, "predict")
            ):
                mean_abs_error = np.abs(y - fitted_estimator.predict(X)).mean()
                msg += f" With a mean absolute error of {mean_abs_error}."
            print(msg)

            return False

To make this callback an auto-propagated one, simply add a ``max_propagation_depth``
property to the class::

    class MyPropagatedCallback:

        def __init__(self, max_propagation_depth):
            self.max_propagation_depth = max_propagation_depth

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
            if (
                fitted_estimator is not None
                and X is not None
                and y is not None
                and hasattr(fitted_estimator, "predict")
            ):
                mean_abs_error = np.abs(y - fitted_estimator.predict(X)).mean()
                msg += f" With a mean absolute error of {mean_abs_error}."
            print(msg)

            return False

.. TODO: add a link to a more detailed gallery example of a custom callback.
