.. _callbacks_support:

===========================
Supporting the Callback API
===========================

.. currentmodule:: sklearn.callback

Scikit-learn offers a callback API to use built-in or custom callbacks with compatible
estimators. This section is intended for developers who wish to implement callback
support in estimators.

A callback is an object that can be registered to an estimator. The callback holds
methods, called hooks, which will get called at different specific steps during the
fitting process of the estimator it is attached to.

There are two types of callbacks : propagated and non-propagated. When using
meta-estimators, a propagated callback must be registered to the outermost
meta-estimator, and it will be propagated to the sub-estimators, meaning its hooks will
be called in both the :term:`fit` function of the meta-estimator and the ones of its
sub-estimators. By contrast a non-propagated callback can be registered to any meta- or
sub-estimator, and its hooks will only be called in that estimator's ``fit`` function.

Supporting the callback API in an estimator boils down to enabling the registration of
callbacks, defining the different steps (referred to as "tasks") of its ``fit``
method, and calling the right hooks at the right time for each defined task. Doing so
requires the use of three types of objects:

- :class:`~CallbackSupportMixin`, which enables callback registration and initializes
  callback handling at the beginning of fit.

- :class:`~CallbackContext`, which represents tasks and is the central object to manage
  callbacks during fit.

- :func:`~with_callbacks`, to guarantee proper callback tear down at the end of fit.

CallbackContext
---------------

The :class:`~CallbackContext` objects are responsible for invoking the callbacks at the
right time during fit. They track the different tasks of the estimator, with one context
representing each task, and capture the tree structure of the tasks involved in the
execution of the fit method.

The task tree
^^^^^^^^^^^^^

During its ``fit`` process, an estimator dynamically builds a tree of tasks. The root
task node represents the whole ``fit`` function itself, and each child node represents a
sub-task, i.e. an arbitrary unit of work. Such child task is tipycally one step of a
loop, with nested loops corresponding to nested sub-tasks. The callbacks attached to an
estimator will have hooks called for each task defined in the estimator. To allow
callbacks to be generic and reusable across estimators, the innermost tasks, i.e. the
leaves of the task tree, must correspond to operations on the full input data (or
batches for incremental estimators).

As an example, KMeans has two nested loops: the outer loop is controlled by `n_init` and
the inner loop is controlled by `max_iter`. Its task tree looks like this::

    KMeans fit (root)
    ├── init 0
    │   ├── iter 0
    │   ├── iter 1
    │   ├── ...
    │   └── iter n
    ├── init 1
    │   ├── iter 0
    │   ├── ...
    │   └── iter n
    └── init 2
        ├── iter 0
        ├── ...
        └── iter n

where each innermost ``iter j`` task corresponds to the computation of the labels and
centers for the full dataset. A callback attached to a KMeans estimator thus will have
its hooks called in the ``fit`` node, each of the ``init i`` nodes and each of the
``iter j`` nodes.

When the estimator is a meta-estimator, a task leaf usually corresponds to fitting a
sub-estimator. Therefore, this leaf and the root task of the sub-estimator actually
represent the same task. In this case the leaf task of the meta-estimator and the root
task of the sub-estimator are merged into a single task.

For instance, a `Pipeline` would have a task tree that looks like this::

    Pipeline fit (root)
    ├── step 0 | StandardScaler fit
    │   └── <insert StandardScaler task tree here>
    └── step 1 | LogisticRegression fit
        └── <insert LogisticRegression task tree here>

In that case, a callback attached to the meta-estimator will have its hooks called
either only in the task nodes of the met-estimator (here ``fit``, ``setp 1`` and ``step
2``) if the callback is not propagated; or in both the tasks of the meta-estimator and
those of its sub-estimators (here the tasks of the ``StandardScaler`` and
``LogisticRegression`` task trees) if the callback is propagated.

The callback context's methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`~CallbackContext` objects are to be used through four methods :

.. method:: CallbackContext.call_on_fit_task_begin(estimator, **kwargs) -> None
    :noindex:

    This method must be called at the beginning of the context's task and wil call the
    callbacks' :meth:`~FitCallback.on_fit_task_begin` hooks.

    It takes as arguments the estimator of the current task, and optional keyword
    arguments to provide extra information to be used by the callbacks. The list of
    possible keys and values for these kwargs is described in the :ref:`related section
    <callback_context_kwargs>` below.

.. method:: CallbackContext.call_on_fit_task_end(estimator, **kwargs) -> bool
    :noindex:

    This method must be called at the end of the context's task and wil call the
    callbacks' :meth:`~FitCallback.on_fit_task_end` hooks. Each of these hooks returns a
    boolean to request the ``fit`` process to stop, and this method will return ``True``
    if any of the hooks does. The estimator may thus use this boolean value to condition
    the interruption of its fitting.

    The arguments are the estimator of the current task, and the same options for
    additional keyword arguments as :meth:`~FitCallback.on_fit_task_begin`, listed in
    the :ref:`related section <callback_context_kwargs>` below.

.. method:: CallbackContext.subcontext(task_name="", task_id=None, max_subtasks=0, sequential_subtasks=True) -> CallbackContext instance
    :noindex:

    This method is used to generate a sub-context for a sub-task of the task the current
    context is responsible for.


    The generated sub-context is returned, which can then be used to create sub-sub-contexts, and so on.

.. method:: CallbackContext.propagate_callback_context(sub_estimator) -> self
    :noindex:

    This method is used in meta-estimators to propagate the callbacks to a sub-estimator.
    It will also propagate the context, merging it with the root context of the
    sub-estimator's ``fit`` function, as the current context's task and the
    sub-estimator's root task represent the same task.

    The current context is returned.

CallbackSupportMixin
--------------------

To support callbacks, an estimator must inherit from the :class:`~CallbackSupportMixin`
class, which exposes the following methods:

- :meth:`~CallbackSupportMixin.set_callbacks`, a public method to be called by the user
  to register callbacks on the estimator.

- :meth:`~CallbackSupportMixin._init_callback_context`, which should be called at the
  beginning of fit to create the root `~CallbackContext`. This context represents the task
  of running fit as a whole.


    The generated :class:`~CallbackContext` instance is returned.


The ``with_callbacks`` decorator
--------------------------------

The ``fit`` function of the estimator must be decorated with the :func:`~with_callbacks`
decorator. This decorator takes care of calling the
:meth:`~_base._BaseCallback.teardown` hook of the callbacks after ``fit`` finishes, and
does so in a ``try`` / ``finally`` block, in order to always perform the callbacks
teardown, even if ``fit`` exits on an error.

Minimal example
---------------

Here is a minimal example of a custom estimator supporting callbacks::

    from sklearn.base import BaseEstimator
    from sklearn.callback import CallbackSupportMixin, with_callbacks


    class MyEstimator(CallbackSupportMixin, BaseEstimator):
        def __init__(self, max_iter):
            self.max_iter = max_iter

        @with_callbacks
        def fit(self, X, y):
            callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
            callback_ctx.call_on_fit_task_begin(estimator=self, X=X, y=y)

            for i in range(self.max_iter):
                subcontext = callback_ctx.subcontext(task_name="iteration")
                subcontext.call_on_fit_task_begin(estimator=self, X=X, y=y)

                # Do something

                subcontext.call_on_fit_task_end(estimator=self, X=X, y=y)

            callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

            return self


For a more detailed example of how to make a custom estimator or meta-estimator
compatible with scikit-learn's callback API, you can refer to this example :
:ref:`sphx_glr_auto_examples_callbacks_plot_callback_support.py`.

.. _callback_context_kwargs:

Contextual information passed to callbacks
------------------------------------------

Additional contextual information can be forwarded to the callback hooks from the
:meth:`~CallbackContext.call_on_fit_task_begin` and
:meth:`~CallbackContext.call_on_fit_task_end` methods, taking the form optional keyword
arguments given to these methods. The content of these ``kwargs`` will depend on the
estimator's implementation, as an estimator might not be able to produce all the
possible values. The list of the possible keys and corresponding values for these
``kwargs`` is as follows.

- `X`: array-like
    The training data of the estimator.

- `y`: array-like
    The training targets of the estimator.

- `metadata`: dict
    A dictionary containing the training and validation metadata.

- `reconstruction_attributes`: dict
    A dictionary of the fitted attributes of the estimator, used by the callback context
    to generate a `fitted_estimator`, i.e. an estimator instance ready to predict,
    transform, etc ... as if the fit had stopped at the end of this task. The
    `fitted_estimator` is the object which will be passed to the callback hooks if
    required.

For each of these kwargs, a callable can be provided instead of the actual value. When
it is the case, if a callback hook requires the kwarg, the callback context will execute
the callable and forward the returned value to the hook. This mechanism enables
lazy-loading the kwarg values, to avoid potentially costly computations when no callback
require a kwarg value.

Here is a list of all the built-in callbacks in scikit-learn, and which ones of these
``kwargs`` must be given to the :class:`~CallbackContext` methods for the callback
to function.

============================== ====================== =========================== ====================
Callback                       Description            Required kwargs for         Required kwargs for
                                                      call_on_fit_task_begin      call_on_fit_task_end
============================== ====================== =========================== ====================
:class:`~ProgressBar`          Display progress bars. None                        None
============================== ====================== =========================== ====================
