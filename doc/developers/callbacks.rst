.. _callbacks_api:

================================
Developing with the Callback API
================================

.. currentmodule:: sklearn.callback

Scikit-learn offers a callback API to use built-in or custom callbacks with compatible
estimators. This section is intended for developers who wish to implement callback
support in estimators, or develop custom callbacks.

Callbacks
---------

In scikit-learn, callbacks are classes following the :class:`~Callback` protocol.
This protocol defines three methods (referred to as callback hooks) which are invoked at
different steps of the estimator's fit process. The three hooks of a callback are :

.. method:: Callback.on_fit_begin(estimator) -> None
    :noindex:

    Called only once at the start of the ``fit`` method and responsible for the
    callback's setup.

.. method:: Callback.on_fit_task_end(estimator, context, **kwargs) -> bool
    :noindex:

    Called at the end of each task in ``fit`` and implementing the main work of the
    callback. It returns a boolean which can be used by the estimator to stop its fit
    process. In addition to the estimator the callback is attached to, it takes a
    ``context`` argument, a :class:`~sklearn.callback.CallbackContext` object holding
    contextual information and described in more detailed in the
    :ref:`callback_context_section` below. Optional keyword arguments can be provided
    depending on the estimator, with the list of possible keys and values described in
    the :ref:`related section <contextual_info_section>` below.

.. method:: Callback.on_fit_end(estimator) -> None
    :noindex:

    Called at the end of the ``fit`` method and responsible for the tear down of the
    callback. It is called even if the fit method does not complete.

.. TODO: add a link to an example doc implementing a custom callback.

Task tree
---------

During its ``fit`` process, an estimator dynamically builds a tree of tasks. The root
task node represents the whole ``fit`` function itself, and each child node represents a
sub-task, i.e. an arbitrary unit of work. Such child task is tipycally one step of a
loop, with nested loops corresponding to nested sub-tasks.

For instance, KMeans has two nested loops: the outer loop is controlled by `n_init` and
the inner loop is controlled by `max_iter`. Its task tree looks like this::

    KMeans fit
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
centers for the full dataset.

When the estimator is a meta-estimator, a task leaf usually corresponds to fitting a
sub-estimator. Therefore, this leaf and the root task of the sub-estimator actually
represent the same task. In this case the leaf task of the meta-estimator and the root
task of the sub-estimator are merged into a single task.

For instance, a `Pipeline` would have a task tree that looks like this::

    Pipeline fit
    ├── step 0 | StandardScaler fit
    │   └── <insert StandardScaler task tree here>
    └── step 1 | LogisticRegression fit
        └── <insert LogisticRegression task tree here>


The callback hook :meth:`~Callback.on_fit_task_end` is called at the end of each task
node of the tree. To allow callbacks to be generic and reusable across estimators, the
innermost tasks, i.e. the leaves of the task tree, must correspond to operations on the
full input data (or batches for incremental estimators).

Concretely, the tree structure is created dynamically and abstracted in an object named
:class:`~CallbackContext`. There is one context for each task, and the context is
responsible for calling the callback hooks for its task and creating contexts for the
child tasks.

Auto-propagated callbacks
-------------------------

In addition to :class:`~Callback`, a second protocol :class:`~AutoPropagatedCallback`
inherits from the first one but adding a
:meth:`~AutoPropagatedCallback.max_estimator_depth` property. This protocol identifies
callbacks that should be propagated through the task tree. A non-propagated callback can
be attached to an estimator at any depth of the task tree, and it will have its hooks
called only in this estimator's task nodes (and its siblings if the estimator is cloned
in a loop, such as in a grid search). Whereas an auto-propagated callback requires to be
registered to the root estimator of the task tree, and its
:meth:`~Callback.on_fit_task_end` hook will be called in all task node of the tree, up
to a depth of :meth:`~AutoPropagatedCallback.max_estimator_depth`. The
:meth:`~Callback.on_fit_begin` and :meth:`~Callback.on_fit_end` hooks of an
auto-propagated callback will still be only called once, in the root task.

An auto-propagated callback can thus aggregate information across the task tree. For
example, the :class:`sklearn.callback.ProgressBar` callback is auto-propagated in order
to show nested progress bars for the sub-estimators in addition to the meta-estimator,
producing outputs like the following::

  MetaEstimator - fit                            ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
    MetaEstimator - outer #0                     ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      MetaEstimator - inner | Estimator - fit #0 ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      MetaEstimator - inner | Estimator - fit #1 ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
    MetaEstimator - outer #1                     ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      MetaEstimator - inner | Estimator - fit #0 ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      MetaEstimator - inner | Estimator - fit #1 ━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00

.. _callback_context_section:

CallbackContext
---------------

The :class:`~CallbackContext` is the object that must be used in the implementation of
an estimator to support callbacks. It is responsible for calling the callback hooks at
the right time, creating child contexts for child tasks and propagating the context and
callbacks to sub-estimators.

The root context object of a task tree is created at the beginning of ``fit`` through
the :meth:`~CallbackSupportMixin._init_callback_context` method of the root estimator.
This method comes from the :class:`~CallbackSupportMixin` mixin that estimtators must
inherit from to support callbacks. The children context objects are then instantiated in
each sub-task through their parent context's :meth:`~CallbackContext.subcontext` method.

The callback hooks are invoked by the context object through its
:meth:`~CallbackContext.eval_on_fit_begin`,
:meth:`~CallbackContext.eval_on_fit_task_end` and
:meth:`~CallbackContext.eval_on_fit_end` methods.
:meth:`~CallbackContext.eval_on_fit_begin` and
:meth:`~CallbackContext.eval_on_fit_task_end` are called during fit, respectively once
at the beginning of ``fit`` and at each end of a task, as the names implie.
:meth:`~CallbackContext.eval_on_fit_end` is automatically called after ``fit``, even if
the ``fit`` call did not complete. This is achieved by decorating the ``fit`` method
with the :func:`~CallbackContext.with_callback_context` decorator, which also takes care
of the context tear-down.

Here is a minimal example of how a callback context should be handled in an estimator::


    class MyEstimator(CallbackSupportMixin, BaseEstimator):
        def __init__(self, max_iter):
            self.max_iter = max_iter

        @with_callback_context
        def fit(self, X, y):
            callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
            callback_ctx.eval_on_fit_begin(estimator=self)

            for i in range(self.max_iter):
                subcontext = callback_ctx.subcontext(task_id=i)

                # Do something

                subcontext.eval_on_fit_task_end(
                    estimator=self,
                    data={"X_train": X, "y_train": y},
                )

            return self


For a more detailed example of how to make a custom estimator or meta-estimator
compatible with scikit-learn's callback API, you can refer to this example :
:ref:`sphx_glr_auto_examples_callbacks_plot_callback_support.py`.

.. _contextual_info_section:

Contextual information passed to callbacks
------------------------------------------

The callback context is also an object that is passed to the
:meth:`~CallbackContext.eval_on_fit_task_end` hooks to give contextual
information about the task being executed and its position in the task tree.

Some extra information is also forwarded to the
:meth:`~CallbackContext.eval_on_fit_task_end` hooks as optional
keyword arguments. The content of these ``kwargs`` will depend on the estimator's
implementation, as an estimator might not be able to produce all the possible values.
The list of the possible keys and corresponding values for these ``kwargs`` is as
follows.

- `data`: dict
    A dictionary containing the training and validation data. The possible keys are
    "X_train", "y_train", "sample_weight_train", "X_val", "y_val", "sample_weight_val".

- `from_reconstruction_attributes`: callable
    A function returning a ready-to-predict, transform, etc ... estimator as if the fit
    had stopped at the end of this task.

- `fit_state`: dict
    Model-specific quantities updated during fit. This is not meant to be used by
    generic callbacks but by a callback designed for a specific estimator.


And here is a list of all the built-in callbacks in scikit-learn, and which ones of the
``kwargs`` they require to function.

============================== ====================== ==================================
Callback                       Description            Required kwargs
============================== ====================== ==================================
:class:`~ProgressBar`          Display progress bars. None
============================== ====================== ==================================

Finally, here is a list of the planned callbacks to be added to the built-in ones.

============================== ====================== ==================================
Planned callback               Description            Required kwargs
============================== ====================== ==================================
MetricMonitor                  Log a metric score     data,
                               between predictions    from_reconstruction_attributes
                               and targets.
============================== ====================== ==================================
