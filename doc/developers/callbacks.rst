.. _callbacks_api:

=================================
Developing with the Callbacks API
=================================

.. currentmodule:: sklearn.callback

Scikit-learn offers a callback API to use built-in or custom callbacks with compatible
estimators. This section is intended for developers who wish to develop or maintain
callbacks or estimators that support callbacks.

Callbacks
---------

In scikit-learn, callbacks are classes following the :class:`~_base.Callback` protocol.
This protocol defines three methods, callback hooks, which are invoked at different
steps of the estimator's fit process. The three hooks of a callback are
:meth:`~sklearn.callback._base.Callback.on_fit_begin`,
:meth:`~_base.Callback.on_fit_task_end` and :meth:`~_base.Callback.on_fit_end`, which
are respectively called at the start of the ``fit`` method, at the end of each task in
``fit`` and at the end of the ``fit`` method.

A second protocol, :class:`~_base.AutoPropagatedCallback`, inherits from the first one
but adds a :meth:`~_base.AutoPropagatedCallback.max_estimator_depth` property. This
protocol identifies callbacks that can be propagated to sub-estimators when registered
to :term:`meta-estimators`.

.. TODO: add a link to an example doc implementing a custom callback.

Task tree
---------

The callback hook :meth:`~callback._base.Callback.on_fit_task_end` is called at the end
of what we call a task in ``fit``. This subsection develops this concept of task.

In scikit-learn estimators, a ``fit`` task is usually one step of a loop, with nested
loops corresponding to nested tasks. In general, a task can be whatever unit of work the
estimator's developer wants it to be. However, if the defined tasks are too unusual,
this might compromise the compatibility of the estimator with scikit-learn's builtin
callbacks.

Tasks have an inherent tree structure, where each task can be decomposed into subtasks,
and so on. The root of the tree represents the whole fit task.

Each loop in the estimator represents a parent task, and each iteration of that loop
represents a child task. To allow callbacks to be generic and reusable across
estimators, the innermost tasks, i.e. the leaves of the task tree, correspond to
operations on the full input data (or batches for incremental estimators).

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
    ├── step 0 | preprocessor fit
    │   └── <insert preprocessor task tree here>
    └── step 1 | estimator fit
        └── <insert estimator task tree here>

Concretely, the tree structure is created dynamically and abstracted in an object named
:class:`~callback._callback_context.CallbackContext`. There is one context for each
task, and the context is responsible for calling the callback hooks for its task and
creating contexts for the child tasks.

A callbcak that is not an :class:`~_base.AutoPropagatedCallback` can be attached to an
estimator at any depth of the task tree, and its hooks will be called only for this task
node (and its siblings if the estimator is cloned in a loop, such as in a grid search).
An :class:`~_base.AutoPropagatedCallback` callback needs to be registered on the
estimator at the root of the task tree. The
:meth:`~sklearn.callback._base.Callback.on_fit_begin` and
:meth:`~_base.Callback.on_fit_end` hooks of that callback will be only called in the
root task, but the :meth:`~_base.Callback.on_fit_task_end` will be called in each node
of the tree (thus in each sub-estimator).

CallbackContext
---------------

The :class:`~callback._callback_context.CallbackContext` is the object that must be
used in the implementation of an estimator to support callbacks. It is responsible for
calling the callback hooks at the right time, hence the need to hold the task tree
structure. A context is created at the beginning of fit, and then sub-contexts are
created for each child task. The callback hooks are invoked through the
:meth:`~callback._callback_context.CallbackContext.eval_on_fit_begin`,
:meth:`~callback._callback_context.CallbackContext.eval_on_fit_task_end` and
:meth:`~callback._callback_context.CallbackContext.eval_on_fit_end` methods.

Here is a minimal example of how a callback context can be handled in an estimator::


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


For a detailed example of how to make a custom estimator or meta-estimator compatible
with scikit-learn's callback API, you can refer to `this example
<https://scikit-learn.org/stable/auto_examples/callbcaks/plot_callback_support.html>`__.

The callback context is also an object that is passed to the callback hooks to give them
contextual information about the task being executed and its position in the task tree.
In particular, this extra information takes the form of keyword arguments passed to
:meth:`~callback._callback_context.CallbackContext.eval_on_fit_task_end`. These
``kwargs`` are optional. However, for an estimator to be compatible with the largest number
of callbacks, it should provide all the values it is capable or producing during fit.
The list of the possible keys and corresponding values for these ``kwargs`` is as
follows.

- `data`: dict
    A dictionary containing the training and validation data. The possible keys are
    "X_train", "y_train", "sample_weight_train", "X_val", "y_val", "sample_weight_val".

- `stopping_criterion`: float
    Usually, iterations stop when `stopping_criterion <= tol`. This is only provided at
    the innermost level of iterations.

- `tol`: float
    Tolerance for the stopping criterion. This is only provided at the innermost level
    of iterations.

- `from_reconstruction_attributes`: callable
    A function returning a ready-to-predict, transform, etc ... estimator as if the fit
    had stopped at the end of this task.

- `fit_state`: dict
    Model-specific quantities updated during fit. This is not meant to be used by
    generic callbacks but by a callback designed for a specific estimator.


And here is a list of all the built-in callbacks in scikit-learn, and which ones of the
``kwargs`` they require to function.

========================================= ====================== =======================
Callback                                  Description            Required kwargs
========================================= ====================== =======================
:class:`~_progressbar.ProgressBar`        Display progress bars. None
========================================= ====================== =======================

Finally, here is a list of the planned callbacks to be added to the built-in ones.

========================================= ====================== =======================
Planned callback                          Description            Required kwargs
========================================= ====================== =======================
MetricMonitor                             Log a metric score     data
                                          between predictions
                                          and targets.
========================================= ====================== =======================
