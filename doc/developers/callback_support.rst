.. _callbacks_support:

===========================================
Implementing callback support in estimators
===========================================

.. currentmodule:: sklearn.callback

Adding callback support in an estimator boils down to enabling the registration of
callbacks, expressing :term:`fit` as a tree of tasks, and invoking the callbacks at the
beginning and end of each of these tasks. To achieve this, scikit-learn provides the
following helpers from the :mod:`~sklearn.callback` module:

- :class:`~CallbackSupportMixin`, which enables callback registration and initializes
  callback handling at the beginning of fit.

- :class:`~CallbackContext`, which represents tasks and is the central object to manage
  callbacks during fit.

- :func:`~with_callbacks`, to guarantee proper callback teardown at the end of fit.

The CallbackSupportMixin class
------------------------------

To support callbacks, an estimator must inherit from the :class:`~CallbackSupportMixin`
class, which exposes the following methods:

- :meth:`~CallbackSupportMixin.set_callbacks`, a public method to be called by the user
  to register callbacks on the estimator.

- :meth:`~CallbackSupportMixin._init_callback_context`, which should be called at the
  beginning of fit to create the root :class:`~CallbackContext`, corresponding to the
  task that represents the entire execution of `fit`. This method also sets up the
  callbacks that are registered on the estimator.

  .. note::

      While the leading underscore signals that
      :meth:`~CallbackSupportMixin._init_callback_context` is intended for internal
      use and should not appear in auto-completion suggestions for end users, it is
      made available to developers building third-party estimators and should be
      considered part of the public API contract.

The CallbackContext class
-------------------------

The :class:`~CallbackContext` objects are responsible for invoking the callbacks at the
right time during fit. They track the different tasks of the estimator, with one context
representing each task, and capture the tree structure of the tasks involved in the
execution of the fit method.

.. _callback_task_definition:

A task is an arbitrary unit of work defined by the estimator. Usually, a task
corresponds to an iteration of the estimator's learning algorithm. They can also
correspond to steps of a pipeline, cross-validation folds, etc. As tasks can be
decomposed into subtasks, the tasks (and therefore callback contexts) have a natural
tree structure, with the root task being the whole fit task.

The callback context objects follow this tree structure, holding references to their
parent and children contexts, dynamically built during `fit`. The root context must be
created by the :meth:`~CallbackSupportMixin._init_callback_context` method.

.. dropdown:: examples of task / context trees

    As an example, KMeans has two nested loops: the outer loop is controlled by the
    `n_init` parameter, and the inner loop is controlled by the `max_iter` parameter.
    Therefore its task tree looks like this::

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

    where each innermost `iter j` task corresponds to the computation of the labels
    and centers for the full dataset. A callback registered on a KMeans estimator thus
    will be invoked at the beginning and end of the `fit` task, each of the `init i`
    tasks and each of the `iter j` tasks.

    By convention, for performance reasons and consistency across estimators, the
    innermost tasks of scikit-learn estimators, i.e. the leaves of the task tree,
    correspond to operations on the full input data (or batches for incremental
    estimators).

    When the estimator is a meta-estimator, a task leaf usually corresponds to fitting a
    sub-estimator. Therefore, this leaf and the root task of the sub-estimator actually
    represent the same task. In this case the leaf task of the meta-estimator and the
    root task of the sub-estimator are merged into a single task. The task trees of the
    meta-estimator and the sub-estimator are combined into a single task tree. For
    instance, a `Pipeline` would have a task tree that looks like this::

        Pipeline fit (root)
        ├── step 0 | StandardScaler fit
        │   └── <insert StandardScaler task tree here>
        └── step 1 | LogisticRegression fit
            └── <insert LogisticRegression task tree here>

To dynamically build the context tree and manage the callbacks during fit, the
:class:`~CallbackContext` class exposes the following methods:

- :meth:`~CallbackContext.subcontext`

  This method allows to create a context for a subtask. Callback contexts should not be
  created directly but through this method (or `_init_callback_context` for the root
  context).

- :meth:`~CallbackContext.call_on_fit_task_begin` and
  :meth:`~CallbackContext.call_on_fit_task_end`

  .. code-block:: python

      def call_on_fit_task_begin(
          self, *, estimator, X=None, y=None, metadata=None, reconstruction_attributes=None
      ) -> None: ...

      def call_on_fit_task_end(
          self, *, estimator, X=None, y=None, metadata=None, reconstruction_attributes=None
      ) -> bool: ...

  These two methods must be called respectively at the beginning and end of the task
  that the context is responsible for. As their name suggests, they call the
  `on_fit_task_begin` and `on_fit_task_end` methods of the callbacks registered on the
  estimator.

  In addition to the callback context that is implicitly passed to the callbacks, the
  keyword arguments are used to pass additional information about the state of the
  fitting process at the task. It is not expected to provide a value for all of them at
  the beginning and end of every task. Estimators are expected to provide all the values
  that they are capable to produce. Callbacks then adapt their behavior based on the
  provided values for a given task.

  .. dropdown:: The `reconstruction_attributes` kwarg

      When `call_on_fit_task_begin/end` is called, the state of the estimator at this
      task is likely to be incomplete and thus unable to predict, transform, etc ... The
      `reconstruction_attributes` kwarg expects a dictionary containing the missing
      attributes to set on the estimator to make it ready as if the fit had stopped at
      this task.

      The callback context will copy the state of the estimator at this task, set the
      reconstruction attributes and pass the resulting estimator to the callbacks as
      `fitted_estimator`.

      If no additional attributes are needed to make the estimator ready, an empty
      dictionary should be passed instead of leaving the default value otherwise the
      callback context won't pass a `fitted_estimator` to the callbacks.

  .. dropdown:: Lazy-loading of the kwargs

      For each of these kwargs, a callable can be provided instead of the actual value.
      When it is the case, if a callback requires the kwarg, the callback context will
      evaluate the callable and forward the returned value to the callback. This
      mechanism enables lazy-loading the kwarg values, to avoid potentially costly
      computations when no callback requires a kwarg value.

      To prevent performance degradations, estimators should lazily pass quantities
      that are expensive to compute.

  .. dropdown:: Interrupting `fit`

      The `call_on_fit_task_end` method returns a boolean, which can be used to
      interrupt the current level of iterations, to implement early stopping for
      instance. It returns `True` if any callback signaled to stop the `fit` process at
      the end of this task and `False` otherwise.

- :meth:`~CallbackContext.propagate_callback_context`.

  This method allows to combine the context trees of individual estimators and
  meta-estimators in estimator compositions (e.g. a `GridSearchCV` on a
  `LogisticRegression`) into a single context tree, rooted at the fit of the top level
  estimator.

  It should be used in a meta-estimator, on a context corresponding to the task of
  fitting a sub-estimator. This task is both a leaf task of the meta-estimator and the
  root task of the sub-estimator. Their corresponding contexts are thus merged into a
  single context in the combined tree.

  In addition, `propagate_callback_context` is a context manager that propagates the
  auto-propagated callbacks from the meta-estimator to the sub-estimator such that they
  are called at the tasks of the sub-estimator as well. It also clears the propagated
  callbacks on exit such that the fitted sub-estimator no longer holds locally
  registered callbacks.

The with_callbacks decorator
------------------------------

For third-party estimators implementing callback support, the `fit` method should be
decorated with the :func:`~with_callbacks` decorator. This decorator guarantees that the
callbacks are torn down after `fit` finishes, even if it exits on an error.

For scikit-learn's built-in estimators, the :func:`~sklearn.base._fit_context` decorator
already takes care of the callbacks teardown, thus `with_callbacks` should not be used.

Minimal example
---------------

Here is a typical implementation of callback support in a custom estimator:

.. code-block:: python

    from sklearn.callback import CallbackSupportMixin, with_callbacks


    class MyEstimator(CallbackSupportMixin):
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

                if subcontext.call_on_fit_task_end(estimator=self, X=X, y=y):
                    break

            callback_ctx.call_on_fit_task_end(estimator=self, X=X, y=y)

            return self


.. TODO: add link to example
.. For a more detailed example of how to make a custom estimator or meta-estimator
.. compatible with scikit-learn's callback API, you can refer to this example :
