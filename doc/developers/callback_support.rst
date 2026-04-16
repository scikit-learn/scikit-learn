.. _callbacks_support:

===========================
Supporting the Callback API
===========================

.. currentmodule:: sklearn.callback

Scikit-learn offers a callback API to use built-in or custom callbacks with compatible
estimators. This section is intended for developers who wish to implement callback
support in estimators.

For a general introduction to callbacks and how they are used in sciki-learn, see the
user guide.

.. TODO: add link to callback user guide

Adding callback support in an estimator boils down to enabling the registration of
callbacks, expressing :term:`fit` as a tree of tasks, and invoking the callbacks at the
beginning and end of each of these tasks. To achieve this, scikit-learn provides the
following helpers from the :mod:`~sklearn.callback` module:

- :class:`~CallbackSupportMixin`, which enables callback registration and initializes
  callback handling at the beginning of fit.

- :class:`~CallbackContext`, which represents tasks and is the central object to manage
  callbacks during fit.

- :func:`~with_callbacks`, to guarantee proper callback teardown at the end of fit.

CallbackSupportMixin
--------------------

To support callbacks, an estimator must inherit from the :class:`~CallbackSupportMixin`
class, which exposes the following methods:

- :meth:`~CallbackSupportMixin.set_callbacks`, a public method to be called by the user
  to register callbacks on the estimator.

- :meth:`~CallbackSupportMixin._init_callback_context`, which should be called at the
  beginning of fit to create the root :class:`~CallbackContext`. This context
  corresponds to the task representing the entire execution of ``fit``.

  .. note::

      While the leading underscore signals that
      :meth:`~CallbackSupportMixin._init_callback_context` is intended for internal
      use and should not appear in auto-completion suggestions for end users, it is
      made available to developers building third-party estimators and should be
      considered part of the public API contract.

CallbackContext
---------------

The :class:`~CallbackContext` objects are responsible for invoking the callbacks at the
right time during fit. They track the different tasks of the estimator, with one context
representing each task, and capture the tree structure of the tasks involved in the
execution of the fit method.

A task is an arbitrary unit of work defined by the estimator. Usually, a task
corresponds to an iteration of the estimator's learning algorithm, a step of a pipeline,
a cross validation fold, etc. As tasks can be decomposed into subtasks, the tasks (and
therefore callback contexts) have a natural tree structure, with the root task being the
whole fit task.

.. dropdown:: The task tree

    During its ``fit`` process, an estimator dynamically builds a tree of tasks. The
    root task node represents the whole ``fit`` function itself, and each child node
    represents a sub-task. Such child task is tipycally one step of a loop, with nested
    loops corresponding to nested sub-tasks. The callbacks attached to an estimator will
    be invoked for each task defined in the estimator. To allow callbacks to be generic
    and reusable across estimators, the innermost tasks, i.e. the leaves of the task
    tree, must correspond to operations on the full input data (or batches for
    incremental estimators).

    As an example, KMeans has two nested loops: the outer loop is controlled by `n_init`
    and the inner loop is controlled by `max_iter`. Its task tree looks like this::

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

    where each innermost ``iter j`` task corresponds to the computation of the labels
    and centers for the full dataset. A callback attached to a KMeans estimator thus
    will be invoked in the ``fit`` node, each of the ``init i`` nodes and each of the
    ``iter j`` nodes.

    When the estimator is a meta-estimator, a task leaf usually corresponds to fitting a
    sub-estimator. Therefore, this leaf and the root task of the sub-estimator actually
    represent the same task. In this case the leaf task of the meta-estimator and the
    root task of the sub-estimator are merged into a single task.

    For instance, a `Pipeline` would have a task tree that looks like this::

        Pipeline fit (root)
        ├── step 0 | StandardScaler fit
        │   └── <insert StandardScaler task tree here>
        └── step 1 | LogisticRegression fit
            └── <insert LogisticRegression task tree here>

    In that case, a callback attached to the meta-estimator will be invoked either only
    in the task nodes of the met-estimator (here ``fit``, ``setp 1`` and ``step 2``) if
    the callback is not propagated; or in both the tasks of the meta-estimator and those
    of its sub-estimators (here the tasks of the ``StandardScaler`` and
    ``LogisticRegression`` task trees) if the callback is propagated.

.. _callback_context_methods:

The callback context's methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`~CallbackContext` objects are to be used through four methods.

- The first two methods are :meth:`~CallbackContext.call_on_fit_task_begin` and
  :meth:`~CallbackContext.call_on_fit_task_end`:

  .. code-block:: python

      CallbackContext.call_on_fit_task_begin(
          estimator, X=None, y=None, metadata=None, reconstruction_attributesNone
      ) -> None

      CallbackContext.call_on_fit_task_end(
          estimator, X=None, y=None, metadata=None, reconstruction_attributesNone
      ) -> bool

  These two methods must be called respectively at the beginning and end of the
  context's task. They invoke the corresponding method
  (:meth:`~FitCallback.on_fit_task_begin` or :meth:`~FitCallback.on_fit_task_end`) of
  each callback attached to the estimator.

  The optional keyword arguments can be used to provide additional contextual
  information to the callbacks. Although optional, certain kwargs might be necessary for
  certain callbacks to function, and not providing them will result in these callbacks
  not doing anything. Therefore, in each task, an estimator should provide all the
  values it is capable to produce to be compatible with most callbacks. These kwargs
  are:

  - `X`: array-like
      The training data of the estimator.

  - `y`: array-like
      The training targets of the estimator.

  - `metadata`: dict
      A dictionary containing the training metadata.

  - `reconstruction_attributes`: dict
      A dictionary of the fitted attributes of the estimator, used by the callback
      context to generate a `fitted_estimator`, i.e. an estimator instance ready to
      predict, transform, etc ... as if the fit had stopped at the beginning/end of this
      task. The `fitted_estimator` is the object which will be passed to the callbacks,
      if required.

  .. dropdown:: Lazy-loading of the kwargs

      For each of these kwargs, a callable can be provided instead of the actual value.
      When it is the case, if a callback requires the kwarg, the callback context will
      execute the callable and forward the returned value to the callback. This
      mechanism enables lazy-loading the kwarg values, to avoid potentially costly
      computations when no callback require a kwarg value.

  The :meth:`~CallbackContext.call_on_fit_task_end` method returns a boolean, which is
  set to `True` if any callback requests the `fit` process to stop (for example to
  perform early stopping). The estimator may thus use this boolean value to condition
  the interruption of its fitting.

- The third method is :meth:`~CallbackContext.subcontext`:

  .. code-block:: python

      CallbackContext.subcontext(
          task_name="", task_id=None, max_subtasks=0, sequential_subtasks=True
      ) -> CallbackContext instance

  This method is to be used when the context's task has a sub-task. It returns a new
  child context responsible for that sub-task. This sub-context can then be used to
  create sub-sub-contexts, and so on.

- The fourth method is :meth:`~CallbackContext.propagate_callback_context`:

  .. code-block:: python

      CallbackContext.propagate_callback_context(sub_estimator) -> self

  When composing estimators, the trees of individual estimators are combined into a
  single tree, rooted at the fit of the top level estimator. This method is used to do
  the tree combination, thus it is to be used in a meta-estimator, on a context
  corresponding to the task of fitting a sub-estimator. This task is both a leaf task of
  the meta-estimator's task tree, and the root task of the sub-estimator's task tree.
  This method will propagate the context and callbacks from the meta-estimator's leaf
  task to the sub-estimator's root task, effectively merging these two tasks into the
  same task in the combined tree.


The `with_callbacks` decorator
--------------------------------

The `fit` function of the estimator must be decorated with the :func:`~with_callbacks`
decorator. This decorator takes care of the callbacks' teardown after `fit` finishes,
even if it exits on an error. For scikit-learn's built-in estimators, the
:func:`~sklearn.base._fit_context` decorator already takes care of the callbacks
teardown, thus :func:`~with_callbacks` should not be used.

Minimal example
---------------

Here is a minimal example of a custom estimator supporting callbacks:

.. code-block:: python

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

.. TODO: add link to example
