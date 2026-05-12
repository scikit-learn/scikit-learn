.. currentmodule:: sklearn.callback

.. _callbacks_user:

Callbacks
=========

.. note::

  The callback API is experimental, and is not yet implemented for all estimators.
  Please refer to the :ref:`list of callback-compatible estimators
  <callback_compatible_estimators>` for more information. It may change without the
  usual deprecation cycle.

This guide demonstrates how to use scikit-learn's callbacks on compatible estimators.
For information about how to implement the callback API, you can refer to the following
sections of the developer's guide:

- the :ref:`Implementing callback support in estimators <callbacks_support>` section for
  making estimators compatible with callbacks.

- the :ref:`Developing callbacks <developing_callbacks>` section for how to implement a
  new callback.

In scikit-learn, callbacks are objects from the :mod:`~sklearn.callback` module that can
be registered on an estimator to insert custom logic like monitoring progress or
metrics, without modifying the underlying learning algorithm. The registered callbacks
are invoked at specific steps of the fitting process.

Registering callbacks
*********************

Estimators that support callbacks expose a :meth:`~CallbackSupportMixin.set_callbacks`
method to register callbacks on them. The following example shows how to register a
:class:`~ProgressBar` callback on a :class:`~sklearn.linear_model.LogisticRegression`::

    >>> from sklearn.callback import ProgressBar
    >>> from sklearn.linear_model import LogisticRegression
    >>> progress_bar = ProgressBar()
    >>> logreg = LogisticRegression(max_iter=200)
    >>> logreg.set_callbacks(progress_bar)
    LogisticRegression(max_iter=200)

Now that the progress bar is registered on the estimator, calling its `fit` method will
display a progress bar::

    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> logreg.fit(X, y)
    LogisticRegression - fit ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
    LogisticRegression(max_iter=200)

Multiple callbacks can be registered on the same estimator, for example a
:class:`~ScoringMonitor` callback can be registered in addition to the
:class:`~ProgressBar`::

    >>> from sklearn.callback import ScoringMonitor
    >>> scoring_monitor = ScoringMonitor(scoring="accuracy")
    >>> logreg.set_callbacks(progress_bar, scoring_monitor)
    LogisticRegression(max_iter=200)

Callback invocation
*******************

During `fit`, the callbacks are invoked at the start and end of each task, where tasks
are arbitrary units of work defined by the estimator. Usually, tasks correspond to
iterations of the estimator's learning algorithm, but they can also correspond to more
abstract operations like fitting an estimator, steps of a pipeline, cross-validation
folds, etc. Within `fit`, tasks are divided into subtasks, which can themselves be
divided and so on, giving them a natural :ref:`tree structure <example_task_tree>` where
fitting the estimator is the root task.

This tree structure will usually be reflected in a callback's generated objects and
tasks will be identified by their name, id, and a reference to their parent task.
For tasks that have a natural ordering, like the iterations of a learning algorithm, the
ids are consecutive integers starting from 0. Some callbacks may provide additional
contextual information about the tasks. Here's an example of the logs of the
:class:`~ScoringMonitor`::

    >>> logreg.fit(X, y)
    LogisticRegression - fit ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
    LogisticRegression(max_iter=200)
    >>> scoring_monitor.get_logs().data_as_pandas[["task_name", "task_id", "accuracy"]]
          task_name  task_id  accuracy
    0           fit        0  0.973...
    1    lbfgs-iter        0  0.333...
    2    lbfgs-iter        1  0.333...
    3    lbfgs-iter        2  0.666...
    4    lbfgs-iter        3  0.926...
    ...

.. The number of rows can vary depending on numerical specificities of the platform, so
   we only report the first lines here.

Usage with meta-estimators
**************************

When using callbacks in estimator compositions, involving estimators and
:term:`meta-estimators`, we distinguish two types of callbacks: regular and
auto-propagated. They serve different purposes and are meant to be registered on
different estimators or meta-estimators in the composition.

Regular callbacks
-----------------

Regular callbacks are meant to be invoked within the `fit` of a given estimator. Their
goal is usually to track the learning process of that estimator.
:class:`~ScoringMonitor`, for example, records the scores at each iteration of a model.
A regular callback can be registered on an estimator at any level of a composition. If a
regular callback is registered on an estimator that is :term:`cloned` by a
meta-estimator, possibly multiple times, that callback will be invoked in each one of
the `fit` executions of the clones.

For example, when tuning the hyperparameters of a
:class:`~sklearn.linear_model.LogisticRegression` using a
:class:`~sklearn.model_selection.GridSearchCV`, a :class:`~ScoringMonitor` can be
registered on the :class:`~sklearn.linear_model.LogisticRegression` to monitor the
scores of the logistic regression model for each parameter combination and each fold of
the grid search::

    >>> from sklearn.model_selection import GridSearchCV
    >>> scoring_monitor = ScoringMonitor(scoring="accuracy")
    >>> logreg = LogisticRegression(max_iter=200).set_callbacks(scoring_monitor)
    >>> grid_search = GridSearchCV(logreg, {"C": [10, 1, 0.1]})
    >>> grid_search.fit(X, y)
    GridSearchCV(estimator=LogisticRegression(max_iter=200),
                 param_grid={'C': [10, 1, 0.1]})
    >>> log = scoring_monitor.get_logs().data_as_pandas
    >>> # show the scores at the end of each fit of the search
    >>> log[log["parent_task_id_path"] == (0,0)][["task_name", "task_id", "accuracy"]]
       task_name  task_id  accuracy
    1        fit        0  0.975...
    2        fit        1  0.966...
    3        fit        2  0.958...
    4        fit        3  0.975...
    5        fit        4  0.966...
    6        fit        5  0.958...
    7        fit        6  0.991...
    8        fit        7  0.983...
    9        fit        8  0.958...
    10       fit        9  0.991...
    11       fit       10  0.983...
    12       fit       11  0.958...
    13       fit       12  0.975...
    14       fit       13  0.975...
    15       fit       14  0.941...

.. TODO(callbacks): link to an example of how to use and plot the scores from the logs

Auto-propagated callbacks
-------------------------

Auto-propagated callbacks are meant to be invoked within the `fit` of all
(meta-)estimators in an estimator composition. Their goal is usually to report more
general information about the status at each step of the composition.
:class:`~ProgressBar` for instance displays nested progress bars for the
meta-estimators, their sub-estimators and so on. When registered on a meta-estimator,
an auto-propagated callback will automatically be registered on all its sub-estimators
that support callbacks.

.. dropdown:: registration restrictions

    Auto-propagated callbacks are designed to be registered on the top-level
    meta-estimator of an estimator composition. If some of its sub-estimators already
    have auto-propagated callbacks registered on them, an error will be raised.

    If the top-level meta-estimator doesn't itself support callbacks, then its
    sub-estimators are allowed to have auto-propagated callbacks registered on them.
    However, be aware that this is not optimal and that callbacks might not be able to
    deliver their full capabilities.

Let's add progress bars to the example of the previous section, tuning the
hyperparameters of a :class:`~sklearn.linear_model.LogisticRegression`. The
:class:`~ProgressBar` callback needs to be registered on the grid search::

    >>> grid_search.set_callbacks(ProgressBar())
    GridSearchCV(estimator=LogisticRegression(max_iter=200),
                 param_grid={'C': [10, 1, 0.1]})
    >>> grid_search.fit(X, y) # doctest: +SKIP
    GridSearchCV - fit                                                           ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      GridSearchCV - search #0                                                   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #0  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #3  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #6  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #9  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #12 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #1  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #4  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #7  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #10 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #13 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #2  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #5  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #8  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #11 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        GridSearchCV - candidate-split-evaluation | LogisticRegression - fit #14 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
      GridSearchCV - refit-with-best-params | LogisticRegression - fit #1        ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
    GridSearchCV(estimator=LogisticRegression(max_iter=200),
                 param_grid={'C': [10, 1, 0.1]})

.. The doctest skip is here because the CI uses much shorter lines which produces an
   output that crops out the progress bars.

.. dropdown:: Control the propagation depth

    Auto-propagated callbacks can be distinguished from regular callbacks by the fact
    that they have a `max_propagation_depth` attribute. This attribute indicates the
    maximum depth of nested estimators at which the callback will be propagated. In the
    previous example, the grid search is at depth 0 and the logistic regression is at
    depth 1 for instance.

    If `max_propagation_depth` is set to 0, the callback will not be propagated to any
    sub-estimators and will only be invoked for the (meta-)estimator it is registered
    on. If it is set to `None`, it will be propagated to all nested levels of the
    estimator composition.

Scikit-learn's built-in callbacks
*********************************

The built-in callbacks currently available in scikit-learn are the following:

============================== =========================================================
Callback                       Description

============================== =========================================================
:class:`~ProgressBar`          Display progress bars.

:class:`~ScoringMonitor`       Log a scoring metric at the end of each task during
                               fit.
============================== =========================================================

.. _callback_compatible_estimators:

Callback Support Status
***********************

The development of support for callbacks in estimators is in progress. Here is a list of
the estimators that support callbacks:

- :class:`sklearn.linear_model.LogisticRegression`
- :class:`sklearn.model_selection.GridSearchCV`
- :class:`sklearn.model_selection.HalvingGridSearchCV`
- :class:`sklearn.model_selection.HalvingRandomSearchCV`
- :class:`sklearn.model_selection.RandomizedSearchCV`
- :class:`sklearn.pipeline.Pipeline`
- :class:`sklearn.preprocessing.StandardScaler`
