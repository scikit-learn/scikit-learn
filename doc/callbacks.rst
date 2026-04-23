.. currentmodule:: sklearn.callback

.. _callbacks_user:

Callbacks
=========

.. note::

  The callback API is experimental, and is not yet implemented for all estimators.
  Please refer to the :ref:`list of callback compatible estimators
  <callback_compatible_estimators>` for more information. It may change without the
  usual deprecation cycle.

This guide demonstrates how to use scikit-learn's callbacks on compatible estimators. If
you are developing a scikit-learn compatible estimator or meta-estimator and want to
make it compatible with callbacks, you can check our related developer guide. If you
want to develop a callback to use with callback compatible estimators, you can refer
to the related section of the developer guide.

.. TODO add refs to callback developer doc

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
    >>> estimator = LogisticRegression()
    >>> estimator.set_callbacks(progress_bar) # doctest: +SKIP

.. TODO: remove the doctest skip

Now that the progress bar is registered on the estimator, calling its `fit` method will
display a bar monitoring the progress of the `fit` execution::

    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> estimator.fit(X, y) # doctest: +SKIP
    LogisticRegression - fit ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:03:12

.. TODO: remove the doctest skip

Multiple callbacks can be registered on the same estimator, for example a
:class:`~ScoringMonitor` callback can be registered in addition to the
:class:`~ProgressBar` one::

    >>> from sklearn.callback import ScoringMonitor # doctest: +SKIP
    >>> scoring_monitor = ScoringMonitor(scoring="precision") # doctest: +SKIP
    >>> estimator.set_callbacks(progress_bar, scoring_monitor) # doctest: +SKIP

.. TODO: remove the doctest skips

Usage with meta-estimators
**************************

There are two types of callbacks, regular callbacks and auto-propagated ones, reflecting
two different usages in the context of an estimator composition, when using a
:term:`meta-estimator`.

Regular callbacks
-----------------

Regular callbacks are meant to track the learning process of one estimator. A regular
callback can be registered on an estimator at any level of the composition, and it will
be invoked only within the `fit` of that estimator's level. If a regular callback is
registered on a sub-estimator which is fitted multiple times by a meta-estimator, or is
cloned by the meta-estimator and the clones get fitted, that callback will be
invoked in each one of these `fit` executions.

For example in the case of a :class:`~sklearn.model_selection.GridSearchCV` applied on a
:class:`~sklearn.linear_model.LogisticRegression`, the regular callback
:class:`~ScoringMonitor` can be registered on the
:class:`~sklearn.linear_model.LogisticRegression`. Then the scores will be logged for
each logistic regression fit of each parameter combination and each fold of the grid
search::

    >>> from sklearn.model_selection import GridSearchCV
    >>> parameters = {"C": [10, 1, 0.1], "l1_ratio": [0, 1]}
    >>> scoring = ScoringMonitor(scoring="precision") # doctest: +SKIP
    >>> sub_estimator = LogisticRegression()
    >>> sub_estimator.set_callbacks(scoring) # doctest: +SKIP
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters)
    >>> meta_estimator.fit(X, y)
    GridSearchCV(estimator=LogisticRegression(),
                 param_grid={'C': [10, 1, 0.1], 'l1_ratio': [0, 1]})
    >>> scoring.get_logs() # doctest: +SKIP

.. TODO: add printed output, or maybe use the logs in a plot
.. TODO: remove the doctest skips

Auto-propagated callbacks
-------------------------

Auto-propagated callbacks are :class:`~AutoPropagatedCallback` objects, and they are
meant to aggregate information from the different levels of an estimator composition.
When registered on a meta-estimator, an auto-propagated callback can be propagated from
the meta-estimator to its sub-estimators, meaning it will get invoked in the `fit` of
both the meta-estimator and the sub-estimators. If the sub-estimators are also
meta-estimators, the auto-propagated callback can be propagated to the
sub-sub-estimators too, and so on. This propagation is controlled through the
`max_propagation_depth` argument of the callback. This argument indicates the estimator
depth up to which the callback will be propagated in an estimator composition.

.. dropdown:: Auto-propagated callbacks registration restriction

    Auto-propagated callbacks cannot be registered on an estimator that is a
    sub-estimator of a callback compatible meta-estimator. In other words, in an
    estimator composition an auto-propagated callback must be registered to the
    outermost meta-estimator that supports callbacks. Otherwise an error will be raised
    when that callback compatible meta-estimator will have its `fit` executed.

For example in the case of a :class:`~sklearn.model_selection.GridSearchCV` applied on a
:class:`~sklearn.linear_model.LogisticRegression`, the auto-propagated callback
:class:`~ProgressBar` must be registered on the
:class:`~sklearn.model_selection.GridSearchCV`, and having its `max_propagation_depth`
set to 1 means that the progressbar will show progress or both the grid search and each
logistic regression fit::

    >>> from sklearn.model_selection import GridSearchCV
    >>> parameters = {"l1_ratio": [0, 1], "fit_intercept": [True, False]}
    >>> sub_estimator = LogisticRegression()
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters)
    >>> meta_estimator.set_callbacks(ProgressBar(max_propagation_depth=1)) # doctest: +SKIP
    >>> meta_estimator.fit(X, y) # doctest: +SKIP

.. TODO: add printed output
.. TODO: remove the doctest skips

Changing the `max_propagation_depth` argument to 0 will make the auto-propagated
callback only be registered on the top-level estimator, in that case it means that only
one progress bar for the :class:`~sklearn.model_selection.GridSearchCV` will be
displayed::

    >>> meta_estimator.set_callbacks(ProgressBar(max_propagation_depth=0)) # doctest: +SKIP
    >>> meta_estimator.fit(X, y) # doctest: +SKIP

.. TODO: add printed output, or maybe use the logs in a plot
.. TODO: remove the doctest skip

Fit tasks
*********

During a callback compatible estimator's `fit`, the callbacks are invoked at the start
and end of each task. The tasks are arbitrary units of work defined by the estimator.
Usually, a task corresponds to an iteration of the estimator's learning algorithm. They
can also correspond to steps of a pipeline, cross-validation folds, etc. As tasks can be
decomposed into subtasks, they have a natural tree structure which can be reflected in
a callback's generated objects, such as logs. For more details about the tasks and their
tree structure, please see the corresponding section in the developer documentation.

.. TODO: add link to the dev doc

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

The development of support for callbacks in estimators is in progress, here is a list of
the estimators which support callbacks.

- :class:`~sklearn.linear_model.LogisticRegression`
- :class:`~sklearn.model_selection.GridSearchCV`
