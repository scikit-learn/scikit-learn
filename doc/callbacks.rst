.. currentmodule:: sklearn.callback

.. _callbacks:

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

In scikit-learn, callbacks are objects from the :mod:`callbacks` module that can be
registered on an estimator to insert custom logic like monitoring progress or metrics,
without modifying the underlying learning algorithm. The registered callbacks are
called at specific steps of the fitting process.

Registering callbacks
*********************

Estimators that support callbacks expose a :meth:`~CallbackSupportMixin.set_callbacks`
method to register callbacks on them. The following example shows how to register a :class:`~ProgressBar` callback on a :class:`~sklearn.linear_model.LogisticRegression`::

    >>> from sklearn.callback import ProgressBar
    >>> from sklearn.linear_model import LogisticRegression
    >>> progress_bar = ProgressBar()
    >>> estimator = LogisticRegression()
    >>> estimator.set_callbacks(progress_bar)

Now that the progress bar is registered on the estimator, calling its `fit` method will
display a bar monitoring the progress of the `fit` execution::

    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> estimator.fit(X, y)
    LogisticRegression - fit ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:03:12

Multiple callbacks can be registered to the same estimator by providing a list to
:meth:`~CallbackSupportMixin.set_callbacks`, for example a :class:`~ScoringMonitor`
callback can be registered in addition to the :class:`~ProgressBar` one::

    >>> from sklearn.callback import ScoringMonitor
    >>> scoring_monitor = ScoringMonitor(scoring="precision")
    >>> estimator.set_callbacks(progress_bar, scoring_monitor)

Usage with meta-estimators
**************************

:term:`Meta-estimators` are estimator that, in their `fit`, will fit other estimators
(designated as sub-estimators). For example a
:class:`~sklearn.model_selection.GridSearchCV` can be used to find the best parameters
for a :class:`~sklearn.linear_model.LogisticRegression` by fitting multiple instances of
that sub-estimator. Callbacks will behave differently in such estimator composition,
depending whether the callback is an :class:`~AutoPropagatedCallback` or not.

Auto-propagated callbacks
-------------------------

An :class:`~AutoPropagatedCallback`, such as :class:`~ProgressBar`, is a callback that
can be propagated to sub-estimators when registered to a meta-estimator. Thus its
methods can be called in the `fit` executions of both the meta-estimator and the
sub-estimators, allowing to aggregate information about the `fit` process of both types
of estimators. If the sub-estimators are also meta-estimator, the auto-propagated
callback can be propagated to the sub-sub-estimators too, and so on. This callback
propagation is controlled through the `max_propagation_depth` argument of the callback.
This argument indicates the estimator depth up to which the callback will be propagated
in an estimator composition.

Auto-propagated callbacks must be registered to the top level meta-estimator of an
estimator composition. Trying to register it to any sub-estimator will raise an error
when fitting the top level meta-estimator.

AutoPropagatedCallback examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of a
:class:`~sklearn.model_selection.GridSearchCV` applied on a
:class:`~sklearn.linear_model.LogisticRegression`, the auto-propagated callback
:class:`~ProgressBar` must be registered to the :class:`~sklearn.model_selection.GridSearchCV`::

    >>> from sklearn.model_selection import GridSearchCV
    >>> parameters = {"l1_ratio": [0, 1], "fit_intercept": [True, False]}
    >>> sub_estimator = LogisticRegression()
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters)
    >>> meta_estimator.set_callbacks(ProgressBar(max_propagation_depth=1))
    >>> meta_estimator.fit(X, y)
    #TODO add printed output

Changing the `max_propagation_depth` argument to 0 will make the auto-propagated
callback only be registered to the top-level estimator, in that case it means that only
one progressbar for the :class:`~sklearn.model_selection.GridSearchCV` will be
displayed::

    >>> meta_estimator.set_callbacks(ProgressBar(max_propagation_depth=0))
    >>> meta_estimator.fit(X, y)
    #TODO add printed output

Trying to register the auto-propagated callback on the sub-estimator will raise an
error::

    >>> sub_estimator = LogisticRegression().set_callbacks(ProgressBar())
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters)
    >>> meta_estimator.fit(X, y)
    #TODO add printed output


Non auto-propagated callbacks
-----------------------------

A non auto-propagated callback, such as :class:`~ScoringMonitor`, can be registered to
any meta or sub-estimator in an estimator composition, and will be only invoked in the
`fit` of the estimators it was registered to. If a sub-estimator with a non
auto-propagated callback is fitted multiple times by a meta-estimator, or if it is
cloned by a meta-estimator and the clones get fitted, the callback will be invoked in
each one of these `fit` executions.

Non auto-propagated examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of a :class:`~sklearn.model_selection.GridSearchCV` applied on a
:class:`~sklearn.linear_model.LogisticRegression`, the non auto-propagated callback
:class:`~ScoringMonitor` can be registered to the
:class:`~sklearn.linear_model.LogisticRegression`::

    >>> scoring = ScoringMonitor(scoring="precision")
    >>> sub_estimator = LogisticRegression().set_callbacks(scoring)
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters)
    >>> meta_estimator.fit()
    >>> scoring.get_logs()
    #TODO add printed output, or maybe use the logs in a plot

The non auto-propagated callback can also be registered to the
:class:`~sklearn.model_selection.GridSearchCV`::

    >>> sub_estimator = LogisticRegression()
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters).set_callbacks(scoring)
    >>> meta_estimator.fit()
    >>> scoring.get_logs()
    #TODO add printed output, or maybe use the logs in a plot

Or it can be registered to both the :class:`~sklearn.model_selection.GridSearchCV` and
the :class:`~sklearn.linear_model.LogisticRegression`::

    >>> sub_estimator = LogisticRegression().set_callbacks(scoring)
    >>> meta_estimator = GridSearchCV(sub_estimator, parameters).set_callbacks(scoring)
    >>> meta_estimator.fit()
    >>> scoring.get_logs()
    #TODO add printed output, or maybe use the logs in a plot

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
