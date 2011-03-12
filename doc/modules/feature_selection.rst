
.. _feature_selection_doc:

=================
Feature selection
=================



Univariate feature selection
============================

Univariate feature selection works by selecting the best features based on
univariate statistical tests. It can seen as a preprocessing step
to an estimator. The `scikit.learn` exposes feature selection routines
a objects that implement the `transform` method. The k-best features
can be selected based on:

.. autofunction:: scikits.learn.feature_selection.univariate_selection.SelectKBest

or by setting a percentile of features to keep using

.. autofunction:: scikits.learn.feature_selection.univariate_selection.SelectPercentile

or using common statistical quantities:

.. autofunction:: scikits.learn.feature_selection.univariate_selection.SelectFpr
.. autofunction:: scikits.learn.feature_selection.univariate_selection.SelectFdr
.. autofunction:: scikits.learn.feature_selection.univariate_selection.SelectFwe

These objects take as input a scoring function that returns
univariate p-values.

.. topic:: Examples:

    :ref:`example_plot_feature_selection.py`


Feature scoring functions
-------------------------

.. warning::

    Beware not to use a regression scoring function with a classification problem.

For classification
..................

.. autofunction:: scikits.learn.feature_selection.univariate_selection.f_classif

For regression
..............

.. autofunction:: scikits.learn.feature_selection.univariate_selection.f_regression


