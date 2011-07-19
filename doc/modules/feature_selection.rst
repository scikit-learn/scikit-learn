
.. _feature_selection_doc:

=================
Feature selection
=================

The classes in the ``scikits.learn.feature_selection`` module can be used
for feature selection/dimensionality reduction on sample sets, either to
improve estimators' accuracy scores or to boost their performance on very
high-dimensional datasets.


Ï<87>Â² (chi squared) feature selection
=======================================

The :class:`Chi2` transformer class will pass through a user-selected number
of features from a vector of samples, selecting those for which the Ï<87>Â²
statistic yields the highest degree of independence from their target
values/classes.

:class:`Chi2` must be ``fit`` on a training set before use. It is designed for
binary or multinomial data, i.e. arrays of either boolean
occurrence/non-occurrence indicators or occurence counts for each sample.
It handles ``scipy.sparse`` matrices as well as instances of ``numpy.array``.


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


