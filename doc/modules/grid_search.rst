.. _grid_search:

==========================================
Grid Search: setting estimator parameters
==========================================

.. currentmodule:: sklearn

Grid Search is used to optimize the parameters of a model (e.g. ``C``,
``kernel`` and ``gamma`` for Support Vector Classifier, ``alpha`` for
Lasso, etc.) using an internal :ref:`cross_validation` scheme).


GridSearchCV
============

The main class for implementing hyperparameters grid search in
scikit-learn is :class:`grid_search.GridSearchCV`. This class is passed
a base model instance (for example ``sklearn.svm.SVC()``) along with a
grid of potential hyper-parameter values such as::

  [{'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
   {'C': [1, 10, 100, 1000], 'kernel': ['linear']}]

The :class:`grid_search.GridSearchCV` instance implements the usual
estimator API: when "fitting" it on a dataset all the possible
combinations of hyperparameter values are evaluated and the best
combinations is retained.

.. topic:: Model selection: development and evaluation

  Model selection with ``GridSearchCV`` can be seen as a way to use the
  labeled data to "train" the hyper-parameters of the grid.

  When evaluating the resulting model it is important to do it on
  held-out samples that were not seen during the grid search process:
  it is recommended to split the data into a **development set** (to
  be fed to the ``GridSearchCV`` instance) and an **evaluation set**
  to compute performance metrics.

  This can be done by using the :func:`cross_validation.train_test_split`
  utility function.


Examples
========

- See :ref:`example_grid_search_digits.py` for an example of
  Grid Search computation on the digits dataset.

- See :ref:`example_grid_search_text_feature_extraction.py` for an example
  of Grid Search coupling parameters from a text documents feature
  extractor (n-gram count vectorizer and TF-IDF transformer) with a
  classifier (here a linear SVM trained with SGD with either elastic
  net or L2 penalty) using a :class:`pipeline.Pipeline` instance.

.. note::

  Computations can be run in parallel if your OS supports it, by using
  the keyword n_jobs=-1, see function signature for more details.


Alternatives to brute force grid search
=======================================

Model specific cross-validation
-------------------------------


Some models can fit data for a range of value of some parameter almost
as efficiently as fitting the estimator for a single value of the
parameter. This feature can be leveraged to perform a more efficient
cross-validation used for model selection of this parameter.

The most common parameter amenable to this strategy is the parameter
encoding the strength of the regularizer. In this case we say that we
compute the **regularization path** of the estimator.

Here is the list of such models:

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.RidgeCV
   linear_model.RidgeClassifierCV
   linear_model.LarsCV
   linear_model.LassoLarsCV
   linear_model.LassoCV
   linear_model.ElasticNetCV


Information Criterion
---------------------

Some models can offer an information-theoretic closed-form formula of the
optimal estimate of the regularization parameter by computing a single
regularization path (instead of several when using cross-validation).

Here is the list of models benefitting from the Aikike Information
Criterion (AIC) or the Bayesian Information Criterion (BIC) for automated
model selection:

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.LassoLarsIC


.. _out_of_bag:

Out of Bag Estimates
--------------------

When using ensemble methods base upon bagging, i.e. generating new
training sets using sampling with replacement, part of the training set
remains unused.  For each classifier in the ensemble, a different part
of the training set is left out.

This left out portion can be used to estimate the generalization error
without having to rely on a separate validation set.  This estimate
comes "for free" as no addictional data is needed and can be used for
model selection.

This is currently implemented in the following classes:

.. autosummary::
   :toctree: generated/
   :template: class.rst

    ensemble.RandomForestClassifier
    ensemble.RandomForestRegressor
    ensemble.ExtraTreesClassifier
    ensemble.ExtraTreesRegressor
    ensemble.GradientBoostingClassifier
    ensemble.GradientBoostingRegressor
