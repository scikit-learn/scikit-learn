.. _grid_search:

.. currentmodule:: sklearn.grid_search

===============================================
Grid Search: Searching for estimator parameters
===============================================

Parameters that are not directly learnt within estimators can be set by
searching a parameter space for the best :ref:`cross_validation` score.
Typical examples include ``C``, ``kernel`` and ``gamma`` for Support Vector
Classifier, ``alpha`` for Lasso, etc.

Any parameter provided when constructing an estimator may be optimized in this
manner.  Specifically, to find the names and current values for all parameters
for a given estimator, use::

  estimator.get_params()

Such parameters are often referred to as *hyperparameters* (particularly in
Bayesian learning), distinguishing them from the parameters optimised in a
machine learning procedure.

A search consists of:

- an estimator (regressor or classifier such as ``sklearn.svm.SVC()``);
- a parameter space;
- a method for searching or sampling candidates;
- a cross-validation scheme; and
- a :ref:`score function <gridsearch_scoring>`.

Two generic approaches to sampling search candidates are provided in
scikit-learn: for given values, :class:`GridSearchCV` exhaustively considers
all parameter combinations, while :class:`RandomizedSearchCV` can sample a
given number of candidates from a parameter space with a specified
distribution.

Exhaustive Grid Search
======================

The grid search provided by :class:`GridSearchCV` exhaustively generates
candidates from a grid of parameter values specified with the `param_grid`
attribute. For instance, the following `param_grid`::

  param_grid = [
    {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
    {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
   ]

specifies that two grids should be explored: one with a linear kernel and
C values in [1, 10, 100, 1000], and the second one with an RBF kernel,
and the cross-product of C values ranging in [1, 10, 100, 1000] and gamma
values in [0.001, 0.0001].

The :class:`GridSearchCV` instance implements the usual estimator API: when
"fitting" it on a dataset all the possible combinations of parameter values are
evaluated and the best combination is retained.

.. topic:: Model selection: development and evaluation

  Model selection with ``GridSearchCV`` can be seen as a way to use the
  labeled data to "train" the parameters of the grid.

  When evaluating the resulting model it is important to do it on
  held-out samples that were not seen during the grid search process:
  it is recommended to split the data into a **development set** (to
  be fed to the ``GridSearchCV`` instance) and an **evaluation set**
  to compute performance metrics.

  This can be done by using the :func:`cross_validation.train_test_split`
  utility function.

.. currentmodule:: sklearn.grid_search

.. _gridsearch_scoring:

Scoring functions for parameter search
--------------------------------------

By default, :class:`GridSearchCV` uses the ``score`` function of the estimator
to evaluate a parameter setting. These are the
:func:`sklearn.metrics.accuracy_score` for classification and
:func:`sklearn.metrics.r2_score` for regression.  For some applications, other
scoring functions are better suited (for example in unbalanced classification,
the accuracy score is often uninformative). An alternative scoring function
can be specified via the ``scoring`` parameter to :class:`GridSearchCV`.  See
:ref:`scoring_parameter` for more details.

.. topic:: Examples:

    - See :ref:`example_grid_search_digits.py` for an example of
      Grid Search computation on the digits dataset.

    - See :ref:`example_grid_search_text_feature_extraction.py` for an example
      of Grid Search coupling parameters from a text documents feature
      extractor (n-gram count vectorizer and TF-IDF transformer) with a
      classifier (here a linear SVM trained with SGD with either elastic
      net or L2 penalty) using a :class:`pipeline.Pipeline` instance.

.. note::

  Computations can be run in parallel if your OS supports it, by using
  the keyword ``n_jobs=-1``, see function signature for more details.


Randomized Parameter Optimization
=================================
While using a grid of parameter settings is currently the most widely used
method for parameter optimization, other search methods have more
favourable properties.
:class:`RandomizedSearchCV` implements a randomized search over parameters,
where each setting is sampled from a distribution over possible parameter values.
This has two main benefits over an exhaustive search:

* A budget can be chosen independent of the number of parameters and possible values.
* Adding parameters that do not influence the performance does not decrease efficiency.

Specifying how parameters should be sampled is done using a dictionary, very
similar to specifying parameters for :class:`GridSearchCV`. Additionally,
a computation budget, being the number of sampled candidates or sampling
iterations, is specified using the ``n_iter`` parameter.
For each parameter, either a distribution over possible values or a list of
discrete choices (which will be sampled uniformly) can be specified::

  [{'C': scipy.stats.expon(scale=100), 'gamma': scipy.stats.expon(scale=.1),
    'kernel': ['rbf'], 'class_weight':['auto', None]}]

This example uses the ``scipy.stats`` module, which contains many useful
distributions for sampling parameters, such as ``expon``, ``gamma``,
``uniform`` or ``randint``.
In principle, any function can be passed that provides a ``rvs`` (random
variate sample) method to sample a value. A call to the ``rvs`` function should
provide independent random samples from possible parameter values on
consecutive calls.

    .. warning::

        The distributions in ``scipy.stats`` do not allow specifying a random
        state. Instead, they use the global numpy random state, that can be seeded
        via ``np.random.seed`` or set using ``np.random.set_state``.

For continuous parameters, such as ``C`` above, it is important to specify
a continuous distribution to take full advantage of the randomization. This way,
increasing ``n_iter`` will always lead to a finer search.

.. topic:: Examples:

    * :ref:`example_randomized_search.py` compares the usage and efficiency
      of randomized search and grid search.

.. topic:: References:

    * Bergstra, J. and Bengio, Y.,
      Random search for hyper-parameter optimization,
      The Journal of Machine Learning Research (2012)


Alternatives to brute force parameter search
============================================

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

.. currentmodule:: sklearn

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
comes "for free" as no additional data is needed and can be used for
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
