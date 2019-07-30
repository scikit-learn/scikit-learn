

.. currentmodule:: sklearn.model_selection

.. _grid_search:

===========================================
Tuning the hyper-parameters of an estimator
===========================================

Hyper-parameters are parameters that are not directly learnt within estimators.
In scikit-learn they are passed as arguments to the constructor of the
estimator classes. Typical examples include ``C``, ``kernel`` and ``gamma``
for Support Vector Classifier, ``alpha`` for Lasso, etc.

It is possible and recommended to search the hyper-parameter space for the
best :ref:`cross validation <cross_validation>` score.

Any parameter provided when constructing an estimator may be optimized in this
manner. Specifically, to find the names and current values for all parameters
for a given estimator, use::

  estimator.get_params()

A search consists of:

- an estimator (regressor or classifier such as ``sklearn.svm.SVC()``);
- a parameter space;
- a method for searching or sampling candidates;
- a cross-validation scheme; and
- a :ref:`score function <gridsearch_scoring>`.

Two generic approaches to parameter search are provided in
scikit-learn: for given values, :class:`GridSearchCV` exhaustively considers
all parameter combinations, while :class:`RandomizedSearchCV` can sample a
given number of candidates from a parameter space with a specified
distribution. Both these tools have successive halving counterparts
:class:`GridHalvingSearchCV` and :class:`RandomHalvingSearchCV`, which can be
much faster at finding a good parameter combination.

After describing these tools we detail, :ref:`best practices
<grid_search_tips>` applicable to these approaches. Some models allow for
specialized, efficient parameter search strategies, outlined in
:ref:`alternative_cv`.

Note that it is common that a small subset of those parameters can have a large
impact on the predictive or computation performance of the model while others
can be left to their default values. It is recommended to read the docstring of
the estimator class to get a finer understanding of their expected behavior,
possibly by reading the enclosed reference to the literature.  

Exhaustive Grid Search
======================

The grid search provided by :class:`GridSearchCV` exhaustively generates
candidates from a grid of parameter values specified with the ``param_grid``
parameter. For instance, the following ``param_grid``::

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

.. currentmodule:: sklearn.model_selection

.. topic:: Examples:

    - See :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_digits.py` for an example of
      Grid Search computation on the digits dataset.

    - See :ref:`sphx_glr_auto_examples_model_selection_grid_search_text_feature_extraction.py` for an example
      of Grid Search coupling parameters from a text documents feature
      extractor (n-gram count vectorizer and TF-IDF transformer) with a
      classifier (here a linear SVM trained with SGD with either elastic
      net or L2 penalty) using a :class:`pipeline.Pipeline` instance.

    - See :ref:`sphx_glr_auto_examples_model_selection_plot_nested_cross_validation_iris.py`
      for an example of Grid Search within a cross validation loop on the iris
      dataset. This is the best practice for evaluating the performance of a
      model with grid search.

    - See :ref:`sphx_glr_auto_examples_model_selection_plot_multi_metric_evaluation.py`
      for an example of :class:`GridSearchCV` being used to evaluate multiple
      metrics simultaneously.

    - See :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_refit_callable.py`
      for an example of using ``refit=callable`` interface in
      :class:`GridSearchCV`. The example shows how this interface adds certain
      amount of flexibility in identifying the "best" estimator. This interface
      can also be used in multiple metrics evaluation.

.. _randomized_parameter_search:

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

  {'C': scipy.stats.expon(scale=100), 'gamma': scipy.stats.expon(scale=.1),
    'kernel': ['rbf'], 'class_weight':['balanced', None]}

This example uses the ``scipy.stats`` module, which contains many useful
distributions for sampling parameters, such as ``expon``, ``gamma``,
``uniform`` or ``randint``.
In principle, any function can be passed that provides a ``rvs`` (random
variate sample) method to sample a value. A call to the ``rvs`` function should
provide independent random samples from possible parameter values on
consecutive calls.

    .. warning::

        The distributions in ``scipy.stats`` prior to version scipy 0.16
        do not allow specifying a random state. Instead, they use the global
        numpy random state, that can be seeded via ``np.random.seed`` or set
        using ``np.random.set_state``. However, beginning scikit-learn 0.18,
        the :mod:`sklearn.model_selection` module sets the random state provided
        by the user if scipy >= 0.16 is also available.

For continuous parameters, such as ``C`` above, it is important to specify
a continuous distribution to take full advantage of the randomization. This way,
increasing ``n_iter`` will always lead to a finer search.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_model_selection_plot_randomized_search.py` compares the usage and efficiency
      of randomized search and grid search.

.. topic:: References:

    * Bergstra, J. and Bengio, Y.,
      Random search for hyper-parameter optimization,
      The Journal of Machine Learning Research (2012)

.. _successive_halving_user_guide:

Searching optimal parameters with successive halving
====================================================

Scikit-learn also provides the :class:`GridHalvingSearchCV` and
:class:`RandomHalvingSearchCV` estimators that can be used to
search a parameter space using successive halving [1]_ [2]_. Successive
halving is an iterative selection process where all candidates are evaluated
with a small amount of resources at the first iteration. Only some of
these candidates are selected for the next iteration, which will be
allocated more resources. What defines a resource is typically the number of
samples to train on, or the number of trees for a gradient boosting /
decision forest estimator.

As illustrated in the figure below, only a small subset of candidates 'survive'
until the last iteration. These are the candidates that have consistently been
part of the best candidates across all iterations.

.. figure:: ../auto_examples/svm/images/sphx_glr_plot_successive_halving_iterations_001.png
   :target: ../auto_examples/model_selection/plot_successive_halving_iterations.html
   :align: center

The amount of resources ``r_i`` allocated for each candidate at iteration
``i`` is controlled by the parameters ``ratio`` and ``r_min`` as follows::

    r_i = ratio**i * r_min

``r_min`` is the amount of resources used at the first iteration and
``ratio`` defines the proportions of candidates that will be selected for
the next iteration::

    n_candidates_to_keep = n_candidates_at_i // ratio

So in the first iteration, we use ``r_min`` resources ``n_candidates``
times. In the second iteration, we use ``r_min * ratio`` resources
``n_candidates // ratio`` times. The third again multiplies the resources
per candidate and divides the number of candidates. This process stops when
the maximum budget per candidate is reached, or when less than ``ratio``
candidates are left.

Here is an example with ``r_min=3`` and ``ratio=2``, starting with 70
candidates:

+-------------+-----------------------+
| ``r_i``     | ``n_candidates_at_i`` |
+=============+=======================+
| 3 (=r_min)  | 70 (=n_candidates)    |
+-------------+-----------------------+
| 3 * 2 = 6   | 70 // 2 = 35          |
+-------------+-----------------------+
| 6 * 2 = 12  | 35 // 2 = 17          |
+-------------+-----------------------+
| 12 * 2 = 24 | 17 // 2 = 8           |
+-------------+-----------------------+
| 24 * 2 = 48 | 8 // 2 = 4            |
+-------------+-----------------------+
| 48 * 2 = 96 | 4 // 2 = 2            |
+-------------+-----------------------+

At the last iteration, ``ratio`` candidates are evaluated, and we can pick
the best one. Note that each ``r_i`` is a multiple of both ``ratio`` and
``r_min``.

Choosing the budget
-------------------

By default, the budget is defined as the number of samples. That is, each
iteration will use an increasing amount of samples to train on. You can however
manually specify a parameter to use as the budget with the ``budget_on``
parameter. Here is an example where the budget is defined as the number of
iterations of a random forest::

    >>> from sklearn.datasets import make_classification
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.model_selection import GridHalvingSearchCV
    >>> import pandas as pd
    >>>
    >>> param_grid = {'max_depth': [3, 5, 10],
    ...               'min_samples_split': [2, 5, 10]}
    >>> base_estimator = RandomForestClassifier(random_state=0)
    >>> X, y = make_classification(n_samples=1000, random_state=0)
    >>> sh = GridHalvingSearchCV(base_estimator, param_grid, cv=5,
    ...                          ratio=2, budget_on='n_estimators',
    ...                          max_budget=30, random_state=0).fit(X, y)
    >>> sh.best_estimator_
    RandomForestClassifier(max_depth=5, n_estimators=8, random_state=0)

Note that it is not possible to budget on a parameter that is part of the
parameter grid.

Exhausting the budget
---------------------

As mentioned above, the first iteration uses ``r_min`` resources. If you have
a big budget, this may be a waste of resources::

    >>> from sklearn.datasets import make_classification
    >>> from sklearn.svm import SVC
    >>> from sklearn.model_selection import GridHalvingSearchCV
    >>> import pandas as pd
    >>> param_grid= {'kernel': ('linear', 'rbf'),
    ...              'C': [1, 10, 100]}
    >>> base_estimator = SVC(gamma='scale')
    >>> X, y = make_classification(n_samples=1000)
    >>> sh = GridHalvingSearchCV(base_estimator, param_grid, cv=5,
    ...                          ratio=2).fit(X, y)
    >>> results = pd.DataFrame(sh.cv_results_)
    >>> results.groupby('iter')['r_i'].unique()
    iter
    0    [20]
    1    [40]
    2    [80]
    Name: r_i, dtype: object

The search process will only use 80 resources at most, while our maximum budget
is ``n_samples=1000``. Note in this case that ``r_min = r_0 = 20``. In order
for the last iteration to use as many resources as possible, you can use the
``force_exhaust_budget`` parameter.::

    >>> sh = GridHalvingSearchCV(base_estimator, param_grid, cv=5,
    ...                            ratio=2, force_exhaust_budget=True,
    ...                            ).fit(X, y)
    >>> results = pd.DataFrame.from_dict(sh.cv_results_)
    >>> results.groupby('iter')['r_i'].unique()
    iter
    0     [250]
    1     [500]
    2    [1000]
    Name: r_i, dtype: object


`r_min` was here automatically set to 250, which results in the last
iteration using all the budget. Since ``force_exhaust_budget`` chooses an
appropriate ``r_min`` to start with, ``r_min`` must be set to 'auto' (default).

Aggressive elimination of candidates
------------------------------------

Ideally, we want the last iteration to evaluate ``ratio`` candidates. We then
just have to pick the best one. When the budget is small with respect to
the number of candidates, the last iteration may have to evaluate more than
``ratio`` candidates::
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.svm import SVC
    >>> from sklearn.model_selection import GridHalvingSearchCV
    >>> import pandas as pd
    >>>
    >>>
    >>> param_grid = {'kernel': ('linear', 'rbf'),
    ...               'C': [1, 10, 100]}
    >>> base_estimator = SVC(gamma='scale')
    >>> X, y = make_classification(n_samples=1000)
    >>> sh = GridHalvingSearchCV(base_estimator, param_grid, cv=5,
    ...                          ratio=2, max_budget=40,
    ...                          aggressive_elimination=False).fit(X, y)
    >>> results = pd.DataFrame.from_dict(sh.cv_results_)
    >>> results.groupby('iter').r_i.unique()
    iter
    0    [20]
    1    [40]
    Name: r_i, dtype: object
    >>> results.groupby('iter').r_i.count()  # number of candidates used at each iteration
    iter
    0    6
    1    3
    Name: r_i, dtype: int64

Since we cannot use more than ``max_budget=40`` resources, the process has to
stop at the second iteration which evaluates more than ``ratio=2`` candidates.

Using the ``aggressive_elimination`` parameter, you can force the search
process to end up with less than ``ratio`` candidates at the last
iteration. To do this, the process will eliminate as many candidates as
necessary using ``r_min`` resources::

    >>> sh = GridHalvingSearchCV(base_estimator, param_grid, cv=5,
    ...                            ratio=2,
    ...                            max_budget=40,
    ...                            aggressive_elimination=True,
    ...                            ).fit(X, y)
    >>> results = pd.DataFrame.from_dict(sh.cv_results_)
    >>> results.groupby('iter').r_i.unique()
    iter
    0    [20]
    1    [20]
    2    [40]
    Name: r_i, dtype: object
    >>> results.groupby('iter').r_i.count()  # number of candidates used at each iteration
    iter
    0    6
    1    3
    2    2
    Name: r_i, dtype: int64

Notice that we end with 2 candidates at the last iteration since we have
eliminated enough candidates during the first iterations, using ``r_i = r_min =
20``.


.. topic:: References:

    .. [1] K. Jamieson, A. Talwalkar,
       `Non-stochastic Best Arm Identification and Hyperparameter
       Optimization <http://proceedings.mlr.press/v51/jamieson16.html>`_, in
       proc. of Machine Learning Research, 2016.
    .. [2] L. Li, K. Jamieson, G. DeSalvo, A. Rostamizadeh, .A Talwalkar,
       `Hyperband: A Novel Bandit-Based Approach to Hyperparameter Optimization
       <https://arxiv.org/abs/1603.06560>`_, in Machine Learning Research
       18, 2018.

.. _grid_search_tips:

Tips for parameter search
=========================

.. _gridsearch_scoring:

Specifying an objective metric
------------------------------

By default, parameter search uses the ``score`` function of the estimator
to evaluate a parameter setting. These are the
:func:`sklearn.metrics.accuracy_score` for classification and
:func:`sklearn.metrics.r2_score` for regression.  For some applications,
other scoring functions are better suited (for example in unbalanced
classification, the accuracy score is often uninformative). An alternative
scoring function can be specified via the ``scoring`` parameter of most
parameter search tools. See :ref:`scoring_parameter` for more details.

.. _multimetric_grid_search:

Specifying multiple metrics for evaluation
------------------------------------------

:class:`GridSearchCV` and :class:`RandomizedSearchCV` allow specifying
multiple metrics for the ``scoring`` parameter.

Multimetric scoring can either be specified as a list of strings of predefined
scores names or a dict mapping the scorer name to the scorer function and/or
the predefined scorer name(s). See :ref:`multimetric_scoring` for more details.

When specifying multiple metrics, the ``refit`` parameter must be set to the
metric (string) for which the ``best_params_`` will be found and used to build
the ``best_estimator_`` on the whole dataset. If the search should not be
refit, set ``refit=False``. Leaving refit to the default value ``None`` will
result in an error when using multiple metrics.

See :ref:`sphx_glr_auto_examples_model_selection_plot_multi_metric_evaluation.py`
for an example usage.

Composite estimators and parameter spaces
-----------------------------------------

Please refer to :ref:`pipeline` for performing parameter searches over
pipelines.

Model selection: development and evaluation
-------------------------------------------

Model selection by evaluating various parameter settings can be seen as a way
to use the labeled data to "train" the parameters of the grid.

When evaluating the resulting model it is important to do it on
held-out samples that were not seen during the grid search process:
it is recommended to split the data into a **development set** (to
be fed to the :class:`GridSearchCV` instance) and an **evaluation set**
to compute performance metrics.

This can be done by using the :func:`train_test_split`
utility function.

Parallelism
-----------

The parameter search tools evaluate each parameter setting independently.
Computations can be run in parallel if your OS supports it, by using the
keyword ``n_jobs=-1``. See function signature for more details.

Robustness to failure
---------------------

Some parameter settings may result in a failure to ``fit`` one or more folds
of the data.  By default, this will cause the entire search to fail, even if
some parameter settings could be fully evaluated. Setting ``error_score=0``
(or `=np.NaN`) will make the procedure robust to such failure, issuing a
warning and setting the score for that fold to 0 (or `NaN`), but completing
the search.

.. _alternative_cv:

Alternatives to brute force parameter search
============================================

Model specific cross-validation
-------------------------------


Some models can fit data for a range of values of some parameter almost
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

   linear_model.ElasticNetCV
   linear_model.LarsCV
   linear_model.LassoCV
   linear_model.LassoLarsCV
   linear_model.LogisticRegressionCV
   linear_model.MultiTaskElasticNetCV
   linear_model.MultiTaskLassoCV
   linear_model.OrthogonalMatchingPursuitCV
   linear_model.RidgeCV
   linear_model.RidgeClassifierCV


Information Criterion
---------------------

Some models can offer an information-theoretic closed-form formula of the
optimal estimate of the regularization parameter by computing a single
regularization path (instead of several when using cross-validation).

Here is the list of models benefiting from the Akaike Information
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
