.. include:: includes/big_toc_css.rst

.. currentmodule:: sklearn

.. _model_evaluation:

========================================================
Model evaluation: quantifying the quality of predictions
========================================================

There are 3 different APIs for evaluating the quality of a model's
predictions:

* **Estimator score method**: Estimators have a ``score`` method providing a
  default evaluation criterion for the problem they are designed to solve.
  This is not discussed in this section, but in each estimator's documentation.

  Here is an example demonstrating the ``score`` method for Linear Regression:

    >>> from sklearn import linear_model
    >>> from sklearn.model_selection import train_test_split
    >>> est = linear_model.LinearRegression()
    >>> X = [[1, 1], [2, 2], [4, 4], [3, 3],[8.75, 8.65], [10, 10], [5, 5.5]]
    >>> y = [1, 2, 4, 3, 8.7, 10, 5.3]
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    >>> est.fit (X_train, y_train)
    LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1, normalize=False)
    >>> est.score(X_test, y_test)
    0.997...


* **Scoring parameter**: Model-evaluation tools using
  :ref:`cross-validation <cross_validation>` (such as
  :func:`model_selection.cross_val_score` and
  :class:`model_selection.GridSearchCV`) rely on an internal *scoring* strategy.
  This is discussed in the section :ref:`scoring_parameter`.

  This is an example of :func:`model_selection.cross_val_score` on the same
  data as above:

    >>> import numpy as np
    >>> from sklearn.model_selection import cross_val_score
    >>> scores = cross_val_score(est, X, y)
    >>> np.mean(scores)
    0.99992202084853343


* **Metric functions**: The :mod:`metrics` module implements functions
  assessing prediction error for specific purposes. These metrics are detailed
  in sections on

  * :ref:`regression_metrics`
  * :ref:`classification_metrics`
  * :ref:`multilabel_ranking_metrics`
  * :ref:`clustering_metrics`
  * :ref:`dummy_estimators`


  Here are examples of two different classification
  metrics :func:`metrics.accuracy_score` and :func:`metrics.roc_auc_score`:

    >>> from sklearn.metrics import accuracy_score
    >>> y_pred = [0, 2, 1, 3]
    >>> y_true = [0, 1, 2, 3]
    >>> accuracy_score(y_true, y_pred)
    0.5


    >>> from sklearn.metrics import roc_auc_score
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> roc_auc_score(y_true, y_scores)
    0.75


As seen with the above examples we can only use the ``score`` method on
classifiers, but the other two are more flexible, and the **scoring**
parameter is probably the easiest to use.


.. seealso::

   For "pairwise" metrics, between *samples* and not estimators or
   predictions, see the :ref:`metrics` section.

.. toctree::
    :hidden:

    ./model_evaluation_scoring
    ./model_evaluation_classification
    ./model_evaluation_multilabel
    ./model_evaluation_regression
    ./model_evaluation_clustering
    ./model_evaluation_dummy

