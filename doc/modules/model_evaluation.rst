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
  This is not discussed on this page, but in each estimator's documentation.

* **Scoring parameter**: Model-evaluation tools using
  :ref:`cross-validation <cross_validation>` (such as
  :func:`model_selection.cross_val_score` and
  :class:`model_selection.GridSearchCV`) rely on an internal *scoring* strategy.
  This is discussed in the section :ref:`scoring_parameter`.

* **Metric functions**: The :mod:`metrics` module implements functions
  assessing prediction error for specific purposes. These metrics are detailed
  in sections on :ref:`classification_metrics`,
  :ref:`multilabel_ranking_metrics`, :ref:`regression_metrics` and
  :ref:`clustering_metrics`.

Finally, :ref:`dummy_estimators` are useful to get a baseline
value of those metrics for random predictions.

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

