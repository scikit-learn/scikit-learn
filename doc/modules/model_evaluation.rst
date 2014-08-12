.. currentmodule:: sklearn

.. _model_evaluation:

========================================================
Model evaluation: quantifying the quality of predictions
========================================================

There are 3 different approaches to evaluate the quality of predictions of a
model:

* **Estimator score method**: Estimators have a ``score`` method providing a
  default evaluation criterion for the problem they are designed to solve.
  This is not discussed on this page, but in each estimator's documentation.

* **Scoring parameter**: Model-evaluation tools using
  :ref:`cross-validation <cross_validation>` (such as
  :func:`cross_validation.cross_val_score` and
  :class:`grid_search.GridSearchCV`) rely on an internal *scoring* strategy.
  This is discussed on section :ref:`scoring_parameter`.

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

.. _scoring_parameter:

The ``scoring`` parameter: defining model evaluation rules
==========================================================

Model selection and evaluation using tools, such as
:class:`grid_search.GridSearchCV` and
:func:`cross_validation.cross_val_score`, take a ``scoring`` parameter that
controls what metric they apply to estimators evaluated.

Common cases: predefined values
-------------------------------

For the most common usecases, you can simply provide a string as the
``scoring`` parameter. Possible values are:

======================     =================================================
Scoring                    Function
======================     =================================================
**Classification**
'accuracy'                 :func:`sklearn.metrics.accuracy_score`
'average_precision'        :func:`sklearn.metrics.average_precision_score`
'f1'                       :func:`sklearn.metrics.f1_score`
'precision'                :func:`sklearn.metrics.precision_score`
'recall'                   :func:`sklearn.metrics.recall_score`
'roc_auc'                  :func:`sklearn.metrics.roc_auc_score`

**Clustering**
'adjusted_rand_score'      :func:`sklearn.metrics.adjusted_rand_score`

**Regression**
'mean_absolute_error'      :func:`sklearn.metrics.mean_absolute_error`
'mean_squared_error'       :func:`sklearn.metrics.mean_squared_error`
'r2'                       :func:`sklearn.metrics.r2_score`
======================     =================================================

Setting the ``scoring`` parameter to a wrong value should give you a list
of acceptable values::

    >>> from sklearn import svm, cross_validation, datasets
    >>> iris = datasets.load_iris()
    >>> X, y = iris.data, iris.target
    >>> model = svm.SVC()
    >>> cross_validation.cross_val_score(model, X, y, scoring='wrong_choice')
    Traceback (most recent call last):
    ValueError: 'wrong_choice' is not a valid scoring value. Valid options are ['accuracy', 'adjusted_rand_score', 'average_precision', 'f1', 'log_loss', 'mean_absolute_error', 'mean_squared_error', 'precision', 'r2', 'recall', 'roc_auc']

.. note::

    The corresponding scorer objects are stored in the dictionary
    ``sklearn.metrics.SCORERS``.

The above choices correspond to error-metric functions that can be applied to
predicted values. These are detailed below, in the next sections.


.. currentmodule:: sklearn.metrics

.. _scoring:

Defining your scoring strategy from metric functions
-----------------------------------------------------

The module :mod:`sklearn.metric` also exposes a set of simple functions
measuring a prediction error given ground truth and prediction:

- functions ending with ``_score`` return a value to
  maximize (the higher the better).

- functions ending with ``_error`` or ``_loss`` return a
  value to minimize (the lower the better).

Metrics available for various machine learning tasks are detailed in sections
below.

Many metrics are not given names to be used as ``scoring`` values,
sometimes because they require additional parameters, such as
:func:`fbeta_score`. In such cases, you need to generate an appropriate
scoring object.  The simplest way to generate a callable object for scoring
is by using :func:`make_scorer`. That function converts metrics
into callables that can be used for model evaluation.

One typical use case is to wrap an existing metric function from the library
with non default value for its parameters, such as the ``beta`` parameter for
the :func:`fbeta_score` function::

    >>> from sklearn.metrics import fbeta_score, make_scorer
    >>> ftwo_scorer = make_scorer(fbeta_score, beta=2)
    >>> from sklearn.grid_search import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]}, scoring=ftwo_scorer)

The second use case is to build a completely new and custom scorer object
from a simple python function::

    >>> def my_custom_loss_func(ground_truth, predictions):
    ...     diff = np.abs(ground_truth - predictions).max()
    ...     return np.log(1 + diff)
    ...
    >>> my_custom_scorer = make_scorer(my_custom_loss_func, greater_is_better=False)
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]}, scoring=my_custom_scorer)

:func:`make_scorer` takes as parameters:

* the function you want to use

* whether it is a score (``greater_is_better=True``) or a loss
  (``greater_is_better=False``),

* whether the function you provided takes predictions as input
  (``needs_threshold=False``) or needs confidence scores \
  (``needs_threshold=True``)

* any additional parameters, such as ``beta`` in an :func:`f1_score`.


.. _diy_scoring:

Implementing your own scoring object
------------------------------------
You can generate even more flexible model scores by constructing your own
scoring object from scratch, without using the :func:`make_scorer` factory.
For a callable to be a scorer, it needs to meet the protocol specified by
the following two rules:

- It can be called with parameters ``(estimator, X, y)``, where ``estimator``
  is the model that should be evaluated, ``X`` is validation data, and ``y`` is
  the ground truth target for ``X`` (in the supervised case) or ``None`` (in the
  unsupervised case).

- It returns a floating point number that quantifies the quality of
  ``estimator``'s predictions on ``X`` which reference to ``y``.
  Again, higher numbers are better.

.. _classification_metrics:

Classification metrics
=======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
functions to measure classification performance.
Some metrics might require probability estimates of the positive class,
confidence values or binary decisions values.

Some of these are restricted to the binary classification case:

.. autosummary::
   :template: function.rst

   hinge_loss
   matthews_corrcoef
   precision_recall_curve
   roc_curve


Others also work in the multiclass case:

.. autosummary::
   :template: function.rst

   confusion_matrix


And some also work in the multilabel case:

.. autosummary::
   :template: function.rst

   accuracy_score
   classification_report
   f1_score
   fbeta_score
   hamming_loss
   jaccard_similarity_score
   log_loss
   precision_recall_fscore_support
   precision_score
   recall_score
   zero_one_loss

And some work with binary and multilabel indicator format:

.. autosummary::
   :template: function.rst

   average_precision_score
   roc_auc_score


In the following sub-sections, we will describe each of those functions.

Accuracy score
--------------

The :func:`accuracy_score` function computes the
`accuracy <http://en.wikipedia.org/wiki/Accuracy_and_precision>`_, the fraction
(default) or the number of correct predictions.

In multilabel classification, the function returns the subset accuracy: if
the entire set of predicted labels for a sample strictly match with the true
set of labels, then the subset accuracy is 1.0, otherwise it is 0.0.

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the fraction of correct predictions over :math:`n_\text{samples}` is
defined as

.. math::

   \texttt{accuracy}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples}-1} 1(\hat{y}_i = y_i)

where :math:`1(x)` is the `indicator function
<http://en.wikipedia.org/wiki/Indicator_function>`_.

  >>> import numpy as np
  >>> from sklearn.metrics import accuracy_score
  >>> y_pred = [0, 2, 1, 3]
  >>> y_true = [0, 1, 2, 3]
  >>> accuracy_score(y_true, y_pred)
  0.5
  >>> accuracy_score(y_true, y_pred, normalize=False)
  2

In the multilabel case with binary label indicators: ::

  >>> accuracy_score(np.array([[0, 1], [1, 1]]), np.ones((2, 2)))
  0.5

.. topic:: Example:

  * See :ref:`example_feature_selection_plot_permutation_test_for_classification.py`
    for an example of accuracy score usage using permutations of
    the dataset.

Confusion matrix
----------------

The :func:`confusion_matrix` function computes the `confusion matrix
<http://en.wikipedia.org/wiki/Confusion_matrix>`_ to evaluate
the accuracy on a classification problem.

By definition, a confusion matrix :math:`C` is such that :math:`C_{i, j}` is
equal to the number of observations known to be in group :math:`i` but
predicted to be in group :math:`j`. Here an example of such confusion matrix::

  >>> from sklearn.metrics import confusion_matrix
  >>> y_true = [2, 0, 2, 2, 0, 1]
  >>> y_pred = [0, 0, 2, 2, 0, 2]
  >>> confusion_matrix(y_true, y_pred)
  array([[2, 0, 0],
         [0, 0, 1],
         [1, 0, 2]])

Here a visual representation of such confusion matrix (this figure comes
from the :ref:`example_model_selection_plot_confusion_matrix.py` example):

.. image:: ../auto_examples/model_selection/images/plot_confusion_matrix_001.png
   :target: ../auto_examples/model_selection/plot_confusion_matrix.html
   :scale: 75
   :align: center

.. topic:: Example:

  * See :ref:`example_model_selection_plot_confusion_matrix.py`
    for an example of confusion matrix usage to evaluate the quality of the
    output of a classifier.

  * See :ref:`example_classification_plot_digits_classification.py`
    for an example of confusion matrix usage in the classification of
    hand-written digits.

  * See :ref:`example_text_document_classification_20newsgroups.py`
    for an example of confusion matrix usage in the classification of text
    documents.


Classification report
----------------------

The :func:`classification_report` function builds a text report showing the
main classification metrics. Here a small example with custom ``target_names``
and inferred labels::

   >>> from sklearn.metrics import classification_report
   >>> y_true = [0, 1, 2, 2, 0]
   >>> y_pred = [0, 0, 2, 2, 0]
   >>> target_names = ['class 0', 'class 1', 'class 2']
   >>> print(classification_report(y_true, y_pred, target_names=target_names))
                precision    recall  f1-score   support
   <BLANKLINE>
       class 0       0.67      1.00      0.80         2
       class 1       0.00      0.00      0.00         1
       class 2       1.00      1.00      1.00         2
   <BLANKLINE>
   avg / total       0.67      0.80      0.72         5
   <BLANKLINE>

.. topic:: Example:

  * See :ref:`example_classification_plot_digits_classification.py`
    for an example of classification report usage in the classification of the
    hand-written digits.

  * See :ref:`example_text_document_classification_20newsgroups.py`
    for an example of classification report usage in the classification of text
    documents.

  * See :ref:`example_model_selection_grid_search_digits.py`
    for an example of classification report usage in parameter estimation using
    grid search with a nested cross-validation.

Hamming loss
-------------

The :func:`hamming_loss` computes the average Hamming loss or `Hamming
distance <http://en.wikipedia.org/wiki/Hamming_distance>`_ between two sets
of samples.

If :math:`\hat{y}_j` is the predicted value for the :math:`j`-th labels of
a given sample, :math:`y_j` is the corresponding true value and
:math:`n_\text{labels}` is the number of class or labels, then the
Hamming loss :math:`L_{Hamming}` between two samples is defined as:

.. math::

   L_{Hamming}(y, \hat{y}) = \frac{1}{n_\text{labels}} \sum_{j=0}^{n_\text{labels} - 1} 1(\hat{y}_j \not= y_j)

where :math:`1(x)` is the `indicator function
<http://en.wikipedia.org/wiki/Indicator_function>`_. ::

  >>> from sklearn.metrics import hamming_loss
  >>> y_pred = [1, 2, 3, 4]
  >>> y_true = [2, 2, 3, 4]
  >>> hamming_loss(y_true, y_pred)
  0.25

In the multilabel case with binary label indicators: ::

  >>> hamming_loss(np.array([[0, 1], [1, 1]]), np.zeros((2, 2)))
  0.75

.. note::

    In multiclass classification, the Hamming loss correspond to the Hamming
    distance between ``y_true`` and ``y_pred`` which is equivalent to the
    :ref:`zero_one_loss` function.

    In multilabel classification, the Hamming loss is different from the
    zero-one loss. The zero-one loss penalizes any predictions that don't
    exactly match the true required set of labels,
    while Hamming loss will penalize the individual labels.
    So, predicting a subset or superset of the true labels
    will give a Hamming loss strictly between zero and one.

    The Hamming loss is upperbounded by the zero-one loss. When normalized
    over samples, the Hamming loss is always between zero and one.


Jaccard similarity coefficient score
-------------------------------------

The :func:`jaccard_similarity_score` function computes the average (default)
or sum of `Jaccard similarity coefficients
<http://en.wikipedia.org/wiki/Jaccard_index>`_, also called Jaccard index,
between pairs of label sets.

The Jaccard similarity coefficient of the :math:`i`-th samples
with a ground truth label set :math:`y_i` and a predicted label set
:math:`\hat{y}_i`  is defined as

.. math::

    J(y_i, \hat{y}_i) = \frac{|y_i \cap \hat{y}_i|}{|y_i \cup \hat{y}_i|}.

In binary and multiclass classification, the Jaccard similarity coefficient
score is equal to the classification accuracy.

::

  >>> import numpy as np
  >>> from sklearn.metrics import jaccard_similarity_score
  >>> y_pred = [0, 2, 1, 3]
  >>> y_true = [0, 1, 2, 3]
  >>> jaccard_similarity_score(y_true, y_pred)
  0.5
  >>> jaccard_similarity_score(y_true, y_pred, normalize=False)
  2

In the multilabel case with binary label indicators: ::

  >>> jaccard_similarity_score(np.array([[0, 1], [1, 1]]), np.ones((2, 2)))
  0.75

.. _precision_recall_f_measure_metrics:

Precision, recall and F-measures
---------------------------------

The `precision <http://en.wikipedia.org/wiki/Precision_and_recall#Precision>`_
is intuitively the ability of the classifier not to label as
positive a sample that is negative.

The `recall <http://en.wikipedia.org/wiki/Precision_and_recall#Recall>`_ is
intuitively the ability of the classifier to find all the positive samples.

The  `F-measure <http://en.wikipedia.org/wiki/F1_score>`_
(:math:`F_\beta` and :math:`F_1` measures) can be interpreted as a weighted
harmonic mean of the precision and recall. A
:math:`F_\beta` measure reaches its best value at 1 and worst score at 0.
With :math:`\beta = 1`, the :math:`F_\beta` measure leads to the
:math:`F_1` measure, wheres the recall and the precision are equally important.

The :func:`precision_recall_curve` computes a precision-recall curve
from the ground truth label and a score given by the classifier
by varying a decision threshold.

The :func:`average_precision_score` function computes the average precision
(AP) from prediction scores. This score corresponds to the area under the
precision-recall curve.

Several functions allow you to analyze the precision, recall and F-measures
score:

.. autosummary::
   :template: function.rst

   average_precision_score
   f1_score
   fbeta_score
   precision_recall_curve
   precision_recall_fscore_support
   precision_score
   recall_score

Note that the :func:`precision_recall_curve` function is restricted to the
binary case. The :func:`average_precision_score` function works only in
binary classification and multilabel indicator format.


.. topic:: Examples:

  * See :ref:`example_text_document_classification_20newsgroups.py`
    for an example of :func:`f1_score` usage with classification of text
    documents.

  * See :ref:`example_model_selection_grid_search_digits.py`
    for an example of :func:`precision_score` and :func:`recall_score` usage
    in parameter estimation using grid search with a nested cross-validation.

  * See :ref:`example_model_selection_plot_precision_recall.py`
    for an example of precision-Recall metric to evaluate the quality of the
    output of a classifier with :func:`precision_recall_curve`.

  * See :ref:`example_linear_model_plot_sparse_recovery.py`
    for an example of :func:`precision_recall_curve` usage in feature selection
    for sparse linear models.

Binary classification
^^^^^^^^^^^^^^^^^^^^^

In a binary classification task, the terms ''positive'' and ''negative'' refer
to the classifier's prediction and the terms ''true'' and ''false'' refer to
whether that prediction corresponds to the external judgment (sometimes known
as the ''observation''). Given these definitions, we can formulate the
following table:

+-------------------+------------------------------------------------+
|                   |    Actual class (observation)                  |
+-------------------+---------------------+--------------------------+
|   Predicted class | tp (true positive)  | fp (false positive)      |
|   (expectation)   | Correct result      | Unexpected result        |
|                   +---------------------+--------------------------+
|                   | fn (false negative) | tn (true negative)       |
|                   | Missing result      | Correct absence of result|
+-------------------+---------------------+--------------------------+

In this context, we can define the notions of precision, recall and F-measure:

.. math::

   \text{precision} = \frac{tp}{tp + fp},

.. math::

   \text{recall} = \frac{tp}{tp + fn},

.. math::

   F_\beta = (1 + \beta^2) \frac{\text{precision} \times \text{recall}}{\beta^2 \text{precision} + \text{recall}}.

Here some small examples in binary classification::

  >>> from sklearn import metrics
  >>> y_pred = [0, 1, 0, 0]
  >>> y_true = [0, 1, 0, 1]
  >>> metrics.precision_score(y_true, y_pred)
  1.0
  >>> metrics.recall_score(y_true, y_pred)
  0.5
  >>> metrics.f1_score(y_true, y_pred)  # doctest: +ELLIPSIS
  0.66...
  >>> metrics.fbeta_score(y_true, y_pred, beta=0.5)  # doctest: +ELLIPSIS
  0.83...
  >>> metrics.fbeta_score(y_true, y_pred, beta=1)  # doctest: +ELLIPSIS
  0.66...
  >>> metrics.fbeta_score(y_true, y_pred, beta=2) # doctest: +ELLIPSIS
  0.55...
  >>> metrics.precision_recall_fscore_support(y_true, y_pred, beta=0.5)  # doctest: +ELLIPSIS
  (array([ 0.66...,  1.        ]), array([ 1. ,  0.5]), array([ 0.71...,  0.83...]), array([2, 2]...))


  >>> import numpy as np
  >>> from sklearn.metrics import precision_recall_curve
  >>> from sklearn.metrics import average_precision_score
  >>> y_true = np.array([0, 0, 1, 1])
  >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
  >>> precision, recall, threshold = precision_recall_curve(y_true, y_scores)
  >>> precision  # doctest: +ELLIPSIS
  array([ 0.66...,  0.5       ,  1.        ,  1.        ])
  >>> recall
  array([ 1. ,  0.5,  0.5,  0. ])
  >>> threshold
  array([ 0.35,  0.4 ,  0.8 ])
  >>> average_precision_score(y_true, y_scores)  # doctest: +ELLIPSIS
  0.79...



Multiclass and multilabel classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In multiclass and multilabel classification task, the notions of precision,
recall and F-measures can be applied to each label independently.
There are a few ways to combine results across labels,
specified by the ``average`` argument to the
:func:`average_precision_score` (multilabel only), :func:`f1_score`,
:func:`fbeta_score`, :func:`precision_recall_fscore_support`,
:func:`precision_score` and :func:`recall_score` functions:

* ``"micro"``: calculate metrics globally by counting the total true
  positives, false negatives and false positives. Except in the multi-label
  case this implies that precision, recall and :math:`F` are equal.
* ``"samples"``: calculate metrics for each sample, comparing sets of
  labels assigned to each, and find the mean across all samples.
  This is only meaningful and available in the multilabel case.
* ``"macro"``: calculate metrics for each label, and find their mean.
  This does not take label imbalance into account.
* ``"weighted"``: calculate metrics for each label, and find their average
  weighted by the number of occurrences of the label in the true data.
  This alters ``"macro"`` to account for label imbalance; it may produce an
  F-score that is not between precision and recall.
* ``None``: calculate metrics for each label and do not average them.

To make this more explicit, consider the following notation:

* :math:`y` the set of *predicted* :math:`(sample, label)` pairs
* :math:`\hat{y}` the set of *true* :math:`(sample, label)` pairs
* :math:`L` the set of labels
* :math:`S` the set of samples
* :math:`y_s` the subset of :math:`y` with sample :math:`s`,
  i.e. :math:`y_s := \left\{(s', l) \in y | s' = s\right\}`
* :math:`y_l` the subset of :math:`y` with label :math:`l`
* similarly, :math:`\hat{y}_s` and :math:`\hat{y}_l` are subsets of
  :math:`\hat{y}`
* :math:`P(A, B) := \frac{\left| A \cap B \right|}{\left|A\right|}`
* :math:`R(A, B) := \frac{\left| A \cap B \right|}{\left|B\right|}`
  (Conventions vary on handling :math:`B = \emptyset`; this implementation uses
  :math:`R(A, B):=0`, and similar for :math:`P`.)
* :math:`F_\beta(A, B) := \left(1 + \beta^2\right) \frac{P(A, B) \times R(A, B)}{\beta^2 P(A, B) + R(A, B)}`

Then the metrics are defined as:

+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
|``average``    | Precision                                                                                                        | Recall                                                                                                           | F\_beta                                                                                                              |
+===============+==================================================================================================================+==================================================================================================================+======================================================================================================================+
|``"micro"``    | :math:`P(y, \hat{y})`                                                                                            | :math:`R(y, \hat{y})`                                                                                            | :math:`F_\beta(y, \hat{y})`                                                                                          |
+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
|``"samples"``  | :math:`\frac{1}{\left|S\right|} \sum_{s \in S} P(y_s, \hat{y}_s)`                                                | :math:`\frac{1}{\left|S\right|} \sum_{s \in S} R(y_s, \hat{y}_s)`                                                | :math:`\frac{1}{\left|S\right|} \sum_{s \in S} F_\beta(y_s, \hat{y}_s)`                                              |
+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
|``"macro"``    | :math:`\frac{1}{\left|L\right|} \sum_{l \in L} P(y_l, \hat{y}_l)`                                                | :math:`\frac{1}{\left|L\right|} \sum_{l \in L} R(y_l, \hat{y}_l)`                                                | :math:`\frac{1}{\left|L\right|} \sum_{l \in L} F_\beta(y_l, \hat{y}_l)`                                              |
+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
|``"weighted"`` | :math:`\frac{1}{\sum_{l \in L} \left|\hat{y}_l\right|} \sum_{l \in L} \left|\hat{y}_l\right| P(y_l, \hat{y}_l)`  | :math:`\frac{1}{\sum_{l \in L} \left|\hat{y}_l\right|} \sum_{l \in L} \left|\hat{y}_l\right| R(y_l, \hat{y}_l)`  | :math:`\frac{1}{\sum_{l \in L} \left|\hat{y}_l\right|} \sum_{l \in L} \left|\hat{y}_l\right| F_\beta(y_l, \hat{y}_l)`|
+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+
|``None``       | :math:`\langle P(y_l, \hat{y}_l) | l \in L \rangle`                                                              | :math:`\langle R(y_l, \hat{y}_l) | l \in L \rangle`                                                              | :math:`\langle F_\beta(y_l, \hat{y}_l) | l \in L \rangle`                                                            |
+---------------+------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------+

  >>> from sklearn import metrics
  >>> y_true = [0, 1, 2, 0, 1, 2]
  >>> y_pred = [0, 2, 1, 0, 0, 1]
  >>> metrics.precision_score(y_true, y_pred, average='macro')  # doctest: +ELLIPSIS
  0.22...
  >>> metrics.recall_score(y_true, y_pred, average='micro')
  ... # doctest: +ELLIPSIS
  0.33...
  >>> metrics.f1_score(y_true, y_pred, average='weighted')  # doctest: +ELLIPSIS
  0.26...
  >>> metrics.fbeta_score(y_true, y_pred, average='macro', beta=0.5)  # doctest: +ELLIPSIS
  0.23...
  >>> metrics.precision_recall_fscore_support(y_true, y_pred, beta=0.5, average=None)
  ... # doctest: +ELLIPSIS
  (array([ 0.66...,  0.        ,  0.        ]), array([ 1.,  0.,  0.]), array([ 0.71...,  0.        ,  0.        ]), array([2, 2, 2]...))


Hinge loss
-----------

The :func:`hinge_loss` function computes the average
`hinge loss function <http://en.wikipedia.org/wiki/Hinge_loss>`_. The hinge
loss is used in maximal margin classification as support vector machines.

If the labels are encoded with +1 and -1,  :math:`y`: is the true
value and :math:`w` is the predicted decisions as output by
``decision_function``, then the hinge loss is defined as:

.. math::

  L_\text{Hinge}(y, w) = \max\left\{1 - wy, 0\right\} = \left|1 - wy\right|_+

Here a small example demonstrating the use of the :func:`hinge_loss` function
with a svm classifier::

  >>> from sklearn import svm
  >>> from sklearn.metrics import hinge_loss
  >>> X = [[0], [1]]
  >>> y = [-1, 1]
  >>> est = svm.LinearSVC(random_state=0)
  >>> est.fit(X, y)
  LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
       intercept_scaling=1, loss='l2', max_iter=1000, multi_class='ovr',
       penalty='l2', random_state=0, tol=0.0001, verbose=0)
  >>> pred_decision = est.decision_function([[-2], [3], [0.5]])
  >>> pred_decision  # doctest: +ELLIPSIS
  array([-2.18...,  2.36...,  0.09...])
  >>> hinge_loss([-1, 1, 1], pred_decision)  # doctest: +ELLIPSIS
  0.3...


Log loss
--------

The log loss, also called logistic regression loss or cross-entropy loss,
is a loss function defined on probability estimates.
It is commonly used in (multinomial) logistic regression and neural networks,
as well as some variants of expectation-maximization,
and can be used to evaluate the probability outputs (``predict_proba``)
of a classifier, rather than its discrete predictions.

For binary classification with a true label :math:`y \in \{0,1\}`
and a probability estimate :math:`p = \operatorname{Pr}(y = 1)`,
the log loss per sample is the negative log-likelihood
of the classifier given the true label:

.. math::

    L_{\log}(y, p) = -\log \operatorname{Pr}(y|p) = -(y \log p) + (1 - y) \log (1 - p))

This extends to the multiclass case as follows.
Let the true labels for a set of samples
be encoded as a 1-of-K binary indicator matrix :math:`Y`,
i.e. :math:`y_{i,k} = 1` if sample :math:`i` has label :math:`k`
taken from a set of :math:`K` labels.
Let :math:`P` be a matrix of probability estimates,
with :math:`p_{i,k} = \operatorname{Pr}(t_{i,k} = 1)`.
Then the log loss of the whole set is

.. math::

    L_{\log}(Y, P) = -\log \operatorname{Pr}(Y|P) = - \frac{1}{N} \sum_{i=0}^{N-1} \sum_{k=0}^{K-1} y_{i,k} \log p_{i,k}

To see how this generalizes the binary log loss given above,
note that in the binary case,
:math:`p_{i,0} = 1 - p_{i,1}` and :math:`y_{i,0} = 1 - y_{i,1}`,
so expanding the inner sum over :math:`y_{i,k} \in \{0,1\}`
gives the binary log loss.

The function :func:`log_loss` computes log loss
given a list of ground-truth labels and a probability matrix,
as returned by an estimator's ``predict_proba`` method.

    >>> from sklearn.metrics import log_loss
    >>> y_true = [0, 0, 1, 1]
    >>> y_pred = [[.9, .1], [.8, .2], [.3, .7], [.01, .99]]
    >>> log_loss(y_true, y_pred)    # doctest: +ELLIPSIS
    0.1738...

The first ``[.9, .1]`` in ``y_pred``
denotes 90% probability that the first sample has label 0.
The log loss is non-negative.


Matthews correlation coefficient
---------------------------------

The :func:`matthews_corrcoef` function computes the Matthew's correlation
coefficient (MCC) for binary classes (quoting the `Wikipedia article on the
Matthew's correlation coefficient
<http://en.wikipedia.org/wiki/Matthews_correlation_coefficient>`_):

    "The Matthews correlation coefficient is used in machine learning as a
    measure of the quality of binary (two-class) classifications. It takes
    into account true and false positives and negatives and is generally
    regarded as a balanced measure which can be used even if the classes are
    of very different sizes. The MCC is in essence a correlation coefficient
    value between -1 and +1. A coefficient of +1 represents a perfect
    prediction, 0 an average random prediction and -1 an inverse prediction.
    The statistic is also known as the phi coefficient."

If :math:`tp`, :math:`tn`, :math:`fp` and :math:`fn` are respectively the
number of true positives, true negatives, false positives ans false negatives,
the MCC coefficient is defined as

.. math::

  MCC = \frac{tp \times tn - fp \times fn}{\sqrt{(tp + fp)(tp + fn)(tn + fp)(tn + fn)}}.

Here a small example illustrating the usage of the :func:`matthews_corrcoef`
function:

    >>> from sklearn.metrics import matthews_corrcoef
    >>> y_true = [+1, +1, +1, -1]
    >>> y_pred = [+1, -1, +1, +1]
    >>> matthews_corrcoef(y_true, y_pred)  # doctest: +ELLIPSIS
    -0.33...

.. _roc_metrics:

Receiver operating characteristic (ROC)
---------------------------------------

The function :func:`roc_curve` computes the `receiver operating characteristic
curve, or ROC curve (quoting
Wikipedia) <http://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_:

  "A receiver operating characteristic (ROC), or simply ROC curve, is a
  graphical plot which illustrates the performance of a binary classifier
  system as its discrimination threshold is varied. It is created by plotting
  the fraction of true positives out of the positives (TPR = true positive
  rate) vs. the fraction of false positives out of the negatives (FPR = false
  positive rate), at various threshold settings. TPR is also known as
  sensitivity, and FPR is one minus the specificity or true negative rate."

This function requires the true binary
value and the target scores, which can either be probability estimates of the
positive class, confidence values, or binary decisions.
Here a small example of how to use the :func:`roc_curve` function::

    >>> import numpy as np
    >>> from sklearn.metrics import roc_curve
    >>> y = np.array([1, 1, 2, 2])
    >>> scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = roc_curve(y, scores, pos_label=2)
    >>> fpr
    array([ 0. ,  0.5,  0.5,  1. ])
    >>> tpr
    array([ 0.5,  0.5,  1. ,  1. ])
    >>> thresholds
    array([ 0.8 ,  0.4 ,  0.35,  0.1 ])

The following figure shows an example of such ROC curve.

.. image:: ../auto_examples/model_selection/images/plot_roc_001.png
   :target: ../auto_examples/model_selection/plot_roc.html
   :scale: 75
   :align: center

The :func:`roc_auc_score` function computes the area under the receiver
operating characteristic (ROC) curve, which is also denoted by
AUC or AUROC.  By computing the
area under the roc curve, the curve information is summarized in one number.
For more information see the `Wikipedia article on AUC
<http://en.wikipedia.org/wiki/Receiver_operating_characteristic#Area_under_curve>`_.

  >>> import numpy as np
  >>> from sklearn.metrics import roc_auc_score
  >>> y_true = np.array([0, 0, 1, 1])
  >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
  >>> roc_auc_score(y_true, y_scores)
  0.75

In multi-label classification, the :func:`roc_auc_score` function is
extended by averaging over the labels:

* ``"micro"``: computes the area under the ROC curve globally obtained
  by considering each element of the label indicator matrix as a label.
* ``"samples"``: computes the area under the ROC curve on each sample,
  comparing the set of labels and scores assigned to each, and find the mean
  across all samples.
* ``"macro"``: computes the area under the ROC curve for each label, and find
  their mean.
* ``"weighted"``: computes the area under the ROC curve for each label, and
  find their average weighted by the number of occurrences of the label in the
  true data.
* ``None``: this returns an array of scores with scores with shape (n_classes,)
  instead of an aggregate scalar score.

Compared to metrics such as the subset accuracy, the hamming loss or the
F1 score, ROC AUC doesn't require to optimize a threshold for each label. The
:func:`roc_auc_score` function can also be used in multi-class classification
if predicted outputs have been binarized.


.. image:: ../auto_examples/model_selection/images/plot_roc_002.png
   :target: ../auto_examples/model_selection/plot_roc.html
   :scale: 75
   :align: center

.. topic:: Examples:

  * See :ref:`example_model_selection_plot_roc.py`
    for an example of receiver operating characteristic (ROC) metric to
    evaluate the quality of the output of a classifier.

  * See :ref:`example_model_selection_plot_roc_crossval.py`
    for an example of receiver operating characteristic (ROC) metric to
    evaluate the quality of the output of a classifier using cross-validation.

  * See :ref:`example_applications_plot_species_distribution_modeling.py`
    for an example of receiver operating characteristic (ROC) metric to
    model species distribution.

.. _zero_one_loss:

Zero one loss
--------------

The :func:`zero_one_loss` function computes the sum or the average of the 0-1
classification loss (:math:`L_{0-1}`) over :math:`n_{\text{samples}}`. By
defaults, the function normalizes over the sample. To get the sum of the
:math:`L_{0-1}`, set ``normalize``  to ``False``.

In multilabel classification, the :func:`zero_one_loss` function corresponds
to the subset zero-one loss: the subset of labels must be correctly predict.

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the 0-1 loss :math:`L_{0-1}` is defined as:

.. math::

   L_{0-1}(y_i, \hat{y}_i) = 1(\hat{y}_i \not= y_i)

where :math:`1(x)` is the `indicator function
<http://en.wikipedia.org/wiki/Indicator_function>`_.


  >>> from sklearn.metrics import zero_one_loss
  >>> y_pred = [1, 2, 3, 4]
  >>> y_true = [2, 2, 3, 4]
  >>> zero_one_loss(y_true, y_pred)
  0.25
  >>> zero_one_loss(y_true, y_pred, normalize=False)
  1

In the multilabel case with binary label indicators: ::

  >>> zero_one_loss(np.array([[0, 1], [1, 1]]), np.ones((2, 2)))
  0.5


.. topic:: Example:

  * See :ref:`example_feature_selection_plot_rfe_with_cross_validation.py`
    for an example of the zero one loss usage to perform recursive feature
    elimination with cross-validation.


.. _multilabel_ranking_metrics:

Multilabel ranking metrics
==========================

.. currentmodule:: sklearn.metrics

In multilabel learning, each sample can have any number of ground truth labels
associated with it. The goal is to give high scores and better rank to
the ground truth labels.

Label ranking average precision
-------------------------------

The :func:`label_ranking_average_precision_score` function
implements the label ranking average precision (LRAP). This metric is linked to
the :func:`average_precision_score` function, but is based on the notion of
label ranking instead of precision and recall.

Label ranking average precision (LRAP) is the average over each ground truth
label assigned to each sample, of the ratio of true vs. total labels with lower
score. This metric will yield better score if you are able to give better rank
to the labels associated to each sample. The obtained score is always strictly
greater than 0 and the best value is 1. If there is exactly one relevant
label per sample, label ranking average precision is equivalent to the `mean
reciprocal rank <http://en.wikipedia.org/wiki/Mean_reciprocal_rank>`.

Formally, given a binary indicator matrix of the ground truth labels
:math:`y \in \mathcal{R}^{n_\text{samples} \times n_\text{labels}}` and the
score associated to each label
:math:`\hat{f} \in \mathcal{R}^{n_\text{samples} \times n_\text{labels}}`,
the average precision is defined as

.. math::
  LRAP(y, \hat{f}) = \frac{1}{n_{\text{samples}}}
    \sum_{i=0}^{n_{\text{samples}} - 1} \frac{1}{|y_i|}
    \sum_{j:y_{ij} = 1} \frac{|\mathcal{L}_{ij}|}{\text{rank}_{ij}}


with :math:`\mathcal{L}_{ij} = \left\{k: y_{ik} = 1, \hat{f}_{ik} \geq \hat{f}_{ij} \right\}`,
:math:`\text{rank}_{ij} = \left|\left\{k: \hat{f}_{ik} \geq \hat{f}_{ij} \right\}\right|`
and :math:`|\cdot|` is the l0 norm or the cardinality of the set.

Here a small example of usage of this function::

    >>> import numpy as np
    >>> from sklearn.metrics import label_ranking_average_precision_score
    >>> y_true = np.array([[1, 0, 0], [0, 0, 1]])
    >>> y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
    >>> label_ranking_average_precision_score(y_true, y_score) # doctest: +ELLIPSIS
    0.416...


.. _regression_metrics:

Regression metrics
===================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
functions to measure regression performance. Some of those have been enhanced
to handle the multioutput case: :func:`mean_absolute_error`,
:func:`mean_absolute_error` and :func:`r2_score`.


Explained variance score
-------------------------

The :func:`explained_variance_score` computes the `explained variance
regression score <http://en.wikipedia.org/wiki/Explained_variation>`_.

If :math:`\hat{y}` is the estimated target output
and :math:`y` is the corresponding (correct) target output, then the explained
variance is  estimated  as follow:

.. math::

  \texttt{explained\_{}variance}(y, \hat{y}) = 1 - \frac{Var\{ y - \hat{y}\}}{Var\{y\}}

The best possible score is 1.0, lower values are worse.

Here a small example of usage of the :func:`explained_variance_score`
function::

    >>> from sklearn.metrics import explained_variance_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> explained_variance_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.957...

Mean absolute error
-------------------

The :func:`mean_absolute_error` function computes the `mean absolute
error <http://en.wikipedia.org/wiki/Mean_absolute_error>`_, which is a risk
function corresponding to the expected value of the absolute error loss or
:math:`l1`-norm loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean absolute error
(MAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MAE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}-1} \left| y_i - \hat{y}_i \right|.

Here a small example of usage of the :func:`mean_absolute_error` function::

  >>> from sklearn.metrics import mean_absolute_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> mean_absolute_error(y_true, y_pred)
  0.5
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> mean_absolute_error(y_true, y_pred)
  0.75



Mean squared error
-------------------

The :func:`mean_squared_error` function computes the `mean square
error <http://en.wikipedia.org/wiki/Mean_squared_error>`_, which is a risk
function corresponding to the expected value of the squared error loss or
quadratic loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean squared error
(MSE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MSE}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (y_i - \hat{y}_i)^2.

Here a small example of usage of the :func:`mean_squared_error`
function::

  >>> from sklearn.metrics import mean_squared_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> mean_squared_error(y_true, y_pred)
  0.375
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> mean_squared_error(y_true, y_pred)  # doctest: +ELLIPSIS
  0.7083...

.. topic:: Examples:

  * See :ref:`example_ensemble_plot_gradient_boosting_regression.py`
    for an example of mean squared error usage to
    evaluate gradient boosting regression.

R² score, the coefficient of determination
-------------------------------------------

The :func:`r2_score` function computes R², the `coefficient of
determination <http://en.wikipedia.org/wiki/Coefficient_of_determination>`_.
It provides a measure of how well future samples are likely to
be predicted by the model.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the score R² estimated
over :math:`n_{\text{samples}}` is defined as

.. math::

  R^2(y, \hat{y}) = 1 - \frac{\sum_{i=0}^{n_{\text{samples}} - 1} (y_i - \hat{y}_i)^2}{\sum_{i=0}^{n_\text{samples} - 1} (y_i - \bar{y})^2}

where :math:`\bar{y} =  \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}} - 1} y_i`.

Here a small example of usage of the :func:`r2_score` function::

  >>> from sklearn.metrics import r2_score
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
  0.948...
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
  0.938...


.. topic:: Example:

  * See :ref:`example_linear_model_plot_lasso_and_elasticnet.py`
    for an example of R² score usage to
    evaluate Lasso and Elastic Net on sparse signals.

.. _clustering_metrics:

Clustering metrics
======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
functions. For more information see the :ref:`clustering_evaluation`
section for instance clustering, and :ref:`biclustering_evaluation` for
biclustering.


.. _dummy_estimators:

Dummy estimators
=================

.. currentmodule:: sklearn.dummy

When doing supervised learning, a simple sanity check consists in comparing
one's estimator against simple rules of thumb. :class:`DummyClassifier`
implements three such simple strategies for classification:

- ``stratified`` generates randomly predictions by respecting the training
  set's class distribution,
- ``most_frequent`` always predicts the most frequent label in the training set,
- ``uniform`` generates predictions uniformly at random.
- ``constant`` always predicts a constant label that is provided by the user.
   A major motivation of this method is F1-scoring when the positive class
   is in the minority.

Note that with all these strategies, the ``predict`` method completely ignores
the input data!

To illustrate :class:`DummyClassifier`, first let's create an imbalanced
dataset::

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.cross_validation import train_test_split
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> y[y != 1] = -1
  >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

Next, let's compare the accuracy of ``SVC`` and ``most_frequent``::

  >>> from sklearn.dummy import DummyClassifier
  >>> from sklearn.svm import SVC
  >>> clf = SVC(kernel='linear', C=1).fit(X_train, y_train)
  >>> clf.score(X_test, y_test) # doctest: +ELLIPSIS
  0.63...
  >>> clf = DummyClassifier(strategy='most_frequent',random_state=0)
  >>> clf.fit(X_train, y_train)
  DummyClassifier(constant=None, random_state=0, strategy='most_frequent')
  >>> clf.score(X_test, y_test)  # doctest: +ELLIPSIS
  0.57...

We see that ``SVC`` doesn't do much better than a dummy classifier. Now, let's
change the kernel::

  >>> clf = SVC(kernel='rbf', C=1).fit(X_train, y_train)
  >>> clf.score(X_test, y_test)  # doctest: +ELLIPSIS
  0.97...

We see that the accuracy was boosted to almost 100%. For a better estimate
of the accuracy, it is recommended to use a cross validation strategy, if it
is not too CPU costly. For more information see the :ref:`cross_validation`
section. Moreover if you want to optimize over the parameter space, it is
highly recommended to use an appropriate methodology see the :ref:`grid_search`
section.

More generally, when the accuracy of a classifier is too close to random
classification, it probably means that something went wrong: features are not
helpful, a hyper parameter is not correctly tuned, the classifier is suffering
from class imbalance, etc...

:class:`DummyRegressor` also implements three simple rules of thumb for regression:

- ``mean`` always predicts the mean of the training targets.
- ``median`` always predicts the median of the training targests.
- ``constant`` always predicts a constant value that is provided by the user.

In all these strategies, the ``predict`` method completely ignores
the input data.
