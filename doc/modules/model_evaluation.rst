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

.. _scoring_parameter:

The ``scoring`` parameter: defining model evaluation rules
==========================================================

Model selection and evaluation using tools, such as
:class:`model_selection.GridSearchCV` and
:func:`model_selection.cross_val_score`, take a ``scoring`` parameter that
controls what metric they apply to the estimators evaluated.

Common cases: predefined values
-------------------------------

For the most common use cases, you can designate a scorer object with the
``scoring`` parameter; the table below shows all possible values.
All scorer objects follow the convention that higher return values are better
than lower return values.  Thus the returns from mean_absolute_error
and mean_squared_error, which measure the distance between the model
and the data, are negated.


========================     =======================================     ==================================
Scoring                      Function                                    Comment
========================     =======================================     ==================================
**Classification**
'accuracy'                   :func:`metrics.accuracy_score`
'average_precision'          :func:`metrics.average_precision_score`
'f1'                         :func:`metrics.f1_score`                    for binary targets
'f1_micro'                   :func:`metrics.f1_score`                    micro-averaged
'f1_macro'                   :func:`metrics.f1_score`                    macro-averaged
'f1_weighted'                :func:`metrics.f1_score`                    weighted average
'f1_samples'                 :func:`metrics.f1_score`                    by multilabel sample
'log_loss'                   :func:`metrics.log_loss`                    requires ``predict_proba`` support
'precision' etc.             :func:`metrics.precision_score`             suffixes apply as with 'f1'
'recall' etc.                :func:`metrics.recall_score`                suffixes apply as with 'f1'
'roc_auc'                    :func:`metrics.roc_auc_score`

**Clustering**
'adjusted_rand_score'        :func:`metrics.adjusted_rand_score`

**Regression**
'mean_absolute_error'        :func:`metrics.mean_absolute_error`
'mean_squared_error'         :func:`metrics.mean_squared_error`
'median_absolute_error'      :func:`metrics.median_absolute_error`
'r2'                         :func:`metrics.r2_score`
========================     =======================================     ==================================

Usage examples:

    >>> from sklearn import svm, datasets
    >>> from sklearn.model_selection import cross_val_score
    >>> iris = datasets.load_iris()
    >>> X, y = iris.data, iris.target
    >>> clf = svm.SVC(probability=True, random_state=0)
    >>> cross_val_score(clf, X, y, scoring='log_loss') # doctest: +ELLIPSIS
    array([-0.07..., -0.16..., -0.06...])
    >>> model = svm.SVC()
    >>> cross_val_score(model, X, y, scoring='wrong_choice')
    Traceback (most recent call last):
    ValueError: 'wrong_choice' is not a valid scoring value. Valid options are ['accuracy', 'adjusted_rand_score', 'average_precision', 'f1', 'f1_macro', 'f1_micro', 'f1_samples', 'f1_weighted', 'log_loss', 'mean_absolute_error', 'mean_squared_error', 'median_absolute_error', 'precision', 'precision_macro', 'precision_micro', 'precision_samples', 'precision_weighted', 'r2', 'recall', 'recall_macro', 'recall_micro', 'recall_samples', 'recall_weighted', 'roc_auc']

.. note::

    The values listed by the ValueError exception correspond to the functions measuring
    prediction accuracy described in the following sections.
    The scorer objects for those functions are stored in the dictionary
    ``sklearn.metrics.SCORERS``.

.. currentmodule:: sklearn.metrics

.. _scoring:

Defining your scoring strategy from metric functions
-----------------------------------------------------

The module :mod:`sklearn.metric` also exposes a set of simple functions
measuring a prediction error given ground truth and prediction:

- functions ending with ``_score`` return a value to
  maximize, the higher the better.

- functions ending with ``_error`` or ``_loss`` return a
  value to minimize, the lower the better.  When converting
  into a scorer object using :func:`make_scorer`, set
  the ``greater_is_better`` parameter to False (True by default; see the
  parameter description below).

Metrics available for various machine learning tasks are detailed in sections
below.

Many metrics are not given names to be used as ``scoring`` values,
sometimes because they require additional parameters, such as
:func:`fbeta_score`. In such cases, you need to generate an appropriate
scoring object.  The simplest way to generate a callable object for scoring
is by using :func:`make_scorer`. That function converts metrics
into callables that can be used for model evaluation.

One typical use case is to wrap an existing metric function from the library
with non-default values for its parameters, such as the ``beta`` parameter for
the :func:`fbeta_score` function::

    >>> from sklearn.metrics import fbeta_score, make_scorer
    >>> ftwo_scorer = make_scorer(fbeta_score, beta=2)
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]}, scoring=ftwo_scorer)

The second use case is to build a completely custom scorer object
from a simple python function using :func:`make_scorer`, which can
take several parameters:

* the python function you want to use (``my_custom_loss_func``
  in the example below)

* whether the python function returns a score (``greater_is_better=True``,
  the default) or a loss (``greater_is_better=False``).  If a loss, the output
  of the python function is negated by the scorer object, conforming to
  the cross validation convention that scorers return higher values for better models.

* for classification metrics only: whether the python function you provided requires continuous decision
  certainties (``needs_threshold=True``).  The default value is
  False.

* any additional parameters, such as ``beta`` or ``labels`` in :func:`f1_score`.

Here is an example of building custom scorers, and of using the
``greater_is_better`` parameter::

    >>> import numpy as np
    >>> def my_custom_loss_func(ground_truth, predictions):
    ...     diff = np.abs(ground_truth - predictions).max()
    ...     return np.log(1 + diff)
    ...
    >>> # loss_func will negate the return value of my_custom_loss_func,
    >>> #  which will be np.log(2), 0.693, given the values for ground_truth
    >>> #  and predictions defined below.
    >>> loss  = make_scorer(my_custom_loss_func, greater_is_better=False)
    >>> score = make_scorer(my_custom_loss_func, greater_is_better=True)
    >>> ground_truth = [[1, 1]]
    >>> predictions  = [0, 1]
    >>> from sklearn.dummy import DummyClassifier
    >>> clf = DummyClassifier(strategy='most_frequent', random_state=0)
    >>> clf = clf.fit(ground_truth, predictions)
    >>> loss(clf,ground_truth, predictions) # doctest: +ELLIPSIS
    -0.69...
    >>> score(clf,ground_truth, predictions) # doctest: +ELLIPSIS
    0.69...


.. _diy_scoring:

Implementing your own scoring object
------------------------------------
You can generate even more flexible model scorers by constructing your own
scoring object from scratch, without using the :func:`make_scorer` factory.
For a callable to be a scorer, it needs to meet the protocol specified by
the following two rules:

- It can be called with parameters ``(estimator, X, y)``, where ``estimator``
  is the model that should be evaluated, ``X`` is validation data, and ``y`` is
  the ground truth target for ``X`` (in the supervised case) or ``None`` (in the
  unsupervised case).

- It returns a floating point number that quantifies the
  ``estimator`` prediction quality on ``X``, with reference to ``y``.
  Again, by convention higher numbers are better, so if your scorer
  returns loss, that value should be negated.


.. _classification_metrics:

Classification metrics
=======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` module implements several loss, score, and utility
functions to measure classification performance.
Some metrics might require probability estimates of the positive class,
confidence values, or binary decisions values.
Most implementations allow each sample to provide a weighted contribution
to the overall score, through the ``sample_weight`` parameter.

Some of these are restricted to the binary classification case:

.. autosummary::
   :template: function.rst

   matthews_corrcoef
   precision_recall_curve
   roc_curve


Others also work in the multiclass case:

.. autosummary::
   :template: function.rst

   confusion_matrix
   hinge_loss


Some also work in the multilabel case:

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

And some work with binary and multilabel (but not multiclass) problems:

.. autosummary::
   :template: function.rst

   average_precision_score
   roc_auc_score


In the following sub-sections, we will describe each of those functions,
preceded by some notes on common API and metric definition.

From binary to multiclass and multilabel
----------------------------------------

Some metrics are essentially defined for binary classification tasks (e.g.
:func:`f1_score`, :func:`roc_auc_score`). In these cases, by default
only the positive label is evaluated, assuming by default that the positive
class is labelled ``1`` (though this may be configurable through the
``pos_label`` parameter).

.. _average:

In extending a binary metric to multiclass or multilabel problems, the data
is treated as a collection of binary problems, one for each class.
There are then a number of ways to average binary metric calculations across
the set of classes, each of which may be useful in some scenario.
Where available, you should select among these using the ``average`` parameter.

* ``"macro"`` simply calculates the mean of the binary metrics,
  giving equal weight to each class.  In problems where infrequent classes
  are nonetheless important, macro-averaging may be a means of highlighting
  their performance. On the other hand, the assumption that all classes are
  equally important is often untrue, such that macro-averaging will
  over-emphasize the typically low performance on an infrequent class.
* ``"weighted"`` accounts for class imbalance by computing the average of
  binary metrics in which each class's score is weighted by its presence in the
  true data sample.
* ``"micro"`` gives each sample-class pair an equal contribution to the overall
  metric (except as a result of sample-weight). Rather than summing the
  metric per class, this sums the dividends and divisors that make up the the
  per-class metrics to calculate an overall quotient.
  Micro-averaging may be preferred in multilabel settings, including
  multiclass classification where a majority class is to be ignored.
* ``"samples"`` applies only to multilabel problems. It does not calculate a
  per-class measure, instead calculating the metric over the true and predicted
  classes for each sample in the evaluation data, and returning their
  (``sample_weight``-weighted) average.
* Selecting ``average=None`` will return an array with the score for each
  class.

While multiclass data is provided to the metric, like binary targets, as an
array of class labels, multilabel data is specified as an indicator matrix,
in which cell ``[i, j]`` has value 1 if sample ``i`` has label ``j`` and value
0 otherwise.

.. _accuracy_score:

Accuracy score
--------------

The :func:`accuracy_score` function computes the
`accuracy <https://en.wikipedia.org/wiki/Accuracy_and_precision>`_, either the fraction
(default) or the count (normalize=False) of correct predictions.


In multilabel classification, the function returns the subset accuracy. If
the entire set of predicted labels for a sample strictly match with the true
set of labels, then the subset accuracy is 1.0; otherwise it is 0.0.

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the fraction of correct predictions over :math:`n_\text{samples}` is
defined as

.. math::

   \texttt{accuracy}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples}-1} 1(\hat{y}_i = y_i)

where :math:`1(x)` is the `indicator function
<https://en.wikipedia.org/wiki/Indicator_function>`_.

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

.. _cohen_kappa:

Cohen's kappa
-------------

The function :func:`cohen_kappa_score` computes Cohen's kappa statistic.
This measure is intended to compare labelings by different human annotators,
not a classifier versus a ground truth.

The kappa score (see docstring) is a number between -1 and 1.
Scores above .8 are generally considered good agreement;
zero or lower means no agreement (practically random labels).

Kappa scores can be computed for binary or multiclass problems,
but not for multilabel problems (except by manually computing a per-label score)
and not for more than two annotators.

.. _confusion_matrix:

Confusion matrix
----------------

The :func:`confusion_matrix` function evaluates
classification accuracy by computing the `confusion matrix
<https://en.wikipedia.org/wiki/Confusion_matrix>`_.

By definition, entry :math:`i, j` in a confusion matrix is
the number of observations actually in group :math:`i`, but
predicted to be in group :math:`j`. Here is an example::

  >>> from sklearn.metrics import confusion_matrix
  >>> y_true = [2, 0, 2, 2, 0, 1]
  >>> y_pred = [0, 0, 2, 2, 0, 2]
  >>> confusion_matrix(y_true, y_pred)
  array([[2, 0, 0],
         [0, 0, 1],
         [1, 0, 2]])

Here is a visual representation of such a confusion matrix (this figure comes
from the :ref:`example_model_selection_plot_confusion_matrix.py` example):

.. image:: ../auto_examples/model_selection/images/plot_confusion_matrix_001.png
   :target: ../auto_examples/model_selection/plot_confusion_matrix.html
   :scale: 75
   :align: center

.. topic:: Example:

  * See :ref:`example_model_selection_plot_confusion_matrix.py`
    for an example of using a confusion matrix to evaluate classifier output
    quality.

  * See :ref:`example_classification_plot_digits_classification.py`
    for an example of using a confusion matrix to classify
    hand-written digits.

  * See :ref:`example_text_document_classification_20newsgroups.py`
    for an example of using a confusion matrix to classify text
    documents.

.. _classification_report:

Classification report
----------------------

The :func:`classification_report` function builds a text report showing the
main classification metrics. Here is a small example with custom ``target_names``
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
    for an example of classification report usage for
    hand-written digits.

  * See :ref:`example_text_document_classification_20newsgroups.py`
    for an example of classification report usage for text
    documents.

  * See :ref:`example_model_selection_grid_search_digits.py`
    for an example of classification report usage for
    grid search with nested cross-validation.

.. _hamming_loss:

Hamming loss
-------------

The :func:`hamming_loss` computes the average Hamming loss or `Hamming
distance <https://en.wikipedia.org/wiki/Hamming_distance>`_ between two sets
of samples.

If :math:`\hat{y}_j` is the predicted value for the :math:`j`-th label of
a given sample, :math:`y_j` is the corresponding true value, and
:math:`n_\text{labels}` is the number of classes or labels, then the
Hamming loss :math:`L_{Hamming}` between two samples is defined as:

.. math::

   L_{Hamming}(y, \hat{y}) = \frac{1}{n_\text{labels}} \sum_{j=0}^{n_\text{labels} - 1} 1(\hat{y}_j \not= y_j)

where :math:`1(x)` is the `indicator function
<https://en.wikipedia.org/wiki/Indicator_function>`_. ::

  >>> from sklearn.metrics import hamming_loss
  >>> y_pred = [1, 2, 3, 4]
  >>> y_true = [2, 2, 3, 4]
  >>> hamming_loss(y_true, y_pred)
  0.25

In the multilabel case with binary label indicators: ::

  >>> hamming_loss(np.array([[0, 1], [1, 1]]), np.zeros((2, 2)))
  0.75

.. note::

    In multiclass classification, the Hamming loss corresponds to the Hamming
    distance between ``y_true`` and ``y_pred`` which is similar to the
    :ref:`zero_one_loss` function.  However, while zero-one loss penalizes
    prediction sets that do not strictly match true sets, the Hamming loss
    penalizes individual labels.  Thus the Hamming loss, upper bounded by the zero-one
    loss, is always between zero and one, inclusive; and predicting a proper subset
    or superset of the true labels will give a Hamming loss between
    zero and one, exclusive.

.. _jaccard_similarity_score:

Jaccard similarity coefficient score
-------------------------------------

The :func:`jaccard_similarity_score` function computes the average (default)
or sum of `Jaccard similarity coefficients
<https://en.wikipedia.org/wiki/Jaccard_index>`_, also called the Jaccard index,
between pairs of label sets.

The Jaccard similarity coefficient of the :math:`i`-th samples,
with a ground truth label set :math:`y_i` and predicted label set
:math:`\hat{y}_i`, is defined as

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

Intuitively, `precision
<https://en.wikipedia.org/wiki/Precision_and_recall#Precision>`_ is the ability
of the classifier not to label as positive a sample that is negative, and
`recall <https://en.wikipedia.org/wiki/Precision_and_recall#Recall>`_ is the
ability of the classifier to find all the positive samples.

The  `F-measure <https://en.wikipedia.org/wiki/F1_score>`_
(:math:`F_\beta` and :math:`F_1` measures) can be interpreted as a weighted
harmonic mean of the precision and recall. A
:math:`F_\beta` measure reaches its best value at 1 and its worst score at 0.
With :math:`\beta = 1`,  :math:`F_\beta` and
:math:`F_1`  are equivalent, and the recall and the precision are equally important.

The :func:`precision_recall_curve` computes a precision-recall curve
from the ground truth label and a score given by the classifier
by varying a decision threshold.

The :func:`average_precision_score` function computes the average precision
(AP) from prediction scores. This score corresponds to the area under the
precision-recall curve. The value is between 0 and 1 and higher is better.
With random predictions, the AP is the fraction of positive samples.

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
    for an example of :func:`f1_score` usage to classify  text
    documents.

  * See :ref:`example_model_selection_grid_search_digits.py`
    for an example of :func:`precision_score` and :func:`recall_score` usage
    to estimate parameters using grid search with nested cross-validation.

  * See :ref:`example_model_selection_plot_precision_recall.py`
    for an example of :func:`precision_recall_curve` usage to evaluate
    classifier output quality.

  * See :ref:`example_linear_model_plot_sparse_recovery.py`
    for an example of :func:`precision_recall_curve` usage to select
    features for sparse linear models.

Binary classification
^^^^^^^^^^^^^^^^^^^^^

In a binary classification task, the terms ''positive'' and ''negative'' refer
to the classifier's prediction, and the terms ''true'' and ''false'' refer to
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

Here are some small examples in binary classification::

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
recall, and F-measures can be applied to each label independently.
There are a few ways to combine results across labels,
specified by the ``average`` argument to the
:func:`average_precision_score` (multilabel only), :func:`f1_score`,
:func:`fbeta_score`, :func:`precision_recall_fscore_support`,
:func:`precision_score` and :func:`recall_score` functions, as described
:ref:`above <average>`. Note that for "micro"-averaging in a multiclass setting
with all labels included will produce equal precision, recall and :math:`F`,
while "weighted" averaging may produce an F-score that is not between
precision and recall.

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

For multiclass classification with a "negative class", it is possible to exclude some labels:

  >>> metrics.recall_score(y_true, y_pred, labels=[1, 2], average='micro')
  ... # excluding 0, no labels were correctly recalled
  0.0

Similarly, labels not present in the data sample may be accounted for in macro-averaging.

  >>> metrics.precision_score(y_true, y_pred, labels=[0, 1, 2, 3], average='macro')
  ... # doctest: +ELLIPSIS
  0.166...

.. _hinge_loss:

Hinge loss
----------

The :func:`hinge_loss` function computes the average distance between
the model and the data using
`hinge loss <https://en.wikipedia.org/wiki/Hinge_loss>`_, a one-sided metric
that considers only prediction errors. (Hinge
loss is used in maximal margin classifiers such as support vector machines.)

If the labels are encoded with +1 and -1,  :math:`y`: is the true
value, and :math:`w` is the predicted decisions as output by
``decision_function``, then the hinge loss is defined as:

.. math::

  L_\text{Hinge}(y, w) = \max\left\{1 - wy, 0\right\} = \left|1 - wy\right|_+

If there are more than two labels, :func:`hinge_loss` uses a multiclass variant
due to Crammer & Singer.
`Here <http://jmlr.csail.mit.edu/papers/volume2/crammer01a/crammer01a.pdf>`_ is
the paper describing it.

If :math:`y_w` is the predicted decision for true label and :math:`y_t` is the
maximum of the predicted decisions for all other labels, where predicted
decisions are output by decision function, then multiclass hinge loss is defined
by:

.. math::

  L_\text{Hinge}(y_w, y_t) = \max\left\{1 + y_t - y_w, 0\right\}

Here a small example demonstrating the use of the :func:`hinge_loss` function
with a svm classifier in a binary class problem::

  >>> from sklearn import svm
  >>> from sklearn.metrics import hinge_loss
  >>> X = [[0], [1]]
  >>> y = [-1, 1]
  >>> est = svm.LinearSVC(random_state=0)
  >>> est.fit(X, y)
  LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
       intercept_scaling=1, loss='squared_hinge', max_iter=1000,
       multi_class='ovr', penalty='l2', random_state=0, tol=0.0001,
       verbose=0)
  >>> pred_decision = est.decision_function([[-2], [3], [0.5]])
  >>> pred_decision  # doctest: +ELLIPSIS
  array([-2.18...,  2.36...,  0.09...])
  >>> hinge_loss([-1, 1, 1], pred_decision)  # doctest: +ELLIPSIS
  0.3...

Here is an example demonstrating the use of the :func:`hinge_loss` function
with a svm classifier in a multiclass problem::

  >>> X = np.array([[0], [1], [2], [3]])
  >>> Y = np.array([0, 1, 2, 3])
  >>> labels = np.array([0, 1, 2, 3])
  >>> est = svm.LinearSVC()
  >>> est.fit(X, Y)
  LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
       intercept_scaling=1, loss='squared_hinge', max_iter=1000,
       multi_class='ovr', penalty='l2', random_state=None, tol=0.0001,
       verbose=0)
  >>> pred_decision = est.decision_function([[-1], [2], [3]])
  >>> y_true = [0, 2, 3]
  >>> hinge_loss(y_true, pred_decision, labels)  #doctest: +ELLIPSIS
  0.56...

.. _log_loss:

Log loss
--------

Log loss, also called logistic regression loss or
cross-entropy loss, is defined on probability estimates.  It is
commonly used in (multinomial) logistic regression and neural networks, as well
as in some variants of expectation-maximization, and can be used to evaluate the
probability outputs (``predict_proba``) of a classifier instead of its
discrete predictions.

For binary classification with a true label :math:`y \in \{0,1\}`
and a probability estimate :math:`p = \operatorname{Pr}(y = 1)`,
the log loss per sample is the negative log-likelihood
of the classifier given the true label:

.. math::

    L_{\log}(y, p) = -\log \operatorname{Pr}(y|p) = -(y \log (p) + (1 - y) \log (1 - p))

This extends to the multiclass case as follows.
Let the true labels for a set of samples
be encoded as a 1-of-K binary indicator matrix :math:`Y`,
i.e., :math:`y_{i,k} = 1` if sample :math:`i` has label :math:`k`
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

The :func:`log_loss` function computes log loss given a list of ground-truth
labels and a probability matrix, as returned by an estimator's ``predict_proba``
method.

    >>> from sklearn.metrics import log_loss
    >>> y_true = [0, 0, 1, 1]
    >>> y_pred = [[.9, .1], [.8, .2], [.3, .7], [.01, .99]]
    >>> log_loss(y_true, y_pred)    # doctest: +ELLIPSIS
    0.1738...

The first ``[.9, .1]`` in ``y_pred`` denotes 90% probability that the first
sample has label 0.  The log loss is non-negative.

.. _matthews_corrcoef:

Matthews correlation coefficient
---------------------------------

The :func:`matthews_corrcoef` function computes the
`Matthew's correlation coefficient (MCC) <https://en.wikipedia.org/wiki/Matthews_correlation_coefficient>`_
for binary classes.  Quoting Wikipedia:


    "The Matthews correlation coefficient is used in machine learning as a
    measure of the quality of binary (two-class) classifications. It takes
    into account true and false positives and negatives and is generally
    regarded as a balanced measure which can be used even if the classes are
    of very different sizes. The MCC is in essence a correlation coefficient
    value between -1 and +1. A coefficient of +1 represents a perfect
    prediction, 0 an average random prediction and -1 an inverse prediction.
    The statistic is also known as the phi coefficient."

If :math:`tp`, :math:`tn`, :math:`fp` and :math:`fn` are respectively the
number of true positives, true negatives, false positives and false negatives,
the MCC coefficient is defined as

.. math::

  MCC = \frac{tp \times tn - fp \times fn}{\sqrt{(tp + fp)(tp + fn)(tn + fp)(tn + fn)}}.

Here is a small example illustrating the usage of the :func:`matthews_corrcoef`
function:

    >>> from sklearn.metrics import matthews_corrcoef
    >>> y_true = [+1, +1, +1, -1]
    >>> y_pred = [+1, -1, +1, +1]
    >>> matthews_corrcoef(y_true, y_pred)  # doctest: +ELLIPSIS
    -0.33...

.. _roc_metrics:

Receiver operating characteristic (ROC)
---------------------------------------

The function :func:`roc_curve` computes the
`receiver operating characteristic curve, or ROC curve <https://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_.
Quoting Wikipedia :

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
Here is a small example of how to use the :func:`roc_curve` function::

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

This figure shows an example of such an ROC curve:

.. image:: ../auto_examples/model_selection/images/plot_roc_001.png
   :target: ../auto_examples/model_selection/plot_roc.html
   :scale: 75
   :align: center

The :func:`roc_auc_score` function computes the area under the receiver
operating characteristic (ROC) curve, which is also denoted by
AUC or AUROC.  By computing the
area under the roc curve, the curve information is summarized in one number.
For more information see the `Wikipedia article on AUC
<https://en.wikipedia.org/wiki/Receiver_operating_characteristic#Area_under_the_curve>`_.

  >>> import numpy as np
  >>> from sklearn.metrics import roc_auc_score
  >>> y_true = np.array([0, 0, 1, 1])
  >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
  >>> roc_auc_score(y_true, y_scores)
  0.75

In multi-label classification, the :func:`roc_auc_score` function is
extended by averaging over the labels as :ref:`above <average>`.

Compared to metrics such as the subset accuracy, the Hamming loss, or the
F1 score, ROC doesn't require optimizing a threshold for each label. The
:func:`roc_auc_score` function can also be used in multi-class classification,
if the predicted outputs have been binarized.


.. image:: ../auto_examples/model_selection/images/plot_roc_002.png
   :target: ../auto_examples/model_selection/plot_roc.html
   :scale: 75
   :align: center

.. topic:: Examples:

  * See :ref:`example_model_selection_plot_roc.py`
    for an example of using ROC to
    evaluate the quality of the output of a classifier.

  * See :ref:`example_model_selection_plot_roc_crossval.py`
    for an example of using ROC to
    evaluate classifier output quality, using cross-validation.

  * See :ref:`example_applications_plot_species_distribution_modeling.py`
    for an example of using ROC to
    model species distribution.

.. _zero_one_loss:

Zero one loss
--------------

The :func:`zero_one_loss` function computes the sum or the average of the 0-1
classification loss (:math:`L_{0-1}`) over :math:`n_{\text{samples}}`. By
default, the function normalizes over the sample. To get the sum of the
:math:`L_{0-1}`, set ``normalize`` to ``False``.

In multilabel classification, the :func:`zero_one_loss` scores a subset as
one if its labels strictly match the predictions, and as a zero if there
are any errors.  By default, the function returns the percentage of imperfectly
predicted subsets.  To get the count of such subsets instead, set
``normalize`` to ``False``

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the 0-1 loss :math:`L_{0-1}` is defined as:

.. math::

   L_{0-1}(y_i, \hat{y}_i) = 1(\hat{y}_i \not= y_i)

where :math:`1(x)` is the `indicator function
<https://en.wikipedia.org/wiki/Indicator_function>`_.


  >>> from sklearn.metrics import zero_one_loss
  >>> y_pred = [1, 2, 3, 4]
  >>> y_true = [2, 2, 3, 4]
  >>> zero_one_loss(y_true, y_pred)
  0.25
  >>> zero_one_loss(y_true, y_pred, normalize=False)
  1

In the multilabel case with binary label indicators, where the first label
set [0,1] has an error: ::

  >>> zero_one_loss(np.array([[0, 1], [1, 1]]), np.ones((2, 2)))
  0.5

  >>> zero_one_loss(np.array([[0, 1], [1, 1]]), np.ones((2, 2)),  normalize=False)
  1

.. topic:: Example:

  * See :ref:`example_feature_selection_plot_rfe_with_cross_validation.py`
    for an example of zero one loss usage to perform recursive feature
    elimination with cross-validation.


.. _multilabel_ranking_metrics:

Multilabel ranking metrics
==========================

.. currentmodule:: sklearn.metrics

In multilabel learning, each sample can have any number of ground truth labels
associated with it. The goal is to give high scores and better rank to
the ground truth labels.

.. _coverage_error:

Coverage error
--------------

The :func:`coverage_error` function computes the average number of labels that
have to be included in the final prediction such that all true labels
are predicted. This is useful if you want to know how many top-scored-labels
you have to predict in average without missing any true one. The best value
of this metrics is thus the average number of true labels.

Formally, given a binary indicator matrix of the ground truth labels
:math:`y \in \left\{0, 1\right\}^{n_\text{samples} \times n_\text{labels}}` and the
score associated with each label
:math:`\hat{f} \in \mathbb{R}^{n_\text{samples} \times n_\text{labels}}`,
the coverage is defined as

.. math::
  coverage(y, \hat{f}) = \frac{1}{n_{\text{samples}}}
    \sum_{i=0}^{n_{\text{samples}} - 1} \max_{j:y_{ij} = 1} \text{rank}_{ij}

with :math:`\text{rank}_{ij} = \left|\left\{k: \hat{f}_{ik} \geq \hat{f}_{ij} \right\}\right|`.
Given the rank definition, ties in ``y_scores`` are broken by giving the
maximal rank that would have been assigned to all tied values.

Here is a small example of usage of this function::

    >>> import numpy as np
    >>> from sklearn.metrics import coverage_error
    >>> y_true = np.array([[1, 0, 0], [0, 0, 1]])
    >>> y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
    >>> coverage_error(y_true, y_score)
    2.5

.. _label_ranking_average_precision:

Label ranking average precision
-------------------------------

The :func:`label_ranking_average_precision_score` function
implements label ranking average precision (LRAP). This metric is linked to
the :func:`average_precision_score` function, but is based on the notion of
label ranking instead of precision and recall.

Label ranking average precision (LRAP) is the average over each ground truth
label assigned to each sample, of the ratio of true vs. total labels with lower
score. This metric will yield better scores if you are able to give better rank
to the labels associated with each sample. The obtained score is always strictly
greater than 0, and the best value is 1. If there is exactly one relevant
label per sample, label ranking average precision is equivalent to the `mean
reciprocal rank <https://en.wikipedia.org/wiki/Mean_reciprocal_rank>`_.

Formally, given a binary indicator matrix of the ground truth labels
:math:`y \in \mathcal{R}^{n_\text{samples} \times n_\text{labels}}` and the
score associated with each label
:math:`\hat{f} \in \mathcal{R}^{n_\text{samples} \times n_\text{labels}}`,
the average precision is defined as

.. math::
  LRAP(y, \hat{f}) = \frac{1}{n_{\text{samples}}}
    \sum_{i=0}^{n_{\text{samples}} - 1} \frac{1}{|y_i|}
    \sum_{j:y_{ij} = 1} \frac{|\mathcal{L}_{ij}|}{\text{rank}_{ij}}


with :math:`\mathcal{L}_{ij} = \left\{k: y_{ik} = 1, \hat{f}_{ik} \geq \hat{f}_{ij} \right\}`,
:math:`\text{rank}_{ij} = \left|\left\{k: \hat{f}_{ik} \geq \hat{f}_{ij} \right\}\right|`
and :math:`|\cdot|` is the l0 norm or the cardinality of the set.

Here is a small example of usage of this function::

    >>> import numpy as np
    >>> from sklearn.metrics import label_ranking_average_precision_score
    >>> y_true = np.array([[1, 0, 0], [0, 0, 1]])
    >>> y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
    >>> label_ranking_average_precision_score(y_true, y_score) # doctest: +ELLIPSIS
    0.416...

.. _label_ranking_loss:

Ranking loss
------------

The :func:`label_ranking_loss` function computes the ranking loss which
averages over the samples the number of label pairs that are incorrectly
ordered, i.e. true labels have a lower score than false labels, weighted by the
the inverse number of false and true labels. The lowest achievable
ranking loss is zero.

Formally, given a binary indicator matrix of the ground truth labels
:math:`y \in \left\{0, 1\right\}^{n_\text{samples} \times n_\text{labels}}` and the
score associated with each label
:math:`\hat{f} \in \mathbb{R}^{n_\text{samples} \times n_\text{labels}}`,
the ranking loss is defined as

.. math::
  \text{ranking\_loss}(y, \hat{f}) =  \frac{1}{n_{\text{samples}}}
    \sum_{i=0}^{n_{\text{samples}} - 1} \frac{1}{|y_i|(n_\text{labels} - |y_i|)}
    \left|\left\{(k, l): \hat{f}_{ik} < \hat{f}_{il}, y_{ik} = 1, y_{il} = 0 \right\}\right|

where :math:`|\cdot|` is the :math:`\ell_0` norm or the cardinality of the set.

Here is a small example of usage of this function::

    >>> import numpy as np
    >>> from sklearn.metrics import label_ranking_loss
    >>> y_true = np.array([[1, 0, 0], [0, 0, 1]])
    >>> y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
    >>> label_ranking_loss(y_true, y_score) # doctest: +ELLIPSIS
    0.75...
    >>> # With the following prediction, we have perfect and minimal loss
    >>> y_score = np.array([[1.0, 0.1, 0.2], [0.1, 0.2, 0.9]])
    >>> label_ranking_loss(y_true, y_score)
    0.0

.. _regression_metrics:

Regression metrics
===================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` module implements several loss, score, and utility
functions to measure regression performance. Some of those have been enhanced
to handle the multioutput case: :func:`mean_squared_error`,
:func:`mean_absolute_error`, :func:`explained_variance_score` and
:func:`r2_score`.


These functions have an ``multioutput`` keyword argument which specifies the
way the scores or losses for each individual target should be averaged. The
default is ``'uniform_average'``, which specifies a uniformly weighted mean
over outputs. If an ``ndarray`` of shape ``(n_outputs,)`` is passed, then its
entries are interpreted as weights and an according weighted average is
returned. If ``multioutput`` is ``'raw_values'`` is specified, then all
unaltered individual scores or losses will be returned in an array of shape
``(n_outputs,)``.


The :func:`r2_score` and :func:`explained_variance_score` accept an additional
value ``'variance_weighted'`` for the ``multioutput`` parameter. This option
leads to a weighting of each individual score by the variance of the
corresponding target variable. This setting quantifies the globally captured
unscaled variance. If the target variables are of different scale, then this
score puts more importance on well explaining the higher variance variables.
``multioutput='variance_weighted'`` is the default value for :func:`r2_score`
for backward compatibility. This will be changed to ``uniform_average`` in the
future.

.. _explained_variance_score:

Explained variance score
-------------------------

The :func:`explained_variance_score` computes the `explained variance
regression score <https://en.wikipedia.org/wiki/Explained_variation>`_.

If :math:`\hat{y}` is the estimated target output, :math:`y` the corresponding
(correct) target output, and :math:`Var` is `Variance
<https://en.wikipedia.org/wiki/Variance>`_, the square of the standard deviation,
then the explained variance is estimated as follow:

.. math::

  \texttt{explained\_{}variance}(y, \hat{y}) = 1 - \frac{Var\{ y - \hat{y}\}}{Var\{y\}}

The best possible score is 1.0, lower values are worse.

Here is a small example of usage of the :func:`explained_variance_score`
function::

    >>> from sklearn.metrics import explained_variance_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> explained_variance_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.957...
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> explained_variance_score(y_true, y_pred, multioutput='raw_values')
    ... # doctest: +ELLIPSIS
    array([ 0.967...,  1.        ])
    >>> explained_variance_score(y_true, y_pred, multioutput=[0.3, 0.7])
    ... # doctest: +ELLIPSIS
    0.990...

.. _mean_absolute_error:

Mean absolute error
-------------------

The :func:`mean_absolute_error` function computes `mean absolute
error <https://en.wikipedia.org/wiki/Mean_absolute_error>`_, a risk
metric corresponding to the expected value of the absolute error loss or
:math:`l1`-norm loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample,
and :math:`y_i` is the corresponding true value, then the mean absolute error
(MAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MAE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}-1} \left| y_i - \hat{y}_i \right|.

Here is a small example of usage of the :func:`mean_absolute_error` function::

  >>> from sklearn.metrics import mean_absolute_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> mean_absolute_error(y_true, y_pred)
  0.5
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> mean_absolute_error(y_true, y_pred)
  0.75
  >>> mean_absolute_error(y_true, y_pred, multioutput='raw_values')
  array([ 0.5,  1. ])
  >>> mean_absolute_error(y_true, y_pred, multioutput=[0.3, 0.7])
  ... # doctest: +ELLIPSIS
  0.849...

.. _mean_squared_error:

Mean squared error
-------------------

The :func:`mean_squared_error` function computes `mean square
error <https://en.wikipedia.org/wiki/Mean_squared_error>`_, a risk
metric corresponding to the expected value of the squared (quadratic) error loss or
loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample,
and :math:`y_i` is the corresponding true value, then the mean squared error
(MSE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MSE}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (y_i - \hat{y}_i)^2.

Here is a small example of usage of the :func:`mean_squared_error`
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

.. _median_absolute_error:

Median absolute error
---------------------

The :func:`median_absolute_error` is particularly interesting because it is
robust to outliers. The loss is calculated by taking the median of all absolute
differences between the target and the prediction.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the median absolute error
(MedAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MedAE}(y, \hat{y}) = \text{median}(\mid y_1 - \hat{y}_1 \mid, \ldots, \mid y_n - \hat{y}_n \mid).

The :func:`median_absolute_error` does not support multioutput.

Here is a small example of usage of the :func:`median_absolute_error`
function::

  >>> from sklearn.metrics import median_absolute_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> median_absolute_error(y_true, y_pred)
  0.5

.. _r2_score:

R² score, the coefficient of determination
-------------------------------------------

The :func:`r2_score` function computes R², the `coefficient of
determination <https://en.wikipedia.org/wiki/Coefficient_of_determination>`_.
It provides a measure of how well future samples are likely to
be predicted by the model. Best possible score is 1.0 and it can be negative
(because the model can be arbitrarily worse). A constant model that always
predicts the expected value of y, disregarding the input features, would get a
R^2 score of 0.0.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the score R² estimated
over :math:`n_{\text{samples}}` is defined as

.. math::

  R^2(y, \hat{y}) = 1 - \frac{\sum_{i=0}^{n_{\text{samples}} - 1} (y_i - \hat{y}_i)^2}{\sum_{i=0}^{n_\text{samples} - 1} (y_i - \bar{y})^2}

where :math:`\bar{y} =  \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}} - 1} y_i`.

Here is a small example of usage of the :func:`r2_score` function::

  >>> from sklearn.metrics import r2_score
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
  0.948...
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> r2_score(y_true, y_pred, multioutput='variance_weighted')
  ... # doctest: +ELLIPSIS
  0.938...
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> r2_score(y_true, y_pred, multioutput='uniform_average')
  ... # doctest: +ELLIPSIS
  0.936...
  >>> r2_score(y_true, y_pred, multioutput='raw_values')
  ... # doctest: +ELLIPSIS
  array([ 0.965...,  0.908...])
  >>> r2_score(y_true, y_pred, multioutput=[0.3, 0.7])
  ... # doctest: +ELLIPSIS
  0.925...


.. topic:: Example:

  * See :ref:`example_linear_model_plot_lasso_and_elasticnet.py`
    for an example of R² score usage to
    evaluate Lasso and Elastic Net on sparse signals.

.. _clustering_metrics:

Clustering metrics
======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` module implements several loss, score, and utility
functions. For more information see the :ref:`clustering_evaluation`
section for instance clustering, and :ref:`biclustering_evaluation` for
biclustering.


.. _dummy_estimators:


Dummy estimators
=================

.. currentmodule:: sklearn.dummy

When doing supervised learning, a simple sanity check consists of comparing
one's estimator against simple rules of thumb. :class:`DummyClassifier`
implements several such simple strategies for classification:

- ``stratified`` generates random predictions by respecting the training
  set class distribution.
- ``most_frequent`` always predicts the most frequent label in the training set.
- ``prior`` always predicts the class that maximizes the class prior
  (like ``most_frequent`) and ``predict_proba`` returns the class prior.
- ``uniform`` generates predictions uniformly at random.
- ``constant`` always predicts a constant label that is provided by the user.
   A major motivation of this method is F1-scoring, when the positive class
   is in the minority.

Note that with all these strategies, the ``predict`` method completely ignores
the input data!

To illustrate :class:`DummyClassifier`, first let's create an imbalanced
dataset::

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.model_selection import train_test_split
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

We see that the accuracy was boosted to almost 100%.  A cross validation
strategy is recommended for a better estimate of the accuracy, if it
is not too CPU costly. For more information see the :ref:`cross_validation`
section. Moreover if you want to optimize over the parameter space, it is highly
recommended to use an appropriate methodology; see the :ref:`grid_search`
section for details.

More generally, when the accuracy of a classifier is too close to random, it
probably means that something went wrong: features are not helpful, a
hyperparameter is not correctly tuned, the classifier is suffering from class
imbalance, etc...

:class:`DummyRegressor` also implements four simple rules of thumb for regression:

- ``mean`` always predicts the mean of the training targets.
- ``median`` always predicts the median of the training targets.
- ``quantile`` always predicts a user provided quantile of the training targets.
- ``constant`` always predicts a constant value that is provided by the user.

In all these strategies, the ``predict`` method completely ignores
the input data.
