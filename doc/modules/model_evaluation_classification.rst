.. currentmodule:: sklearn.metrics

.. _classification_metrics:

Classification metrics
======================

The :mod:`sklearn.metrics` module implements several loss, score, and utility
functions to measure classification performance.
Some metrics might require probability estimates of the positive class,
confidence values, or binary decisions values.
Most implementations allow each sample to provide a weighted contribution
to the overall score, through the ``sample_weight`` parameter.

Some of these are restricted to the binary classification case:

.. autosummary::
   :template: function.rst

   precision_recall_curve
   roc_curve


Others also work in the multiclass case:

.. autosummary::
   :template: function.rst

   cohen_kappa_score
   confusion_matrix
   hinge_loss
   matthews_corrcoef


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
  metric per class, this sums the dividends and divisors that make up the
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

  * See :ref:`sphx_glr_auto_examples_feature_selection_plot_permutation_test_for_classification.py`
    for an example of accuracy score usage using permutations of
    the dataset.

.. _cohen_kappa:

Cohen's kappa
-------------

The function :func:`cohen_kappa_score` computes `Cohen's kappa
<https://en.wikipedia.org/wiki/Cohen%27s_kappa>`_ statistic.
This measure is intended to compare labelings by different human annotators,
not a classifier versus a ground truth.

The kappa score (see docstring) is a number between -1 and 1.
Scores above .8 are generally considered good agreement;
zero or lower means no agreement (practically random labels).

Kappa scores can be computed for binary or multiclass problems,
but not for multilabel problems (except by manually computing a per-label score)
and not for more than two annotators.

  >>> from sklearn.metrics import cohen_kappa_score
  >>> y_true = [2, 0, 2, 2, 0, 1]
  >>> y_pred = [0, 0, 2, 2, 0, 2]
  >>> cohen_kappa_score(y_true, y_pred)
  0.4285714285714286

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
from the :ref:`sphx_glr_auto_examples_model_selection_plot_confusion_matrix.py` example):

.. image:: ../auto_examples/model_selection/images/sphx_glr_plot_confusion_matrix_001.png
   :target: ../auto_examples/model_selection/plot_confusion_matrix.html
   :scale: 75
   :align: center

For binary problems, we can get counts of true negatives, false positives,
false negatives and true positives as follows::

  >>> y_true = [0, 0, 0, 1, 1, 1, 1, 1]
  >>> y_pred = [0, 1, 0, 1, 0, 1, 0, 1]
  >>> tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
  >>> tn, fp, fn, tp
  (2, 1, 2, 3)

.. topic:: Example:

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_confusion_matrix.py`
    for an example of using a confusion matrix to evaluate classifier output
    quality.

  * See :ref:`sphx_glr_auto_examples_classification_plot_digits_classification.py`
    for an example of using a confusion matrix to classify
    hand-written digits.

  * See :ref:`sphx_glr_auto_examples_text_document_classification_20newsgroups.py`
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
   >>> y_pred = [0, 0, 2, 1, 0]
   >>> target_names = ['class 0', 'class 1', 'class 2']
   >>> print(classification_report(y_true, y_pred, target_names=target_names))
                precision    recall  f1-score   support
   <BLANKLINE>
       class 0       0.67      1.00      0.80         2
       class 1       0.00      0.00      0.00         1
       class 2       1.00      0.50      0.67         2
   <BLANKLINE>
   avg / total       0.67      0.60      0.59         5
   <BLANKLINE>

.. topic:: Example:

  * See :ref:`sphx_glr_auto_examples_classification_plot_digits_classification.py`
    for an example of classification report usage for
    hand-written digits.

  * See :ref:`sphx_glr_auto_examples_text_document_classification_20newsgroups.py`
    for an example of classification report usage for text
    documents.

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_digits.py`
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

  * See :ref:`sphx_glr_auto_examples_text_document_classification_20newsgroups.py`
    for an example of :func:`f1_score` usage to classify  text
    documents.

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_digits.py`
    for an example of :func:`precision_score` and :func:`recall_score` usage
    to estimate parameters using grid search with nested cross-validation.

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_precision_recall.py`
    for an example of :func:`precision_recall_curve` usage to evaluate
    classifier output quality.

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
  0.83...



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


In the binary (two-class) case, :math:`tp`, :math:`tn`, :math:`fp` and
:math:`fn` are respectively the number of true positives, true negatives, false
positives and false negatives, the MCC is defined as

.. math::

  MCC = \frac{tp \times tn - fp \times fn}{\sqrt{(tp + fp)(tp + fn)(tn + fp)(tn + fn)}}.

In the multiclass case, the Matthews correlation coefficient can be `defined
<http://rk.kvl.dk/introduction/index.html>`_ in terms of a
:func:`confusion_matrix` :math:`C` for :math:`K` classes.  To simplify the
definition consider the following intermediate variables:

* :math:`t_k=\sum_{i}^{K} C_{ik}` the number of times class :math:`k` truly occurred,
* :math:`p_k=\sum_{i}^{K} C_{ki}` the number of times class :math:`k` was predicted,
* :math:`c=\sum_{k}^{K} C_{kk}` the total number of samples correctly predicted,
* :math:`s=\sum_{i}^{K} \sum_{j}^{K} C_{ij}` the total number of samples.

Then the multiclass MCC is defined as:

.. math::
    MCC = \frac{
        c \times s - \sum_{k}^{K} p_k \times t_k
    }{\sqrt{
        (s^2 - \sum_{k}^{K} p_k^2) \times
        (s^2 - \sum_{k}^{K} t_k^2)
    }}

When there are more than two labels, the value of the MCC will no longer range
between -1 and +1. Instead the minimum value will be somewhere between -1 and 0
depending on the number and distribution of ground true labels. The maximum
value is always +1.

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

.. image:: ../auto_examples/model_selection/images/sphx_glr_plot_roc_001.png
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


.. image:: ../auto_examples/model_selection/images/sphx_glr_plot_roc_002.png
   :target: ../auto_examples/model_selection/plot_roc.html
   :scale: 75
   :align: center

.. topic:: Examples:

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_roc.py`
    for an example of using ROC to
    evaluate the quality of the output of a classifier.

  * See :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py`
    for an example of using ROC to
    evaluate classifier output quality, using cross-validation.

  * See :ref:`sphx_glr_auto_examples_applications_plot_species_distribution_modeling.py`
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

  * See :ref:`sphx_glr_auto_examples_feature_selection_plot_rfe_with_cross_validation.py`
    for an example of zero one loss usage to perform recursive feature
    elimination with cross-validation.

.. _brier_score_loss:

Brier score loss
----------------

The :func:`brier_score_loss` function computes the
`Brier score <https://en.wikipedia.org/wiki/Brier_score>`_
for binary classes. Quoting Wikipedia:

    "The Brier score is a proper score function that measures the accuracy of
    probabilistic predictions. It is applicable to tasks in which predictions
    must assign probabilities to a set of mutually exclusive discrete outcomes."

This function returns a score of the mean square difference between the actual
outcome and the predicted probability of the possible outcome. The actual
outcome has to be 1 or 0 (true or false), while the predicted probability of
the actual outcome can be a value between 0 and 1.

The brier score loss is also between 0 to 1 and the lower the score (the mean
square difference is smaller), the more accurate the prediction is. It can be
thought of as a measure of the "calibration" of a set of probabilistic
predictions.

.. math::

   BS = \frac{1}{N} \sum_{t=1}^{N}(f_t - o_t)^2

where : :math:`N` is the total number of predictions, :math:`f_t` is the
predicted probability of the actual outcome :math:`o_t`.

Here is a small example of usage of this function:::

    >>> import numpy as np
    >>> from sklearn.metrics import brier_score_loss
    >>> y_true = np.array([0, 1, 1, 0])
    >>> y_true_categorical = np.array(["spam", "ham", "ham", "spam"])
    >>> y_prob = np.array([0.1, 0.9, 0.8, 0.4])
    >>> y_pred = np.array([0, 1, 1, 0])
    >>> brier_score_loss(y_true, y_prob)
    0.055
    >>> brier_score_loss(y_true, 1-y_prob, pos_label=0)
    0.055
    >>> brier_score_loss(y_true_categorical, y_prob, pos_label="ham")
    0.055
    >>> brier_score_loss(y_true, y_prob > 0.5)
    0.0


.. topic:: Example:

  * See :ref:`sphx_glr_auto_examples_calibration_plot_calibration.py`
    for an example of Brier score loss usage to perform probability
    calibration of classifiers.

.. topic:: References:

  * G. Brier, `Verification of forecasts expressed in terms of probability
    <http://docs.lib.noaa.gov/rescue/mwr/078/mwr-078-01-0001.pdf>`_,
    Monthly weather review 78.1 (1950)
