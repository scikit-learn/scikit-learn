.. currentmodule:: sklearn.model_selection

.. _cutoffclassifier:

========================================================
Tuning cut-off decision threshold for classes prediction
========================================================

Classifiers are predictive models: they use statistical learning to predict
outcomes. The outcomes of a classifier takes two forms: a "soft" score for each
sample in relation to each class, and a "hard" categorical prediction (i.e.
class label). Soft predictions are obtained using :term:`predict_proba` or
:term:`decision_function` while hard predictions are obtained using
:term:`predict`.

In scikit-learn, there is a connection between soft and hard prediction. In the
case of a binary classification, hard predictions are obtained by associating
the positive class with probability value greater than 0.5 (obtained with
:term:`predict_proba`) or decision function value greater than 0 (obtained with
:term:`decision_function`).

    >>> from sklearn.datasets import make_classification
    >>> from sklearn.tree import DecisionTreeClassifier
    >>> X, y = make_classification(random_state=0)
    >>> classifier = DecisionTreeClassifier(max_depth=2, random_state=0).fit(X, y)
    >>> classifier.predict_proba(X[:4])
    array([[0.94   , 0.06   ],
           [0.94   , 0.06   ],
           [0.04..., 0.95...],
           [0.04..., 0.95...]])
    >>> classifier.predict(X[:4])
    array([0, 0, 1, 1])


Similar rules apply for other classification problems.

While these approaches are reasonable as default behaviors, they might not be
adapted to some cases. The context and nature of the use case define the
expected behavior of the classifier and thus the strategy to convert soft
predictions into hard predictions. We illustrate this point with an example.

Let's imagine the deployment of a predictive model helping medical doctors to
detect cancers. In a setting where this model would be a tool to discard
obvious cases, doctors might be interested to have a high recall (all cancers
cases should be tagged as such) to not miss any patient with a cancer. However,
it will be at the cost of having more false positive predictions (i.e. lower
precision). Thus, in terms of decision threshold, it would be better to
classify a patient having a cancer for a lower probability than 0.5.

Post-tuning of the decision threshold
=====================================

One solution to address the problem stated in the introduction is to tune the decision
threshold of the classifier once this model has been trained. The
:class:`~sklearn.model_selection.CutOffClassifier` allows to tune this threshold using
an internal cross-validation. The optimum threshold is tuned to maximize a given metric
with or without constraints.

The following image illustrate the tuning of the cut-off point for a gradient
boosting classifier. While the vanilla and tuned classifiers provide the same
Receiver Operating Characteristic (ROC) and Precision-Recall curves, and thus
the same :term:`predict_proba` outputs, the "hard" predictions defer because of
the tuned cut-off point. The vanilla classifier predicts the class of interest
for a probability greater than 0.5 while the tuned classifier predicts the
class of interest for a very low probability (around 0.02). This cut-off point
is maximizes a utility metric defined by the business case (in this case an
insurance company).

.. figure:: ../auto_examples/model_selection/images/sphx_glr_plot_cutoff_tuning_002.png
   :target: ../auto_examples/model_selection/plot_cutoff_tuning.html
   :align: center

Available options to tune the cut-off point
-------------------------------------------

The cut-off point can be tuned with different strategies controlled by the parameter
`objective_metric`.

A straightforward use case is to maximize a pre-defined scikit-learn metric. These
metrics can be found by calling the function :func:`~sklearn.metrics.get_scorer_names`.
We provide an example where we maximize the balanced accuracy.

.. note::

    It is important to notice that these metrics comes with default parameter, notably
    the label of the class of interested (i.e. `pos_label`). Thus, if this label is not
    the right one for your application, you need to define a scorer and pass the right
    `pos_label` (and additional parameters) using the
    :func:`~sklearn.metrics.make_scorer`. You should refer to :ref:`scoring` to get all
    information to define your own scoring function. For instance, we show how to pass
    the information to the scorer that the label of interest is `0` when maximizing the
    :func:`~sklearn.metrics.f1_score`:

        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.metrics import make_scorer, f1_score
        >>> X, y = make_classification(
        ...    n_samples=1_000, weights=[0.1, 0.9], random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
        >>> pos_label = 0
        >>> scorer = make_scorer(f1_score, pos_label=pos_label)
        >>> model = CutOffClassifier(classifier, objective_metric=scorer).fit(
        ...     X_train, y_train)
        >>> scorer(model, X_test, y_test)
        0.82...
        >>> # compare it with the internal score found by cross-validation
        >>> model.objective_score_
        0.86...

A second strategy aims at maximizing a metric while imposing constraints on another
metric. Four pre-defined options exist, 2 that uses the Receiver Operating
Characteristic (ROC) statistic and 2 that uses the Precision-Recall statistic.

- `"max_tpr_at_tnr_constraint"`: maximizes the True Positive Rate (TPR) such that the
  True Negative Rate (TNR) is the closest to a given value.
- `"max_tnr_at_tpr_constraint"`: maximizes the TNR such that the TPR is the closest to
  a given value.
- `"max_precision_at_recall_constraint"`: maximizes the precision such that the recall
    is the closest to a given value.
- `"max_recall_at_precision_constraint"`: maximizes the recall such that the precision
    is the closest to a given value.

For these options, the `constraint_value` parameter needs to be defined. In addition,
you can use the `pos_label` parameter to indicate the label of the class of interest.

The final strategy maximizes a custom utility function. This problem is also known as
cost-sensitive learning. The utility function is defined by providing dictionary
containing the cost-gain associated with the entries of the confusion matrix. The keys
are defined as `{"tn", "fp", "fn", "tp"}`. The class of interest is defined using the
`pos_label` parameter. Refer to :ref:`cost_sensitive_learning_example` for an example
depicting the use of such a utility function.

Important notes regarding the internal cross-validation
-------------------------------------------------------

By default :class:`~sklearn.model_selection.CutOffClassifier` uses a 5-fold stratified
cross-validation to tune the cut-off point. The parameter `cv` allows to control the
cross-validation strategy. It is possible to go around cross-validation by passing
`cv="prefit"` and provide an already fitted classifier. In this case, the cut-off point
is tuned on the data provided to the `fit` method.

However, you should be extremely careful when using this option. You should never use
the same data for training the classifier and tuning the cut-off point at the risk of
overfitting. Refer to :ref:`cutoffclassifier_no_cv` that shows such overfitting. If
you are in a situation where you have limited resources, you should can consider using
a float number that will use a single split internally.

The option `cv="prefit"` should only be used when the provided classifier was already
trained on some data and you want to tune (or re-tune) on a new validation set.

Examples
--------

- See :ref:`sphx_glr_auto_examples_model_selection_plot_cutoff_tuning.py` example for
  an example of tuning the decision threshold of a classifier.
