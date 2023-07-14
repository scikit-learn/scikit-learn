.. currentmodule:: sklearn.model_selection

.. _tunedthresholdclassifier:

========================================================
Tuning cut-off decision threshold for classes prediction
========================================================

Classifiers are predictive models: they use statistical learning to predict outcomes.
The outcomes of a classifier are scores for each sample in relation to each class and
categorical prediction (class label). Scores are obtained from :term:`predict_proba` or
:term:`decision_function`. The former returns posterior probability estimates for each
class, while the latter returns a decision score for each class. The decision score is a
measure of how strongly the sample is predicted to belong to the positive class (e.g.,
the distance to the decision boundary). A decision rule is then defined by thresholding
the scores and obtained the class label for each sample. Those labels are obtained with
:term:`predict`.

For binary classification in scikit-learn, class labels are obtained by associating the
positive class with posterior probability estimates greater than 0.5 (obtained with
:term:`predict_proba`) or decision scores greater than 0 (obtained with
:term:`decision_function`).

Here, we show an example that illustrates the relation between posterior
probability estimates and class labels::

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

While these approaches are reasonable as default behaviors, they are not ideal for
all cases. The context and nature of the use case defines the expected behavior of the
classifier and thus, the strategy to convert soft predictions into hard predictions. We
illustrate this point with an example.

Let's imagine the deployment of a predictive model helping medical doctors to detect
tumour. In a setting where this model was a tool to discard obvious cases and false
positives don't lead to potentially harmful treatments, doctors might be interested in
having a high recall (all cancer cases should be tagged as such) to not miss any patient
with a cancer. However, that is at the cost of having more false positive predictions
(i.e. lower precision). Thus, in terms of decision threshold, it may be better to
classify a patient as having a cancer for a posterior probability estimate lower than
0.5.

Post-tuning of the decision threshold
=====================================

One solution to address the problem stated in the introduction is to tune the decision
threshold of the classifier once the model has been trained. The
:class:`~sklearn.model_selection.TunedThresholdClassifier` tunes this threshold using an
internal cross-validation. The optimum threshold is chosen to maximize a given metric
with or without constraints.

The following image illustrates the tuning of the cut-off point for a gradient boosting
classifier. While the vanilla and tuned classifiers provide the same Receiver Operating
Characteristic (ROC) and Precision-Recall curves, and thus the same
:term:`predict_proba` outputs, the class label predictions differ because of the tuned
decision threshold. The vanilla classifier predicts the class of interest for a
posterior probability greater than 0.5 while the tuned classifier predicts the class of
interest for a very low probability (around 0.02). This cut-off point optimizes a
utility metric defined by the business (in this case an insurance company).

.. figure:: ../auto_examples/model_selection/images/sphx_glr_plot_tuned_threshold_classifier_002.png
   :target: ../auto_examples/model_selection/plot_tuned_threshold_classifier.html
   :align: center

Options to tune the cut-off point
-------------------------------------------

The cut-off point can be tuned through different strategies controlled by the parameter
`objective_metric`.

One way to tune the threshold is by maximizing a pre-defined scikit-learn metric. These
metrics can be found by calling the function :func:`~sklearn.metrics.get_scorer_names`.
In this example, we maximize the balanced accuracy.

.. note::

    It is important to notice that these metrics come with default parameters, notably
    the label of the class of interested (i.e. `pos_label`). Thus, if this label is not
    the right one for your application, you need to define a scorer and pass the right
    `pos_label` (and additional parameters) using the
    :func:`~sklearn.metrics.make_scorer`. Refer to :ref:`scoring` to get 
    information to define your own scoring function. For instance, we show how to pass
    the information to the scorer that the label of interest is `0` when maximizing the
    :func:`~sklearn.metrics.f1_score`:

        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.model_selection import (
        ...     TunedThresholdClassifier, train_test_split
        ... )
        >>> from sklearn.metrics import make_scorer, f1_score
        >>> X, y = make_classification(
        ...    n_samples=1_000, weights=[0.1, 0.9], random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
        >>> pos_label = 0
        >>> scorer = make_scorer(f1_score, pos_label=pos_label)
        >>> base_model = LogisticRegression()
        >>> model = TunedThresholdClassifier(base_model, objective_metric=scorer).fit(
        ...     X_train, y_train)
        >>> scorer(model, X_test, y_test)
        0.82...
        >>> # compare it with the internal score found by cross-validation
        >>> model.objective_score_
        0.86...

A second strategy aims to maximize one metric while imposing constraints on another
metric. There are four pre-defined options, 2 use the Receiver Operating
Characteristic (ROC) statistics and 2 use the Precision-Recall statistics.

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
cost-sensitive learning. The utility function is defined by providing a dictionary
containing the cost-gain associated with the entries of the confusion matrix. The keys
are defined as `{"tn", "fp", "fn", "tp"}`. The class of interest is defined using the
`pos_label` parameter. Refer to :ref:`cost_sensitive_learning_example` for an example
depicting the use of such a utility function.

Important notes regarding the internal cross-validation
-------------------------------------------------------

By default :class:`~sklearn.model_selection.TunedThresholdClassifier` uses a
5-fold stratified cross-validation to tune the cut-off point. The parameter
`cv` allows to control the cross-validation strategy. It is possible to go
around cross-validation by passing `cv="prefit"` and providing a fitted
classifier. In this case, the cut-off point is tuned on the data provided to
the `fit` method.

However, you should be extremely careful when using this option. You should never use
the same data for training the classifier and tuning the cut-off point at the risk of
overfitting. Refer to :ref:`tunedthresholdclassifier_no_cv` for an example. If
you are in a situation where you have limited resources, you should consider using
a float number that will use a single split internally.

The option `cv="prefit"` should only be used when the provided classifier was already
trained, and you just want to find the best cut-off using a new validation set.

Manually setting the decision thresholding
-------------------------------------------

The previous sections discussed strategies to find an optimal decision threshold. It is
also possible to manually set the decision threshold in
:class`~sklearn.model_selection.TunedThresholdClassifier` by setting the parameter
`strategy` to `"constant"` and providing the desired threshold using the parameter
`constant_threshold`.

Examples
--------

- See
  :ref:`sphx_glr_auto_examples_model_selection_plot_tuned_threshold_classifier.py`
  example for an example of tuning the decision threshold of a classifier.
