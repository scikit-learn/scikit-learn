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
adapted to some cases. The context and nature of the use case will define the
expected behavior of the classifier and thus the strategy to convert soft
predictions into hard predictions. We will illustrate this point with an
example.

Let's imagine the deployment of a predictive model helping medical doctors to
detect cancers. In a setting where this model would be a tool to discard
obvious cases, doctors might be interested to have a high recall (all cancers
cases should be tagged as such) to not miss any patient with a cancer. However,
it will be at the cost of having more false positive predictions (i.e. lower
precision). Thus, in terms of decision threshold, it would be better to
classify a patient having a cancer for a lower probability than 0.5.

Post-tuning of the decision threshold
=====================================

One solution to address the problem stated in the introduction is to tune the
decision threshold of the classifier once this model has been trained. The
:class:`CutOffClassifier` allows to tune this threshold using an internal
cross-validation.
