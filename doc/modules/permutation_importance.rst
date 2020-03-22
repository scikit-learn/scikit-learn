
.. _permutation_importance:

Permutation feature importance
==============================

.. currentmodule:: sklearn.inspection

Permutation feature importance is a model inspection technique that can be used
for any :term:`fitted` :term:`estimator` when the data is tabular. This is
especially useful for non-linear or opaque :term:`estimators`. The permutation
feature importance is defined to be the decrease in a model score when a single
feature value is randomly shuffled [1]_. This procedure breaks the relationship
between the feature and the target, thus the drop in the model score is
indicative of how much the model depends on the feature. This technique
benefits from being model agnostic and can be calculated many times with
different permutations of the feature.

The :func:`permutation_importance` function calculates the feature importance
of :term:`estimators` for a given dataset. The ``n_repeats`` parameter sets the
number of times a feature is randomly shuffled and returns a sample of feature
importances.

Let's consider the following trained regression model::

  >>> from sklearn.datasets import load_diabetes
  >>> from sklearn.model_selection import train_test_split
  >>> from sklearn.linear_model import Ridge
  >>> diabetes = load_diabetes()
  >>> X_train, X_val, y_train, y_val = train_test_split(
  ...     diabetes.data, diabetes.target, random_state=0)
  ...
  >>> model = Ridge(alpha=1e-2).fit(X_train, y_train)
  >>> model.score(X_val, y_val)
  0.356...

Its validation performance, measured via the :math:`R^2` score, is
significantly larger than the chance level. This makes it possible to use the
:func:`permutation_importance` function to probe which features are most
predictive::

  >>> from sklearn.inspection import permutation_importance
  >>> r = permutation_importance(model, X_val, y_val,
  ...                            n_repeats=30,
  ...                            random_state=0)
  ...
  >>> for i in r.importances_mean.argsort()[::-1]:
  ...     if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
  ...         print(f"{diabetes.feature_names[i]:<8}"
  ...               f"{r.importances_mean[i]:.3f}"
  ...               f" +/- {r.importances_std[i]:.3f}")
  ...
  s5      0.204 +/- 0.050
  bmi     0.176 +/- 0.048
  bp      0.088 +/- 0.033
  sex     0.056 +/- 0.023

Note that the importance values for the top features represent a large
fraction of the reference score of 0.356.

Permutation importances can be computed either on the training set or on a
held-out testing or validation set. Using a held-out set makes it possible to
highlight which features contribute the most to the generalization power of the
inspected model. Features that are important on the training set but not on the
held-out set might cause the model to overfit.

.. warning::

  Features that are deemed of **low importance for a bad model** (low
  cross-validation score) could be **very important for a good model**.
  Therefore it is always important to evaluate the predictive power of a model
  using a held-out set (or better with cross-validation) prior to computing
  importances. Permutation importance does not reflect to the intrinsic
  predictive value of a feature by itself but **how important this feature is
  for a particular model**.

Outline of the permutation importance algorithm
-----------------------------------------------

- Inputs: fitted predictive model :math:`m`, tabular dataset (training or
  validation) :math:`D`.
- Compute the reference score :math:`s` of the model :math:`m` on data
  :math:`D` (for instance the accuracy for a classifier or the :math:`R^2` for
  a regressor).
- For each feature :math:`j` (column of :math:`D`):

  - For each repetition :math:`k` in :math:`{1, ..., K}`:

    - Randomly shuffle column :math:`j` of dataset :math:`D` to generate a
      corrupted version of the data named :math:`\tilde{D}_{k,j}`.
    - Compute the score :math:`s_{k,j}` of model :math:`m` on corrupted data
      :math:`\tilde{D}_{k,j}`.

  - Compute importance :math:`i_j` for feature :math:`f_j` defined as:

    .. math:: i_j = s - \frac{1}{K} \sum_{k=1}^{K} s_{k,j}

Relation to impurity-based importance in trees
----------------------------------------------

Tree-based models provide an alternative measure of :ref:`feature importances
based on the mean decrease in impurity <random_forest_feature_importance>`
(MDI). Impurity is quantified by the splitting criterion of the decision trees
(Gini, Entropy or Mean Squared Error). However, this method can give high
importance to features that may not be predictive on unseen data when the model
is overfitting. Permutation-based feature importance, on the other hand, avoids
this issue, since it can be computed on unseen data.

Furthermore, impurity-based feature importance for trees are **strongly
biased** and **favor high cardinality features** (typically numerical features)
over low cardinality features such as binary features or categorical variables
with a small number of possible categories.

Permutation-based feature importances do not exhibit such a bias. Additionally,
the permutation feature importance may be computed performance metric on the
model predictions predictions and can be used to analyze any model class (not
just tree-based models).

The following example highlights the limitations of impurity-based feature
importance in contrast to permutation-based feature importance:
:ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`.

Misleading values on strongly correlated features
-------------------------------------------------

When two features are correlated and one of the features is permuted, the model
will still have access to the feature through its correlated feature. This will
result in a lower importance value for both features, where they might
*actually* be important.

One way to handle this is to cluster features that are correlated and only
keep one feature from each cluster. This strategy is explored in the following
example:
:ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`.

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`
  * :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance_multicollinear.py`

.. topic:: References:

   .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
       2001. https://doi.org/10.1023/A:1010933404324
