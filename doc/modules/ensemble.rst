
.. _forest:

RandomForest
============

In random forests (see :class:`RandomForestClassifier` and
:class:`RandomForestRegressor` classes), each tree in the ensemble is built
from a sample drawn with replacement (i.e., a bootstrap sample) from the
training set.

Furthermore, when splitting each node during the construction of a tree, the
best split is found either from all input features or a random subset of size
``max_features``. (See the :ref:`parameter description <forest_parameters>` for
more details).

The purpose of these two sources of randomness is to decrease the variance of
the forest estimator. Indeed, individual decision trees typically exhibit high
variance and tend to overfit. The injected randomness in forests yields
decision trees with somewhat decoupled prediction errors. By taking an average
of those predictions, some errors can cancel out. Random forests achieve a
reduced variance by combining diverse trees, sometimes at the cost of a slight
increase in bias. In practice, the variance reduction is often significant
hence yielding an overall better model.

In contrast to the original publication [B2001]_, the scikit-learn
implementation combines classifiers by averaging their probabilistic
prediction, instead of letting each classifier vote for a single class.

A competitive alternative to random forests is the :ref:`Extra-Trees method
<extra_trees>` (see :class:`ExtraTreesClassifier` and
:class:`ExtraTreesRegressor` classes), which takes randomness one step further:
the splitting thresholds are randomized. Instead of looking for the most
discriminative threshold, thresholds are drawn at random for each candidate
feature and the best of these randomly-generated thresholds is picked as the
splitting rule. This usually allows to reduce the variance of the model a bit
more, at the cost of slightly increasing the bias.

.. _forest_parameters:

Parameters
----------

The main parameters to adjust when using these methods are ``n_estimators`` and
``max_features``. The former is the number of trees in the forest. The larger
the better, but also the longer it will take to compute. In addition, note that
results will stop getting significantly better beyond a critical number of
trees. The latter is the size of the random subsets of features to consider
when splitting a node. The lower the greater the reduction of variance, but also
the greater the increase in bias. Empirical good default values are
``max_features=sqrt(n_features)`` for classification tasks (where
``n_features`` is the number of features in the data) and
``max_features=n_features/3`` for regression tasks.

Good results are often achieved by setting ``max_depth=None`` in combination
with ``min_samples_split=2`` (i.e., when fully developing the trees). Bear in
mind though that these values are usually not optimal, and they might result in
models that consume a lot of RAM. The best parameter values should always be
cross-validated. In addition, note that in random forests, bootstrap samples
are used by default (``bootstrap=True``) while the default strategy for
extra-trees is to use the whole dataset (``bootstrap=False``). When using
bootstrap sampling the generalization accuracy can be estimated on the left out
or *out of bag* samples. This can be enabled by setting ``oob_score=True``.

Another important parameter is the criterion used to split nodes. For both
classification and regression, ``"gini"`` and ``"entropy"`` are supported for
classification, while ``"squared_error"`` (and the deprecated
``"mse"``) and ``"mae"`` are supported for regression. The default for both
Random Forests and Extra-Trees is ``"gini"`` for classification and
``"squared_error"`` for regression.

When the goal is to reduce memory usage, the sparsity of the data can be
exploited: sparse matrices can be accepted as input, as long as they are
formatted as ``csc_matrix``. Additionally, ``max_features`` can be lowered to
reduce memory usage, as fewer features are evaluated at each split.

The parameter ``splitter`` controls the strategy used to choose the split at
each node. Supported strategies are ``"best"`` to choose the best split and
``"random"`` to choose the best random split. When ``"best"`` is used (the
default), the splitter performs a brute-force exhaustive search of possible
split midpoints to find the optimal split for each feature considered at each
node.

.. note::

   The default ``max_features`` differs for :class:`RandomForestClassifier`
   (``sqrt(n_features)``) and :class:`RandomForestRegressor`
   (``n_features``). See the class documentation for details.

Finally, this module also features the parallel construction of the trees and
the parallel computation of the predictions through the ``n_jobs`` parameter.
If ``n_jobs=k`` then computations are partitioned into ``k`` jobs, and run on
``k`` cores of the machine. If ``n_jobs=-1`` then all available cores are
used. Note that because of inter-process communication overhead, the speedup
might not be linear (i.e., using ``k`` jobs will unfortunately not be ``k``
times as fast). Significant speedup can still be achieved though when building
a large number of trees, or when building a single tree requires a fair amount
of work (e.g., on large datasets).

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_ensemble_plot_forest_iris.py`
* :ref:`sphx_glr_auto_examples_ensemble_plot_forest_importances_faces.py`
* :ref:`sphx_glr_auto_examples_ensemble_plot_bias_variance.py`

.. rubric:: References

.. [B2001] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

.. [B1998] L. Breiman, "Arcing Classifiers", Annals of Statistics 1998.

* P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
  Machine Learning, 63(1), 3-42, 2006.

