.. _tree:

==============
Decision Trees
==============

.. currentmodule:: sklearn.tree

**Decision Trees (DTs)** are a non-parametric supervised learning method used
for :ref:`classification <tree_classification>` and :ref:`regression
<tree_regression>`. The goal is to create a model that predicts the value of a
target variable by learning simple decision rules inferred from the data
features.

For instance, in the example below, decision trees learn from data to
approximate a sine curve with a set of if-then-else decision rules. The deeper
the tree, the more complex the decision rules and the fitter the model.

.. figure:: ../auto_examples/tree/images/plot_tree_regression_1.png
   :target: ../auto_examples/tree/plot_tree_regression.html
   :scale: 75
   :align: center

Some advantages of decision trees are:

    - Simple to understand and to interpret. Trees can be visualised.

    - Requires little data preparation. Other techniques often require data
      normalisation, dummy variables need to be created and blank values to
      be removed. Note however that this module does not support missing
      values.

    - The cost of using the tree (i.e., predicting data) is logarithmic in the
      number of data points used to train the tree.

    - Able to handle both numerical and categorical data. Other techniques
      are usually specialised in analysing datasets that have only one type
      of variable. See :ref:`algorithms <tree_algorithms>` for more
      information.

    - Uses a white box model. If a given situation is observable in a model,
      the explanation for the condition is easily explained by boolean logic.
      By constrast, in a black box model (e.g., in an artificial neural
      network), results may be more difficult to interpret.

    - Possible to validate a model using statistical tests. That makes it
      possible to account for the reliability of the model.

    - Performs well even if its assumptions are somewhat violated by
      the true model from which the data were generated.

The disadvantages of decision trees include:

    - Decision-tree learners can create over-complex trees that do not
      generalise the data well. This is called overfitting. Mechanisms
      such as pruning (not currently supported), setting the minimum
      number of samples required at a leaf node or setting the maximum
      depth of the tree are necessary to avoid this problem.

    - Decision trees can be unstable because small variations in the
      data might result in a completely different tree being generated.
      This problem is mitigated by using decision trees within an
      ensemble.

    - The problem of learning an optimal decision tree is known to be
      NP-complete under several aspects of optimality and even for simple
      concepts. Consequently, practical decision-tree learning algorithms
      are based on heuristic algorithms such as the greedy algorithm where
      locally optimal decisions are made at each node. Such algorithms
      cannot guarantee to return the globally optimal decision tree.  This
      can be mitigated by training multiple trees in an ensemble learner,
      where the features and samples are randomly sampled with replacement.

    - There are concepts that are hard to learn because decision trees
      do not express them easily, such as XOR, parity or multiplexer problems.

    - Decision tree learners create biased trees if some classes dominate.
      It is therefore recommended to balance the dataset prior to fitting
      with the decision tree.

.. _tree_classification:

Classification
==============

:class:`DecisionTreeClassifier` is a class capable of performing multi-class
classification on a dataset.

As other classifiers, :class:`DecisionTreeClassifier` take as input two
arrays: an array X of size [n_samples, n_features] holding the training
samples, and an array Y of integer values, size [n_samples], holding
the class labels for the training samples::

    >>> from sklearn import tree
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(X, Y)

After being fitted, the model can then be used to predict new values::

    >>> clf.predict([[2., 2.]])
    array([1])

:class:`DecisionTreeClassifier` is capable of both binary (where the
labels are [-1, 1]) classification and multiclass (where the labels are
[0, ..., K-1]) classification.

Using the Iris dataset, we can construct a tree as follows::

    >>> from sklearn.datasets import load_iris
    >>> from sklearn import tree
    >>> iris = load_iris()
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(iris.data, iris.target)

Once trained, we can export the tree in `Graphviz
<http://www.graphviz.org/>`_ format using the :func:`export_graphviz`
exporter. Below is an example export of a tree trained on the entire
iris dataset::

    >>> from StringIO import StringIO
    >>> out = StringIO()
    >>> out = tree.export_graphviz(clf, out_file=out)

.. only:: html

    .. figure:: ../images/iris.svg
       :align: center

.. only:: latex

    .. figure:: ../images/iris.pdf
       :align: center

After being fitted, the model can then be used to predict new values::

    >>> clf.predict(iris.data[0, :])
    array([0])

.. figure:: ../auto_examples/tree/images/plot_iris_1.png
   :target: ../auto_examples/tree/plot_iris.html
   :align: center
   :scale: 75

.. topic:: Examples:

 * :ref:`example_tree_plot_iris.py`

.. _tree_regression:

Regression
==========

.. figure:: ../auto_examples/tree/images/plot_tree_regression_1.png
   :target: ../auto_examples/tree/plot_tree_regression.html
   :scale: 75
   :align: center

Decision trees can also be applied to regression problems, using the
:class:`DecisionTreeRegressor` class.

As in the classification setting, the fit method will take as argument arrays X
and y, only that in this case y is expected to have floating point values
instead of integer values::

    >>> from sklearn import tree
    >>> X = [[0, 0], [2, 2]]
    >>> y = [0.5, 2.5]
    >>> clf = tree.DecisionTreeRegressor()
    >>> clf = clf.fit(X, y)
    >>> clf.predict([[1, 1]])
    array([ 0.5])


.. topic:: Examples:

 * :ref:`example_tree_plot_tree_regression.py`

Complexity
==========

In general, the run time cost to construct a balanced binary tree is
:math:`O(n_{samples}n_{features}log(n_{samples}))` and query time
:math:`O(log(n_{samples}))`.  Although the tree construction algorithm attempts
to generate balanced trees, they will not always be balanced.  Assuming that the
subtrees remain approximately balanced, the cost at each node consists of
searching through :math:`O(n_{features})` to find the feature that offers the
largest reduction in entropy.  This has a cost of
:math:`O(n_{features}n_{samples}log(n_{samples}))` at each node, leading to a
total cost over the entire trees (by summing the cost at each node) of
:math:`O(n_{features}n_{samples}^{2}log(n_{samples}))`.

Scikit-learn offers a more efficient implementation for the construction of
decision trees.  A naive implementation (as above) would recompute the class
label histograms (for classification) or the means (for regression) at for each
new split point along a given feature. By presorting the feature over all
relevant samples, and retaining a running label count, we reduce the complexity
at each node to :math:`O(n_{features}log(n_{samples}))`, which results in a
total cost of :math:`O(n_{features}n_{samples}log(n_{samples}))`.

This implementation also offers a parameter `min_density` to control the trade-
off  between memory usage and speed within the tree construction process. A
sample mask is used to mask data points that are inactive at a given node, which
avoids the copying of data (important for large datasets or training trees
within an ensemble). Density is defined as the ratio of 'active' data samples to
total samples at a given node.  The minimum density parameter specifies the
level below which fancy indexing (and therefore data copied) and the sample mask
reset.
If `min_density` is 1, then fancy indexing is always used for data partitioning
during the tree building phase. In this case, the size of memory (as a
proportion of the input data :math:`a`) required at a node of depth :math:`n`
can be approximated using a geometric series: :math:`size = a \frac{1 - r^n}{1 -
r}` where :math:`r` is the ratio of samples used at each node.  A best case
analysis shows that the lowest memory requirement (for an infinitely deep tree)
is :math:`2 \times a`, where each partition divides the data in half.  A worst
case analysis shows that the memory requirement can  increase to :math:`n \times
a`. In practise it usually requires 3 to 4 times :math:`a`.
Setting `min_density` to 0 will always use the sample mask to select the subset
of samples at each node.  This results in little to no additional memory being
allocated, making it appropriate for massive datasets or within ensemble
learners, but at the expense of being slower when training deep trees. The
default value for `min_density` is 0.1 which strikes a balance between
potentially excessive allocation of temporary memory and the long training
times.

Tips on practical use
=====================

  * Decision trees tend to overfit on data with a large number of features.
    Getting the right ratio of samples to number of features is important, since
    a tree with few samples in high dimensional space is very likely to overfit.

  * Consider performing  dimensionality reduction (:ref:`PCA <PCA>`,
    :ref:`ICA <ICA>`, or :ref:`feature_selection`) beforehand to
    give your tree a better chance of finding features that are discriminative.

  * Visualise your tree as you are training by using the ``export``
    function.  Use ``max_depth=3`` as an initial tree depth to get a feel for
    how the tree is fitting to your data, and then increase the depth.

  * Remember that the number of samples required to populate the tree doubles
    for each additional level the tree grows to.  Use ``max_depth`` to control
    the size of the tree to prevent overfitting.

  * Use ``min_split`` to control the number of samples at a leaf node.  A very
    small number will usually mean the tree will overfit, whereas a large number
    will prevent the tree from learning the data.  Try ``min_split=5`` as an
    initial value.

  * Balance your dataset before training to prevent the tree from creating
    a tree biased toward the classes that are dominant.

.. _tree_algorithms:

Tree algorithms: ID3, C4.5, C5.0 and CART
==========================================

What are all the various decision tree algorithms and how do they differ
from each other? Which one is implemented in scikit-learn?

ID3_ (Iterative Dichotomiser 3) was developed in 1986 by Ross Quinlan.
The algorithm creates a multiway tree, finding for each node (i.e. in
a greedy manner) the categorical feature that will yield the largest
information gain for categorical targets. Trees are grown to their
maximum size and then a pruning step is usually applied to improve the
ability of the tree to generalise to unseen data.

C4.5 is the successor to ID3 and removed the restriction that features
must be categorical by dynamically defining a discrete attribute (based
on numerical variables) that partitions the continuous attribute value
into a discrete set of intervals. C4.5 converts the trained trees
(i.e. the output of the ID3 algorithm) into sets of if-then rules.
These accuracy of each rule is then evaluated to determine the order
in which they should be applied. Pruning is done by removing a rule's
precondition if the accuracy of the rule improves without it.

C5.0 is Quinlan's latest version release under a proprietary license.
It uses less memory and builds smaller rulesets than C4.5 while being
more accurate.

CART_ (Classification and Regression Trees) is very similar to C4.5, but
it differs in that it supports numerical target variables (regression) and
does not compute rule sets. CART constructs binary trees using the feature
and threshold that yield the largest information gain at each node.

scikit-learn uses an optimised version of the CART algorithm.

.. _ID3: http://en.wikipedia.org/wiki/ID3_algorithm
.. _CART: http://en.wikipedia.org/wiki/Predictive_analytics#Classification_and_regression_trees

.. _tree_mathematical_formulation:

Mathematical formulation
========================

Given training vectors :math:`x_i \in R^n`, i=1,..., l and a label vector
:math:`y \in R^l`, a decision tree recursively partitions the space such
that the samples with the same labels are grouped together.

Let the data at node :math:`m` be represented by :math:`Q`. For
each candidate split :math:`\theta = (j, t_m)` consisting of a
feature :math:`j` and threshold :math:`t_m`, partition the data into
:math:`Q_{left}(\theta)` and :math:`Q_{right}(\theta)` subsets

.. math::

    Q_{left}(\theta) = {(x, y) | x_j <= t_m}

    Q_{right}(\theta) = Q \setminus Q_{left}(\theta)

The impurity at :math:`m` is computed using an impurity function
:math:`H()`, the choice of which depends on the task being solved
(classification or regression)

.. math::

   G(Q, \theta) = \frac{n_{left}}{N_m} H(Q_{left}(\theta))
   + \frac{n_{right}}{N_m} H(Q_{right}(\theta))

Select the parameters that minimises the impurity

.. math::

    \theta^* = argmin_\theta  G(Q, \theta)

Recurse for subsets :math:`Q_{left}(\theta^*)` and
:math:`Q_{right}(\theta^*)` until the maximum allowable depth is reached,
:math:`N_m < min\_samples` or :math:`N_m = 1`.

Classification criteria
-----------------------

If a target is a classification outcome taking on values 0,1,...,K-1,
for node :math:`m`, representing a region :math:`R_m` with :math:`N_m`
observations, let

.. math::

    p_{mk} = 1/ N_m \sum_{x_i \in R_m} I(y_i = k)

be the proportion of class k observations in node :math:`m`

Common measures of impurity are Gini

.. math::

    H(X_m) = \sum_k p_{mk} (1 - p_{mk})

Cross-Entropy

.. math::

    H(X_m) = \sum_k p_{mk} log(p_{mk})

and Misclassification

.. math::

    H(X_m) = 1 - max(p_{mk})

Regression criteria
-------------------

If the target is a continuous value, then for node :math:`m`,
representing a region :math:`R_m` with :math:`N_m` observations, a common
criterion to minimise is the Mean Squared Error

.. math::

    c_m = \frac{1}{N_m} \sum_{i \in N_m} y_i

    H(X_m) = \frac{1}{N_m} \sum_{i \in N_m} (y_i - c_m)^2


.. topic:: References:

    * http://en.wikipedia.org/wiki/Decision_tree_learning

    * http://en.wikipedia.org/wiki/Predictive_analytics

    * L. Breiman, J. Friedman, R. Olshen, and C. Stone. Classification and
      Regression Trees. Wadsworth, Belmont, CA, 1984.

    * J.R. Quinlan. C4. 5: programs for machine learning. Morgan Kaufmann, 1993.

    * T. Hastie, R. Tibshirani and J. Friedman.
      Elements of Statistical Learning, Springer, 2009.
