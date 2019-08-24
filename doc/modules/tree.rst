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

.. figure:: ../auto_examples/tree/images/sphx_glr_plot_tree_regression_001.png
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

    - Able to handle multi-output problems.

    - Uses a white box model. If a given situation is observable in a model,
      the explanation for the condition is easily explained by boolean logic.
      By contrast, in a black box model (e.g., in an artificial neural
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

As with other classifiers, :class:`DecisionTreeClassifier` takes as input two arrays:
an array X, sparse or dense, of size ``[n_samples, n_features]``  holding the
training samples, and an array Y of integer values, size ``[n_samples]``,
holding the class labels for the training samples::

    >>> from sklearn import tree
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(X, Y)

After being fitted, the model can then be used to predict the class of samples::

    >>> clf.predict([[2., 2.]])
    array([1])

Alternatively, the probability of each class can be predicted, which is the
fraction of training samples of the same class in a leaf::

    >>> clf.predict_proba([[2., 2.]])
    array([[0., 1.]])

:class:`DecisionTreeClassifier` is capable of both binary (where the
labels are [-1, 1]) classification and multiclass (where the labels are
[0, ..., K-1]) classification.

Using the Iris dataset, we can construct a tree as follows::

    >>> from sklearn.datasets import load_iris
    >>> from sklearn import tree
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(X, y)

Once trained, you can plot the tree with the plot_tree function::


    >>> tree.plot_tree(clf.fit(iris.data, iris.target)) # doctest: +SKIP

.. figure:: ../auto_examples/tree/images/sphx_glr_plot_iris_dtc_002.png
   :target: ../auto_examples/tree/plot_iris_dtc.html
   :scale: 75
   :align: center

We can also export the tree in `Graphviz
<https://www.graphviz.org/>`_ format using the :func:`export_graphviz`
exporter. If you use the `conda <https://conda.io>`_ package manager, the graphviz binaries  

and the python package can be installed with 

    conda install python-graphviz
   
Alternatively binaries for graphviz can be downloaded from the graphviz project homepage,
and the Python wrapper installed from pypi with `pip install graphviz`. 

Below is an example graphviz export of the above tree trained on the entire
iris dataset; the results are saved in an output file `iris.pdf`::


    >>> import graphviz # doctest: +SKIP
    >>> dot_data = tree.export_graphviz(clf, out_file=None) # doctest: +SKIP
    >>> graph = graphviz.Source(dot_data) # doctest: +SKIP
    >>> graph.render("iris") # doctest: +SKIP

The :func:`export_graphviz` exporter also supports a variety of aesthetic
options, including coloring nodes by their class (or value for regression) and
using explicit variable and class names if desired. Jupyter notebooks also
render these plots inline automatically::

    >>> dot_data = tree.export_graphviz(clf, out_file=None, # doctest: +SKIP
    ...                      feature_names=iris.feature_names,  # doctest: +SKIP
    ...                      class_names=iris.target_names,  # doctest: +SKIP
    ...                      filled=True, rounded=True,  # doctest: +SKIP
    ...                      special_characters=True)  # doctest: +SKIP
    >>> graph = graphviz.Source(dot_data)  # doctest: +SKIP
    >>> graph # doctest: +SKIP

.. only:: html

    .. figure:: ../images/iris.svg
       :align: center

.. only:: latex

    .. figure:: ../images/iris.pdf
       :align: center

.. figure:: ../auto_examples/tree/images/sphx_glr_plot_iris_dtc_001.png
   :target: ../auto_examples/tree/plot_iris_dtc.html
   :align: center
   :scale: 75

Alternatively, the tree can also be exported in textual format with the
function :func:`export_text`. This method doesn't require the installation
of external libraries and is more compact:

    >>> from sklearn.datasets import load_iris
    >>> from sklearn.tree import DecisionTreeClassifier
    >>> from sklearn.tree.export import export_text
    >>> iris = load_iris()
    >>> decision_tree = DecisionTreeClassifier(random_state=0, max_depth=2)
    >>> decision_tree = decision_tree.fit(iris.data, iris.target)
    >>> r = export_text(decision_tree, feature_names=iris['feature_names'])
    >>> print(r)
    |--- petal width (cm) <= 0.80
    |   |--- class: 0
    |--- petal width (cm) >  0.80
    |   |--- petal width (cm) <= 1.75
    |   |   |--- class: 1
    |   |--- petal width (cm) >  1.75
    |   |   |--- class: 2
    <BLANKLINE>

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_tree_plot_iris_dtc.py`
 * :ref:`sphx_glr_auto_examples_tree_plot_unveil_tree_structure.py`

.. _tree_regression:

Regression
==========

.. figure:: ../auto_examples/tree/images/sphx_glr_plot_tree_regression_001.png
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
    array([0.5])

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_tree_plot_tree_regression.py`


.. _tree_multioutput:

Multi-output problems
=====================

A multi-output problem is a supervised learning problem with several outputs
to predict, that is when Y is a 2d array of size ``[n_samples, n_outputs]``.

When there is no correlation between the outputs, a very simple way to solve
this kind of problem is to build n independent models, i.e. one for each
output, and then to use those models to independently predict each one of the n
outputs. However, because it is likely that the output values related to the
same input are themselves correlated, an often better way is to build a single
model capable of predicting simultaneously all n outputs. First, it requires
lower training time since only a single estimator is built. Second, the
generalization accuracy of the resulting estimator may often be increased.

With regard to decision trees, this strategy can readily be used to support
multi-output problems. This requires the following changes:

  - Store n output values in leaves, instead of 1;
  - Use splitting criteria that compute the average reduction across all
    n outputs.

This module offers support for multi-output problems by implementing this
strategy in both :class:`DecisionTreeClassifier` and
:class:`DecisionTreeRegressor`. If a decision tree is fit on an output array Y
of size ``[n_samples, n_outputs]`` then the resulting estimator will:

  * Output n_output values upon ``predict``;

  * Output a list of n_output arrays of class probabilities upon
    ``predict_proba``.


The use of multi-output trees for regression is demonstrated in
:ref:`sphx_glr_auto_examples_tree_plot_tree_regression_multioutput.py`. In this example, the input
X is a single real value and the outputs Y are the sine and cosine of X.

.. figure:: ../auto_examples/tree/images/sphx_glr_plot_tree_regression_multioutput_001.png
   :target: ../auto_examples/tree/plot_tree_regression_multioutput.html
   :scale: 75
   :align: center

The use of multi-output trees for classification is demonstrated in
:ref:`sphx_glr_auto_examples_plot_multioutput_face_completion.py`. In this example, the inputs
X are the pixels of the upper half of faces and the outputs Y are the pixels of
the lower half of those faces.

.. figure:: ../auto_examples/images/sphx_glr_plot_multioutput_face_completion_001.png
   :target: ../auto_examples/plot_multioutput_face_completion.html
   :scale: 75
   :align: center

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_tree_plot_tree_regression_multioutput.py`
 * :ref:`sphx_glr_auto_examples_plot_multioutput_face_completion.py`

.. topic:: References:

 * M. Dumont et al,  `Fast multi-class image annotation with random subwindows
   and multiple output randomized trees
   <http://www.montefiore.ulg.ac.be/services/stochastic/pubs/2009/DMWG09/dumont-visapp09-shortpaper.pdf>`_, International Conference on
   Computer Vision Theory and Applications 2009

.. _tree_complexity:

Complexity
==========

In general, the run time cost to construct a balanced binary tree is
:math:`O(n_{samples}n_{features}\log(n_{samples}))` and query time
:math:`O(\log(n_{samples}))`.  Although the tree construction algorithm attempts
to generate balanced trees, they will not always be balanced.  Assuming that the
subtrees remain approximately balanced, the cost at each node consists of
searching through :math:`O(n_{features})` to find the feature that offers the
largest reduction in entropy.  This has a cost of
:math:`O(n_{features}n_{samples}\log(n_{samples}))` at each node, leading to a
total cost over the entire trees (by summing the cost at each node) of
:math:`O(n_{features}n_{samples}^{2}\log(n_{samples}))`.

Scikit-learn offers a more efficient implementation for the construction of
decision trees.  A naive implementation (as above) would recompute the class
label histograms (for classification) or the means (for regression) at for each
new split point along a given feature. Presorting the feature over all
relevant samples, and retaining a running label count, will reduce the complexity
at each node to :math:`O(n_{features}\log(n_{samples}))`, which results in a
total cost of :math:`O(n_{features}n_{samples}\log(n_{samples}))`. This is an option
for all tree based algorithms. By default it is turned on for gradient boosting,
where in general it makes training faster, but turned off for all other algorithms as
it tends to slow down training when training deep trees.


Tips on practical use
=====================

  * Decision trees tend to overfit on data with a large number of features.
    Getting the right ratio of samples to number of features is important, since
    a tree with few samples in high dimensional space is very likely to overfit.

  * Consider performing  dimensionality reduction (:ref:`PCA <PCA>`,
    :ref:`ICA <ICA>`, or :ref:`feature_selection`) beforehand to
    give your tree a better chance of finding features that are discriminative.

  * :ref:`sphx_glr_auto_examples_tree_plot_unveil_tree_structure.py` will help
    in gaining more insights about how the decision tree makes predictions, which is
    important for understanding the important features in the data.

  * Visualise your tree as you are training by using the ``export``
    function.  Use ``max_depth=3`` as an initial tree depth to get a feel for
    how the tree is fitting to your data, and then increase the depth.

  * Remember that the number of samples required to populate the tree doubles
    for each additional level the tree grows to.  Use ``max_depth`` to control
    the size of the tree to prevent overfitting.

  * Use ``min_samples_split`` or ``min_samples_leaf`` to ensure that multiple
    samples inform every decision in the tree, by controlling which splits will
    be considered. A very small number will usually mean the tree will overfit,
    whereas a large number will prevent the tree from learning the data. Try
    ``min_samples_leaf=5`` as an initial value. If the sample size varies
    greatly, a float number can be used as percentage in these two parameters.
    While ``min_samples_split`` can create arbitrarily small leaves,
    ``min_samples_leaf`` guarantees that each leaf has a minimum size, avoiding
    low-variance, over-fit leaf nodes in regression problems.  For
    classification with few classes, ``min_samples_leaf=1`` is often the best
    choice.

  * Balance your dataset before training to prevent the tree from being biased
    toward the classes that are dominant. Class balancing can be done by
    sampling an equal number of samples from each class, or preferably by
    normalizing the sum of the sample weights (``sample_weight``) for each
    class to the same value. Also note that weight-based pre-pruning criteria,
    such as ``min_weight_fraction_leaf``, will then be less biased toward
    dominant classes than criteria that are not aware of the sample weights,
    like ``min_samples_leaf``.

  * If the samples are weighted, it will be easier to optimize the tree
    structure using weight-based pre-pruning criterion such as
    ``min_weight_fraction_leaf``, which ensure that leaf nodes contain at least
    a fraction of the overall sum of the sample weights.

  * All decision trees use ``np.float32`` arrays internally.
    If training data is not in this format, a copy of the dataset will be made.

  * If the input matrix X is very sparse, it is recommended to convert to sparse
    ``csc_matrix`` before calling fit and sparse ``csr_matrix`` before calling
    predict. Training time can be orders of magnitude faster for a sparse
    matrix input compared to a dense matrix when features have zero values in
    most of the samples.



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

scikit-learn uses an optimised version of the CART algorithm; however, scikit-learn 
implementation does not support categorical variables for now.

.. _ID3: https://en.wikipedia.org/wiki/ID3_algorithm
.. _CART: https://en.wikipedia.org/wiki/Predictive_analytics#Classification_and_regression_trees_.28CART.29


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

    \theta^* = \operatorname{argmin}_\theta  G(Q, \theta)

Recurse for subsets :math:`Q_{left}(\theta^*)` and
:math:`Q_{right}(\theta^*)` until the maximum allowable depth is reached,
:math:`N_m < \min_{samples}` or :math:`N_m = 1`.

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

Entropy

.. math::

    H(X_m) = - \sum_k p_{mk} \log(p_{mk})

and Misclassification

.. math::

    H(X_m) = 1 - \max(p_{mk})

where :math:`X_m` is the training data in node :math:`m`

Regression criteria
-------------------

If the target is a continuous value, then for node :math:`m`,
representing a region :math:`R_m` with :math:`N_m` observations, common
criteria to minimise as for determining locations for future
splits are Mean Squared Error, which minimizes the L2 error
using mean values at terminal nodes, and Mean Absolute Error, which 
minimizes the L1 error using median values at terminal nodes. 

Mean Squared Error:

.. math::

    \bar{y}_m = \frac{1}{N_m} \sum_{i \in N_m} y_i

    H(X_m) = \frac{1}{N_m} \sum_{i \in N_m} (y_i - \bar{y}_m)^2

Mean Absolute Error:

.. math::

    \bar{y}_m = \frac{1}{N_m} \sum_{i \in N_m} y_i

    H(X_m) = \frac{1}{N_m} \sum_{i \in N_m} |y_i - \bar{y}_m|

where :math:`X_m` is the training data in node :math:`m`


.. _minimal_cost_complexity_pruning:

Minimal Cost-Complexity Pruning
===============================

Minimal cost-complexity pruning is an algorithm used to prune a tree to avoid
over-fitting, described in Chapter 3 of [BRE]_. This algorithm is parameterized
by :math:`\alpha\ge0` known as the complexity parameter. The complexity
parameter is used to define the cost-complexity measure, :math:`R_\alpha(T)` of
a given tree :math:`T`:

.. math::

  R_\alpha(T) = R(T) + \alpha|T|

where :math:`|T|` is the number of terminal nodes in :math:`T` and :math:`R(T)`
is traditionally defined as the total misclassification rate of the terminal
nodes. Alternatively, scikit-learn uses the total sample weighted impurity of
the terminal nodes for :math:`R(T)`. As shown above, the impurity of a node
depends on the criterion. Minimal cost-complexity pruning finds the subtree of
:math:`T` that minimizes :math:`R_\alpha(T)`.

The cost complexity measure of a single node is
:math:`R_\alpha(t)=R(t)+\alpha`. The branch, :math:`T_t`, is defined to be a
tree where node :math:`t` is its root. In general, the impurity of a node
is greater than the sum of impurities of its terminal nodes,
:math:`R(T_t)<R(t)`. However, the cost complexity measure of a node,
:math:`t`, and its branch, :math:`T_t`, can be equal depending on
:math:`\alpha`. We define the effective :math:`\alpha` of a node to be the
value where they are equal, :math:`R_\alpha(T_t)=R_\alpha(t)` or
:math:`\alpha_{eff}(t)=\frac{R(t)-R(T_t)}{|T|-1}`. A non-terminal node
with the smallest value of :math:`\alpha_{eff}` is the weakest link and will
be pruned. This process stops when the pruned tree's minimal
:math:`\alpha_{eff}` is greater than the ``ccp_alpha`` parameter.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_tree_plot_cost_complexity_pruning.py`
  
.. topic:: References:

    .. [BRE] L. Breiman, J. Friedman, R. Olshen, and C. Stone. Classification
      and Regression Trees. Wadsworth, Belmont, CA, 1984.

    * https://en.wikipedia.org/wiki/Decision_tree_learning

    * https://en.wikipedia.org/wiki/Predictive_analytics

    * J.R. Quinlan. C4. 5: programs for machine learning. Morgan
      Kaufmann, 1993.

    * T. Hastie, R. Tibshirani and J. Friedman. Elements of Statistical
      Learning, Springer, 2009.
