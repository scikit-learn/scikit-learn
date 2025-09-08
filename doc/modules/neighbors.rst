.. _neighbors:

=================
Nearest Neighbors
=================

.. sectionauthor:: Jake Vanderplas <vanderplas@astro.washington.edu>

.. currentmodule:: sklearn.neighbors

:mod:`sklearn.neighbors` provides functionality for unsupervised and
supervised neighbors-based learning methods.  Unsupervised nearest neighbors
is the foundation of many other learning methods,
notably manifold learning and spectral clustering.  Supervised neighbors-based
learning comes in two flavors: `classification`_ for data with
discrete labels, and `regression`_ for data with continuous labels.

The principle behind nearest neighbor methods is to find a predefined number
of training samples closest in distance to the new point, and
predict the label from these.  The number of samples can be a user-defined
constant (k-nearest neighbor learning), or vary based
on the local density of points (radius-based neighbor learning).
The distance can, in general, be any metric measure: standard Euclidean
distance is the most common choice.
Neighbors-based methods are known as *non-generalizing* machine
learning methods, since they simply "remember" all of its training data
(possibly transformed into a fast indexing structure such as a
:ref:`Ball Tree <ball_tree>` or :ref:`KD Tree <kd_tree>`).

Despite its simplicity, nearest neighbors has been successful in a
large number of classification and regression problems, including
handwritten digits and satellite image scenes. Being a non-parametric method,
it is often successful in classification situations where the decision
boundary is very irregular.

The classes in :mod:`sklearn.neighbors` can handle either NumPy arrays or
`scipy.sparse` matrices as input.  For dense matrices, a large number of
possible distance metrics are supported.  For sparse matrices, arbitrary
Minkowski metrics are supported for searches.

There are many learning routines which rely on nearest neighbors at their
core.  One example is :ref:`kernel density estimation <kernel_density>`,
discussed in the :ref:`density estimation <density_estimation>` section.


.. _unsupervised_neighbors:

Unsupervised Nearest Neighbors
==============================

:class:`NearestNeighbors` implements unsupervised nearest neighbors learning.
It acts as a uniform interface to three different nearest neighbors
algorithms: :class:`BallTree`, :class:`KDTree`, and a
brute-force algorithm based on routines in :mod:`sklearn.metrics.pairwise`.
The choice of neighbors search algorithm is controlled through the keyword
``'algorithm'``, which must be one of
``['auto', 'ball_tree', 'kd_tree', 'brute']``.  When the default value
``'auto'`` is passed, the algorithm attempts to determine the best approach
from the training data.  For a discussion of the strengths and weaknesses
of each option, see `Nearest Neighbor Algorithms`_.

.. warning::

    Regarding the Nearest Neighbors algorithms, if two
    neighbors :math:`k+1` and :math:`k` have identical distances
    but different labels, the result will depend on the ordering of the
    training data.

Finding the Nearest Neighbors
-----------------------------
For the simple task of finding the nearest neighbors between two sets of
data, the unsupervised algorithms within :mod:`sklearn.neighbors` can be
used:

    >>> from sklearn.neighbors import NearestNeighbors
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
    >>> distances, indices = nbrs.kneighbors(X)
    >>> indices
    array([[0, 1],
           [1, 0],
           [2, 1],
           [3, 4],
           [4, 3],
           [5, 4]]...)
    >>> distances
    array([[0.        , 1.        ],
           [0.        , 1.        ],
           [0.        , 1.41421356],
           [0.        , 1.        ],
           [0.        , 1.        ],
           [0.        , 1.41421356]])

Because the query set matches the training set, the nearest neighbor of each
point is the point itself, at a distance of zero.

It is also possible to efficiently produce a sparse graph showing the
connections between neighboring points:

    >>> nbrs.kneighbors_graph(X).toarray()
    array([[1., 1., 0., 0., 0., 0.],
           [1., 1., 0., 0., 0., 0.],
           [0., 1., 1., 0., 0., 0.],
           [0., 0., 0., 1., 1., 0.],
           [0., 0., 0., 1., 1., 0.],
           [0., 0., 0., 0., 1., 1.]])

The dataset is structured such that points nearby in index order are nearby
in parameter space, leading to an approximately block-diagonal matrix of
K-nearest neighbors.  Such a sparse graph is useful in a variety of
circumstances which make use of spatial relationships between points for
unsupervised learning: in particular, see :class:`~sklearn.manifold.Isomap`,
:class:`~sklearn.manifold.LocallyLinearEmbedding`, and
:class:`~sklearn.cluster.SpectralClustering`.

KDTree and BallTree Classes
---------------------------
Alternatively, one can use the :class:`KDTree` or :class:`BallTree` classes
directly to find nearest neighbors.  This is the functionality wrapped by
the :class:`NearestNeighbors` class used above.  The Ball Tree and KD Tree
have the same interface; we'll show an example of using the KD Tree here:

    >>> from sklearn.neighbors import KDTree
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> kdt = KDTree(X, leaf_size=30, metric='euclidean')
    >>> kdt.query(X, k=2, return_distance=False)
    array([[0, 1],
           [1, 0],
           [2, 1],
           [3, 4],
           [4, 3],
           [5, 4]]...)

Refer to the :class:`KDTree` and :class:`BallTree` class documentation
for more information on the options available for nearest neighbors searches,
including specification of query strategies, distance metrics, etc. For a list
of valid metrics use `KDTree.valid_metrics` and `BallTree.valid_metrics`:

    >>> from sklearn.neighbors import KDTree, BallTree
    >>> KDTree.valid_metrics
    ['euclidean', 'l2', 'minkowski', 'p', 'manhattan', 'cityblock', 'l1', 'chebyshev', 'infinity']
    >>> BallTree.valid_metrics
    ['euclidean', 'l2', 'minkowski', 'p', 'manhattan', 'cityblock', 'l1', 'chebyshev', 'infinity', 'seuclidean', 'mahalanobis', 'hamming', 'canberra', 'braycurtis', 'jaccard', 'dice', 'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath', 'haversine', 'pyfunc']

.. _classification:

Nearest Neighbors Classification
================================

Neighbors-based classification is a type of *instance-based learning* or
*non-generalizing learning*: it does not attempt to construct a general
internal model, but simply stores instances of the training data.
Classification is computed from a simple majority vote of the nearest
neighbors of each point: a query point is assigned the data class which
has the most representatives within the nearest neighbors of the point.

scikit-learn implements two different nearest neighbors classifiers:
:class:`KNeighborsClassifier` implements learning based on the :math:`k`
nearest neighbors of each query point, where :math:`k` is an integer value
specified by the user.  :class:`RadiusNeighborsClassifier` implements learning
based on the number of neighbors within a fixed radius :math:`r` of each
training point, where :math:`r` is a floating-point value specified by
the user.

The :math:`k`-neighbors classification in :class:`KNeighborsClassifier`
is the most commonly used technique. The optimal choice of the value :math:`k`
is highly data-dependent: in general a larger :math:`k` suppresses the effects
of noise, but makes the classification boundaries less distinct.

In cases where the data is not uniformly sampled, radius-based neighbors
classification in :class:`RadiusNeighborsClassifier` can be a better choice.
The user specifies a fixed radius :math:`r`, such that points in sparser
neighborhoods use fewer nearest neighbors for the classification.  For
high-dimensional parameter spaces, this method becomes less effective due
to the so-called "curse of dimensionality".

The basic nearest neighbors classification uses uniform weights: that is, the
value assigned to a query point is computed from a simple majority vote of
the nearest neighbors.  Under some circumstances, it is better to weight the
neighbors such that nearer neighbors contribute more to the fit.  This can
be accomplished through the ``weights`` keyword.  The default value,
``weights = 'uniform'``, assigns uniform weights to each neighbor.
``weights = 'distance'`` assigns weights proportional to the inverse of the
distance from the query point.  Alternatively, a user-defined function of the
distance can be supplied to compute the weights.

.. |classification_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_classification_001.png
   :target: ../auto_examples/neighbors/plot_classification.html
   :scale: 75

.. centered:: |classification_1|

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_neighbors_plot_classification.py`: an example of
  classification using nearest neighbors.

.. _regression:

Nearest Neighbors Regression
============================

Neighbors-based regression can be used in cases where the data labels are
continuous rather than discrete variables.  The label assigned to a query
point is computed based on the mean of the labels of its nearest neighbors.

scikit-learn implements two different neighbors regressors:
:class:`KNeighborsRegressor` implements learning based on the :math:`k`
nearest neighbors of each query point, where :math:`k` is an integer
value specified by the user.  :class:`RadiusNeighborsRegressor` implements
learning based on the neighbors within a fixed radius :math:`r` of the
query point, where :math:`r` is a floating-point value specified by the
user.

The basic nearest neighbors regression uses uniform weights: that is,
each point in the local neighborhood contributes uniformly to the
classification of a query point.  Under some circumstances, it can be
advantageous to weight points such that nearby points contribute more
to the regression than faraway points.  This can be accomplished through
the ``weights`` keyword.  The default value, ``weights = 'uniform'``,
assigns equal weights to all points.  ``weights = 'distance'`` assigns
weights proportional to the inverse of the distance from the query point.
Alternatively, a user-defined function of the distance can be supplied,
which will be used to compute the weights.

.. figure:: ../auto_examples/neighbors/images/sphx_glr_plot_regression_001.png
   :target: ../auto_examples/neighbors/plot_regression.html
   :align: center
   :scale: 75

The use of multi-output nearest neighbors for regression is demonstrated in
:ref:`sphx_glr_auto_examples_miscellaneous_plot_multioutput_face_completion.py`. In this example, the inputs
X are the pixels of the upper half of faces and the outputs Y are the pixels of
the lower half of those faces.

.. figure:: ../auto_examples/miscellaneous/images/sphx_glr_plot_multioutput_face_completion_001.png
   :target: ../auto_examples/miscellaneous/plot_multioutput_face_completion.html
   :scale: 75
   :align: center


.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_neighbors_plot_regression.py`: an example of regression
  using nearest neighbors.

* :ref:`sphx_glr_auto_examples_miscellaneous_plot_multioutput_face_completion.py`:
  an example of multi-output regression using nearest neighbors.


Nearest Neighbor Algorithms
===========================

.. _brute_force:

Brute Force
-----------

Fast computation of nearest neighbors is an active area of research in
machine learning. The most naive neighbor search implementation involves
the brute-force computation of distances between all pairs of points in the
dataset: for :math:`N` samples in :math:`D` dimensions, this approach scales
as :math:`O[D N^2]`.  Efficient brute-force neighbors searches can be very
competitive for small data samples.
However, as the number of samples :math:`N` grows, the brute-force
approach quickly becomes infeasible.  In the classes within
:mod:`sklearn.neighbors`, brute-force neighbors searches are specified
using the keyword ``algorithm = 'brute'``, and are computed using the
routines available in :mod:`sklearn.metrics.pairwise`.

.. _kd_tree:

K-D Tree
--------

To address the computational inefficiencies of the brute-force approach, a
variety of tree-based data structures have been invented.  In general, these
structures attempt to reduce the required number of distance calculations
by efficiently encoding aggregate distance information for the sample.
The basic idea is that if point :math:`A` is very distant from point
:math:`B`, and point :math:`B` is very close to point :math:`C`,
then we know that points :math:`A` and :math:`C`
are very distant, *without having to explicitly calculate their distance*.
In this way, the computational cost of a nearest neighbors search can be
reduced to :math:`O[D N \log(N)]` or better. This is a significant
improvement over brute-force for large :math:`N`.

An early approach to taking advantage of this aggregate information was
the *KD tree* data structure (short for *K-dimensional tree*), which
generalizes two-dimensional *Quad-trees* and 3-dimensional *Oct-trees*
to an arbitrary number of dimensions.  The KD tree is a binary tree
structure which recursively partitions the parameter space along the data
axes, dividing it into nested orthotropic regions into which data points
are filed.  The construction of a KD tree is very fast: because partitioning
is performed only along the data axes, no :math:`D`-dimensional distances
need to be computed. Once constructed, the nearest neighbor of a query
point can be determined with only :math:`O[\log(N)]` distance computations.
Though the KD tree approach is very fast for low-dimensional (:math:`D < 20`)
neighbors searches, it becomes inefficient as :math:`D` grows very large:
this is one manifestation of the so-called "curse of dimensionality".
In scikit-learn, KD tree neighbors searches are specified using the
keyword ``algorithm = 'kd_tree'``, and are computed using the class
:class:`KDTree`.


.. dropdown:: References

  * `"Multidimensional binary search trees used for associative searching"
    <https://dl.acm.org/citation.cfm?doid=361002.361007>`_,
    Bentley, J.L., Communications of the ACM (1975)


.. _ball_tree:

Ball Tree
---------

To address the inefficiencies of KD Trees in higher dimensions, the *ball tree*
data structure was developed.  Where KD trees partition data along
Cartesian axes, ball trees partition data in a series of nesting
hyper-spheres.  This makes tree construction more costly than that of the
KD tree, but results in a data structure which can be very efficient on
highly structured data, even in very high dimensions.

A ball tree recursively divides the data into
nodes defined by a centroid :math:`C` and radius :math:`r`, such that each
point in the node lies within the hyper-sphere defined by :math:`r` and
:math:`C`. The number of candidate points for a neighbor search
is reduced through use of the *triangle inequality*:

.. math::   |x+y| \leq |x| + |y|

With this setup, a single distance calculation between a test point and
the centroid is sufficient to determine a lower and upper bound on the
distance to all points within the node.
Because of the spherical geometry of the ball tree nodes, it can out-perform
a *KD-tree* in high dimensions, though the actual performance is highly
dependent on the structure of the training data.
In scikit-learn, ball-tree-based
neighbors searches are specified using the keyword ``algorithm = 'ball_tree'``,
and are computed using the class :class:`BallTree`.
Alternatively, the user can work with the :class:`BallTree` class directly.


.. dropdown:: References

  * `"Five Balltree Construction Algorithms"
    <https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=17ac002939f8e950ffb32ec4dc8e86bdd8cb5ff1>`_,
    Omohundro, S.M., International Computer Science Institute
    Technical Report (1989)

.. dropdown:: Choice of Nearest Neighbors Algorithm

  The optimal algorithm for a given dataset is a complicated choice, and
  depends on a number of factors:

  * number of samples :math:`N` (i.e. ``n_samples``) and dimensionality
    :math:`D` (i.e. ``n_features``).

    * *Brute force* query time grows as :math:`O[D N]`
    * *Ball tree* query time grows as approximately :math:`O[D \log(N)]`
    * *KD tree* query time changes with :math:`D` in a way that is difficult
      to precisely characterise.  For small :math:`D` (less than 20 or so)
      the cost is approximately :math:`O[D\log(N)]`, and the KD tree
      query can be very efficient.
      For larger :math:`D`, the cost increases to nearly :math:`O[DN]`, and
      the overhead due to the tree
      structure can lead to queries which are slower than brute force.

    For small data sets (:math:`N` less than 30 or so), :math:`\log(N)` is
    comparable to :math:`N`, and brute force algorithms can be more efficient
    than a tree-based approach.  Both :class:`KDTree` and :class:`BallTree`
    address this through providing a *leaf size* parameter: this controls the
    number of samples at which a query switches to brute-force.  This allows both
    algorithms to approach the efficiency of a brute-force computation for small
    :math:`N`.

  * data structure: *intrinsic dimensionality* of the data and/or *sparsity*
    of the data. Intrinsic dimensionality refers to the dimension
    :math:`d \le D` of a manifold on which the data lies, which can be linearly
    or non-linearly embedded in the parameter space. Sparsity refers to the
    degree to which the data fills the parameter space (this is to be
    distinguished from the concept as used in "sparse" matrices.  The data
    matrix may have no zero entries, but the **structure** can still be
    "sparse" in this sense).

    * *Brute force* query time is unchanged by data structure.
    * *Ball tree* and *KD tree* query times can be greatly influenced
      by data structure.  In general, sparser data with a smaller intrinsic
      dimensionality leads to faster query times.  Because the KD tree
      internal representation is aligned with the parameter axes, it will not
      generally show as much improvement as ball tree for arbitrarily
      structured data.

    Datasets used in machine learning tend to be very structured, and are
    very well-suited for tree-based queries.

  * number of neighbors :math:`k` requested for a query point.

    * *Brute force* query time is largely unaffected by the value of :math:`k`
    * *Ball tree* and *KD tree* query time will become slower as :math:`k`
      increases.  This is due to two effects: first, a larger :math:`k` leads
      to the necessity to search a larger portion of the parameter space.
      Second, using :math:`k > 1` requires internal queueing of results
      as the tree is traversed.

    As :math:`k` becomes large compared to :math:`N`, the ability to prune
    branches in a tree-based query is reduced.  In this situation, Brute force
    queries can be more efficient.

  * number of query points.  Both the ball tree and the KD Tree
    require a construction phase.  The cost of this construction becomes
    negligible when amortized over many queries.  If only a small number of
    queries will be performed, however, the construction can make up
    a significant fraction of the total cost.  If very few query points
    will be required, brute force is better than a tree-based method.

  Currently, ``algorithm = 'auto'`` selects ``'brute'`` if any of the following
  conditions are verified:

  * input data is sparse
  * ``metric = 'precomputed'``
  * :math:`D > 15`
  * :math:`k >= N/2`
  * ``effective_metric_`` isn't in the ``VALID_METRICS`` list for either
    ``'kd_tree'`` or ``'ball_tree'``

  Otherwise, it selects the first out of ``'kd_tree'`` and ``'ball_tree'`` that
  has ``effective_metric_`` in its ``VALID_METRICS`` list. This heuristic is
  based on the following assumptions:

  * the number of query points is at least the same order as the number of
    training points
  * ``leaf_size`` is close to its default value of ``30``
  * when :math:`D > 15`, the intrinsic dimensionality of the data is generally
    too high for tree-based methods

.. dropdown:: Effect of ``leaf_size``

  As noted above, for small sample sizes a brute force search can be more
  efficient than a tree-based query.  This fact is accounted for in the ball
  tree and KD tree by internally switching to brute force searches within
  leaf nodes.  The level of this switch can be specified with the parameter
  ``leaf_size``.  This parameter choice has many effects:

  **construction time**
    A larger ``leaf_size`` leads to a faster tree construction time, because
    fewer nodes need to be created

  **query time**
    Both a large or small ``leaf_size`` can lead to suboptimal query cost.
    For ``leaf_size`` approaching 1, the overhead involved in traversing
    nodes can significantly slow query times.  For ``leaf_size`` approaching
    the size of the training set, queries become essentially brute force.
    A good compromise between these is ``leaf_size = 30``, the default value
    of the parameter.

  **memory**
    As ``leaf_size`` increases, the memory required to store a tree structure
    decreases.  This is especially important in the case of ball tree, which
    stores a :math:`D`-dimensional centroid for each node.  The required
    storage space for :class:`BallTree` is approximately ``1 / leaf_size`` times
    the size of the training set.

  ``leaf_size`` is not referenced for brute force queries.

.. dropdown:: Valid Metrics for Nearest Neighbor Algorithms

  For a list of available metrics, see the documentation of the
  :class:`~sklearn.metrics.DistanceMetric` class and the metrics listed in
  `sklearn.metrics.pairwise.PAIRWISE_DISTANCE_FUNCTIONS`. Note that the "cosine"
  metric uses :func:`~sklearn.metrics.pairwise.cosine_distances`.

  A list of valid metrics for any of the above algorithms can be obtained by using their
  ``valid_metric`` attribute. For example, valid metrics for ``KDTree`` can be generated by:

      >>> from sklearn.neighbors import KDTree
      >>> print(sorted(KDTree.valid_metrics))
      ['chebyshev', 'cityblock', 'euclidean', 'infinity', 'l1', 'l2', 'manhattan', 'minkowski', 'p']

.. _nearest_centroid_classifier:

Nearest Centroid Classifier
===========================

The :class:`NearestCentroid` classifier is a simple algorithm that represents
each class by the centroid of its members. In effect, this makes it
similar to the label updating phase of the :class:`~sklearn.cluster.KMeans` algorithm.
It also has no parameters to choose, making it a good baseline classifier. It
does, however, suffer on non-convex classes, as well as when classes have
drastically different variances, as equal variance in all dimensions is
assumed. See Linear Discriminant Analysis (:class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis`)
and Quadratic Discriminant Analysis (:class:`~sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis`)
for more complex methods that do not make this assumption. Usage of the default
:class:`NearestCentroid` is simple:

    >>> from sklearn.neighbors import NearestCentroid
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = NearestCentroid()
    >>> clf.fit(X, y)
    NearestCentroid()
    >>> print(clf.predict([[-0.8, -1]]))
    [1]


Nearest Shrunken Centroid
-------------------------

The :class:`NearestCentroid` classifier has a ``shrink_threshold`` parameter,
which implements the nearest shrunken centroid classifier. In effect, the value
of each feature for each centroid is divided by the within-class variance of
that feature. The feature values are then reduced by ``shrink_threshold``. Most
notably, if a particular feature value crosses zero, it is set
to zero. In effect, this removes the feature from affecting the classification.
This is useful, for example, for removing noisy features.

In the example below, using a small shrink threshold increases the accuracy of
the model from 0.81 to 0.82.

.. |nearest_centroid_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nearest_centroid_001.png
   :target: ../auto_examples/neighbors/plot_nearest_centroid.html
   :scale: 50

.. |nearest_centroid_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nearest_centroid_002.png
   :target: ../auto_examples/neighbors/plot_nearest_centroid.html
   :scale: 50

.. centered:: |nearest_centroid_1| |nearest_centroid_2|

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_neighbors_plot_nearest_centroid.py`: an example of
  classification using nearest centroid with different shrink thresholds.

.. _neighbors_transformer:

Nearest Neighbors Transformer
=============================

Many scikit-learn estimators rely on nearest neighbors: Several classifiers and
regressors such as :class:`KNeighborsClassifier` and
:class:`KNeighborsRegressor`, but also some clustering methods such as
:class:`~sklearn.cluster.DBSCAN` and
:class:`~sklearn.cluster.SpectralClustering`, and some manifold embeddings such
as :class:`~sklearn.manifold.TSNE` and :class:`~sklearn.manifold.Isomap`.

All these estimators can compute internally the nearest neighbors, but most of
them also accept precomputed nearest neighbors :term:`sparse graph`,
as given by :func:`~sklearn.neighbors.kneighbors_graph` and
:func:`~sklearn.neighbors.radius_neighbors_graph`. With mode
`mode='connectivity'`, these functions return a binary adjacency sparse graph
as required, for instance, in :class:`~sklearn.cluster.SpectralClustering`.
Whereas with `mode='distance'`, they return a distance sparse graph as required,
for instance, in :class:`~sklearn.cluster.DBSCAN`. To include these functions in
a scikit-learn pipeline, one can also use the corresponding classes
:class:`KNeighborsTransformer` and :class:`RadiusNeighborsTransformer`.
The benefits of this sparse graph API are multiple.

First, the precomputed graph can be reused multiple times, for instance while
varying a parameter of the estimator. This can be done manually by the user, or
using the caching properties of the scikit-learn pipeline:

    >>> import tempfile
    >>> from sklearn.manifold import Isomap
    >>> from sklearn.neighbors import KNeighborsTransformer
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.datasets import make_regression
    >>> cache_path = tempfile.gettempdir()  # we use a temporary folder here
    >>> X, _ = make_regression(n_samples=50, n_features=25, random_state=0)
    >>> estimator = make_pipeline(
    ...     KNeighborsTransformer(mode='distance'),
    ...     Isomap(n_components=3, metric='precomputed'),
    ...     memory=cache_path)
    >>> X_embedded = estimator.fit_transform(X)
    >>> X_embedded.shape
    (50, 3)

Second, precomputing the graph can give finer control on the nearest neighbors
estimation, for instance enabling multiprocessing though the parameter
`n_jobs`, which might not be available in all estimators.

Finally, the precomputation can be performed by custom estimators to use
different implementations, such as approximate nearest neighbors methods, or
implementation with special data types. The precomputed neighbors
:term:`sparse graph` needs to be formatted as in
:func:`~sklearn.neighbors.radius_neighbors_graph` output:

* a CSR matrix (although COO, CSC or LIL will be accepted).
* only explicitly store nearest neighborhoods of each sample with respect to the
  training data. This should include those at 0 distance from a query point,
  including the matrix diagonal when computing the nearest neighborhoods
  between the training data and itself.
* each row's `data` should store the distance in increasing order (optional.
  Unsorted data will be stable-sorted, adding a computational overhead).
* all values in data should be non-negative.
* there should be no duplicate `indices` in any row
  (see https://github.com/scipy/scipy/issues/5807).
* if the algorithm being passed the precomputed matrix uses k nearest neighbors
  (as opposed to radius neighborhood), at least k neighbors must be stored in
  each row (or k+1, as explained in the following note).

.. note::
  When a specific number of neighbors is queried (using
  :class:`KNeighborsTransformer`), the definition of `n_neighbors` is ambiguous
  since it can either include each training point as its own neighbor, or
  exclude them. Neither choice is perfect, since including them leads to a
  different number of non-self neighbors during training and testing, while
  excluding them leads to a difference between `fit(X).transform(X)` and
  `fit_transform(X)`, which is against scikit-learn API.
  In :class:`KNeighborsTransformer` we use the definition which includes each
  training point as its own neighbor in the count of `n_neighbors`. However,
  for compatibility reasons with other estimators which use the other
  definition, one extra neighbor will be computed when `mode == 'distance'`.
  To maximise compatibility with all estimators, a safe choice is to always
  include one extra neighbor in a custom nearest neighbors estimator, since
  unnecessary neighbors will be filtered by following estimators.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_neighbors_approximate_nearest_neighbors.py`:
  an example of pipelining :class:`KNeighborsTransformer` and
  :class:`~sklearn.manifold.TSNE`. Also proposes two custom nearest neighbors
  estimators based on external packages.

* :ref:`sphx_glr_auto_examples_neighbors_plot_caching_nearest_neighbors.py`:
  an example of pipelining :class:`KNeighborsTransformer` and
  :class:`KNeighborsClassifier` to enable caching of the neighbors graph
  during a hyper-parameter grid-search.

.. _nca:

Neighborhood Components Analysis
================================

.. sectionauthor:: William de Vazelhes <william.de-vazelhes@inria.fr>

Neighborhood Components Analysis (NCA, :class:`NeighborhoodComponentsAnalysis`)
is a distance metric learning algorithm which aims to improve the accuracy of
nearest neighbors classification compared to the standard Euclidean distance.
The algorithm directly maximizes a stochastic variant of the leave-one-out
k-nearest neighbors (KNN) score on the training set. It can also learn a
low-dimensional linear projection of data that can be used for data
visualization and fast classification.

.. |nca_illustration_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_illustration_001.png
   :target: ../auto_examples/neighbors/plot_nca_illustration.html
   :scale: 50

.. |nca_illustration_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_illustration_002.png
   :target: ../auto_examples/neighbors/plot_nca_illustration.html
   :scale: 50

.. centered:: |nca_illustration_1| |nca_illustration_2|

In the above illustrating figure, we consider some points from a randomly
generated dataset. We focus on the stochastic KNN classification of point no.
3. The thickness of a link between sample 3 and another point is proportional
to their distance, and can be seen as the relative weight (or probability) that
a stochastic nearest neighbor prediction rule would assign to this point. In
the original space, sample 3 has many stochastic neighbors from various
classes, so the right class is not very likely. However, in the projected space
learned by NCA, the only stochastic neighbors with non-negligible weight are
from the same class as sample 3, guaranteeing that the latter will be well
classified. See the :ref:`mathematical formulation <nca_mathematical_formulation>`
for more details.


Classification
--------------

Combined with a nearest neighbors classifier (:class:`KNeighborsClassifier`),
NCA is attractive for classification because it can naturally handle
multi-class problems without any increase in the model size, and does not
introduce additional parameters that require fine-tuning by the user.

NCA classification has been shown to work well in practice for data sets of
varying size and difficulty. In contrast to related methods such as Linear
Discriminant Analysis, NCA does not make any assumptions about the class
distributions. The nearest neighbor classification can naturally produce highly
irregular decision boundaries.

To use this model for classification, one needs to combine a
:class:`NeighborhoodComponentsAnalysis` instance that learns the optimal
transformation with a :class:`KNeighborsClassifier` instance that performs the
classification in the projected space. Here is an example using the two
classes:

    >>> from sklearn.neighbors import (NeighborhoodComponentsAnalysis,
    ... KNeighborsClassifier)
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.pipeline import Pipeline
    >>> X, y = load_iris(return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ... stratify=y, test_size=0.7, random_state=42)
    >>> nca = NeighborhoodComponentsAnalysis(random_state=42)
    >>> knn = KNeighborsClassifier(n_neighbors=3)
    >>> nca_pipe = Pipeline([('nca', nca), ('knn', knn)])
    >>> nca_pipe.fit(X_train, y_train)
    Pipeline(...)
    >>> print(nca_pipe.score(X_test, y_test))
    0.96190476...

.. |nca_classification_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_classification_001.png
   :target: ../auto_examples/neighbors/plot_nca_classification.html
   :scale: 50

.. |nca_classification_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_classification_002.png
   :target: ../auto_examples/neighbors/plot_nca_classification.html
   :scale: 50

.. centered:: |nca_classification_1| |nca_classification_2|

The plot shows decision boundaries for Nearest Neighbor Classification and
Neighborhood Components Analysis classification on the iris dataset, when
training and scoring on only two features, for visualisation purposes.

.. _nca_dim_reduction:

Dimensionality reduction
------------------------

NCA can be used to perform supervised dimensionality reduction. The input data
are projected onto a linear subspace consisting of the directions which
minimize the NCA objective. The desired dimensionality can be set using the
parameter ``n_components``. For instance, the following figure shows a
comparison of dimensionality reduction with Principal Component Analysis
(:class:`~sklearn.decomposition.PCA`), Linear Discriminant Analysis
(:class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis`) and
Neighborhood Component Analysis (:class:`NeighborhoodComponentsAnalysis`) on
the Digits dataset, a dataset with size :math:`n_{samples} = 1797` and
:math:`n_{features} = 64`. The data set is split into a training and a test set
of equal size, then standardized. For evaluation the 3-nearest neighbor
classification accuracy is computed on the 2-dimensional projected points found
by each method. Each data sample belongs to one of 10 classes.

.. |nca_dim_reduction_1| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_dim_reduction_001.png
   :target: ../auto_examples/neighbors/plot_nca_dim_reduction.html
   :width: 32%

.. |nca_dim_reduction_2| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_dim_reduction_002.png
   :target: ../auto_examples/neighbors/plot_nca_dim_reduction.html
   :width: 32%

.. |nca_dim_reduction_3| image:: ../auto_examples/neighbors/images/sphx_glr_plot_nca_dim_reduction_003.png
   :target: ../auto_examples/neighbors/plot_nca_dim_reduction.html
   :width: 32%

.. centered:: |nca_dim_reduction_1| |nca_dim_reduction_2| |nca_dim_reduction_3|


.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_neighbors_plot_nca_classification.py`
* :ref:`sphx_glr_auto_examples_neighbors_plot_nca_dim_reduction.py`
* :ref:`sphx_glr_auto_examples_manifold_plot_lle_digits.py`

.. _nca_mathematical_formulation:

Mathematical formulation
------------------------

The goal of NCA is to learn an optimal linear transformation matrix of size
``(n_components, n_features)``, which maximises the sum over all samples
:math:`i` of the probability :math:`p_i` that :math:`i` is correctly
classified, i.e.:

.. math::

  \underset{L}{\arg\max} \sum\limits_{i=0}^{N - 1} p_{i}

with :math:`N` = ``n_samples`` and :math:`p_i` the probability of sample
:math:`i` being correctly classified according to a stochastic nearest
neighbors rule in the learned embedded space:

.. math::

  p_{i}=\sum\limits_{j \in C_i}{p_{i j}}

where :math:`C_i` is the set of points in the same class as sample :math:`i`,
and :math:`p_{i j}` is the softmax over Euclidean distances in the embedded
space:

.. math::

  p_{i j} = \frac{\exp(-||L x_i - L x_j||^2)}{\sum\limits_{k \ne
            i} {\exp{-(||L x_i - L x_k||^2)}}} , \quad p_{i i} = 0

.. dropdown:: Mahalanobis distance

  NCA can be seen as learning a (squared) Mahalanobis distance metric:

  .. math::

      || L(x_i - x_j)||^2 = (x_i - x_j)^TM(x_i - x_j),

  where :math:`M = L^T L` is a symmetric positive semi-definite matrix of size
  ``(n_features, n_features)``.


Implementation
--------------

This implementation follows what is explained in the original paper [1]_. For
the optimisation method, it currently uses scipy's L-BFGS-B with a full
gradient computation at each iteration, to avoid to tune the learning rate and
provide stable learning.

See the examples below and the docstring of
:meth:`NeighborhoodComponentsAnalysis.fit` for further information.

Complexity
----------

Training
^^^^^^^^
NCA stores a matrix of pairwise distances, taking ``n_samples ** 2`` memory.
Time complexity depends on the number of iterations done by the optimisation
algorithm. However, one can set the maximum number of iterations with the
argument ``max_iter``. For each iteration, time complexity is
``O(n_components x n_samples x min(n_samples, n_features))``.


Transform
^^^^^^^^^
Here the ``transform`` operation returns :math:`LX^T`, therefore its time
complexity equals ``n_components * n_features * n_samples_test``. There is no
added space complexity in the operation.


.. rubric:: References

.. [1] `"Neighbourhood Components Analysis"
  <https://papers.nips.cc/paper_files/paper/2004/file/42fe880812925e520249e808937738d2-Paper.pdf>`_,
  J. Goldberger, S. Roweis, G. Hinton, R. Salakhutdinov, Advances in
  Neural Information Processing Systems, Vol. 17, May 2005, pp. 513-520.

* `Wikipedia entry on Neighborhood Components Analysis
  <https://en.wikipedia.org/wiki/Neighbourhood_components_analysis>`_
