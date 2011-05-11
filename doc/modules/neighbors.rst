=================
Nearest Neighbors
=================

.. currentmodule:: scikits.learn.neighbors

The principle behind nearest neighbor methods is to find the k
training points closest in euclidean distance to x0, and then classify
using a majority vote among the k neighbors.

Despite its simplicity, nearest neighbors has been successful in a
large number of classification problems, including handwritten digits
or satellite image scenes. It is ofter successful in situation where
the decision boundary is very irregular.

Classification
==============

The :class:`NeighborsClassifier` implements the nearest-neighbors
classification method using a vote heuristic: the class most present
in the k nearest neighbors of a point is assigned to this point.

It is possible to use different nearest neighbor search algorithms by
using the keyword ``algorithm``. Possible values are ``'auto'``,
``'ball_tree'``, ``'brute'`` and ``'brute_inplace'``. ``'ball_tree'``
will create an instance of :class:`BallTree` to conduct the search,
which is usually very efficient in low-dimensional spaces. In higher
dimension, a brute-force approach is prefered thus parameters
``'brute'`` and ``'brute_inplace'`` can be used . Both conduct a
brute-force search, the difference being that ``'brute_inplace'`` does
not perform any precomputations, and thus is better suited for
low-memory settings.  Finally, ``'auto'`` is a simple heuristic that
will guess the best approach based on the current dataset.


.. figure:: ../auto_examples/images/plot_neighbors_1.png
   :target: ../auto_examples/plot_neighbors.html
   :align: center
   :scale: 75


.. topic:: Examples:

  * :ref:`example_plot_neighbors.py`: an example of classification
    using nearest neighbor.


Regression
==========

The :class:`NeighborsRegressor` estimator implements a
nearest-neighbors regression method by weighting the targets of the
k-neighbors. Two different weighting strategies are implemented:
``barycenter`` and ``mean``. ``barycenter`` will apply the weights
that best reconstruct the point from its neighbors while ``mean`` will
apply constant weights to each point. This plot shows the behavior of
both classifier for a simple regression task.

.. figure:: ../auto_examples/images/plot_neighbors_regression_1.png
   :target: ../auto_examples/plot_neighbors_regression.html
   :align: center
   :scale: 75


.. topic:: Examples:

  * :ref:`example_plot_neighbors_regression.py`: an example of regression
    using nearest neighbor.

.. _ball_tree:

Efficient implementation: the ball tree
==========================================

Behind the scenes, nearest neighbor search is done by the object
:class:`BallTree`. This algorithm makes it possible to rapidly find
the nearest neighbors of a point in N-dimensional spaces.

The optimal neighbor search algorithm for any problem
depends on a few factors: the number of candidate points (N), 
the dimensionality of the parameter space (D),
the structure of the data, and the number of neighbors desired (K).
For small N (less than 500 or so), a brute force search is 
generally adequate.  As N increases, however, the computation time 
for a neighbor search grows as O[N].  Various tree methods have
been developed to reduce this to O[log(N)], the simplest of which is the
KD tree.  A KD tree gains this speedup by recursively partitioning the data,
each time along a single dimension, so that clusters of candidate
neighbors can be ruled-out via a single distance calculation.  Unfortunately,
as D grows large compared to log(N), KD trees become very inefficient:
this is one manifestation of the so-called "curse of dimensionality".

Ball trees were developed as an alternative for high-dimensional neighbor
searches.  Compared to a KDTree, the cost is a slightly longer construction 
time, though for repeated queries, this added construction time quickly
becomes insignificant. The advantage of the Ball Tree is that it is 
efficient even as the dimensionality D becomes large.

A ball tree recursively divides the data into
nodes defined by a centroid C and radius R, such that each
point in the node lies within the hyper-sphere of radius R centered
at C.  The number of candidate points for a neighbor search 
is reduced through use of the triangle inequality:

.. math::   |x+y| \leq |x| + |y|

With this setup, a single distance calculation between a test point and
the centroid is sufficient to determine a lower bound on the distance
to all points within the node.
Because of the spherical geometry of the ball tree nodes, its performance
does not degrade at high dimensions, as long as the distribution of
points has a lower-dimensional or sparse structure.
For dense distributions of points that uniformly fill the 
parameter space, the performance gain of a ball tree will be reduced.
Because datasets of interest are usually structured, 
the ball tree is a better choice in practice.

One final performance note is the effect of the 
number of neighbors K.  The performance gain over brute force
for any tree algorithm will degrade as K increases.  Depending on the
size, dimensionality, and structure of the dataset, a search for 
very large K may be performed more efficiently with a brute force algorithm
than with a tree algorithm.

A ball tree can be constructed in a number of ways: the optimal method
depends on the size and structure of the data.  Currently, the
:class:`BallTree` implements only the KD construction algorithm
(see reference below).
This is the simplest and fastest construction algorithm, but is
sub-optimal for some datasets.

.. topic:: References:

   * `"Five balltree construction algorithms"
     <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.91.8209&rep=rep1&type=pdf>`_,
     Omohundro, S.M., International Computer Science Institute
     Technical Report



