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

The :class:`Neighbors` estimators implements the nearest-neighbors
method.

.. topic:: Examples:

  * :ref:`example_plot_neighbors.py`: an example of classification
    using nearest neighbor.

Regression
==========

Nearest neighbor regression is not (yet) implemented, yet it should be
straightforward using the BallTree class.


.. currentmodule:: scikits.learn.ball_tree

Efficient implementation: the ball tree
==========================================

Behind the scenes, nearest neighbor search is done by the object
:class:`BallTree`, which is a fast way to perform neighbor searches in data
sets of very high dimensionality.

This class provides an interface to an optimized BallTree
implementation to rapidly look up the nearest neighbors of any point.
Ball Trees are particularly useful for very high-dimensionality data,
where more traditional tree searches (e.g. KD-Trees) perform poorly.

The cost is a slightly longer construction time, though for repeated
queries, this added construction time quickly becomes insignificant.

A Ball Tree reduces the number of candidate points for a neighbor search
through use of the triangle inequality:

.. math::   |x+y| \leq |x| + |y|

Each node of the :class:`BallTree` defines a centroid, C, and a radius r such
that each point in the node lies within the hyper-sphere of radius r,
centered at C.  With this setup, a single distance calculation between
a test point and the centroid is sufficient to determine a lower bound
on the distance to all points within the node.  Carefully taking
advantage of this property leads to determining neighbors in O[log(N)]
time, as opposed to O[N] time for a brute-force search.

A pure C implementation of brute-force search is also provided in
function :func:`knn_brute`.

.. topic:: References:

   * `"Five balltree construction algorithms"
     <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.91.8209&rep=rep1&type=pdf>`_,
     Omohundro, S.M., International Computer Science Institute
     Technical Report
     


