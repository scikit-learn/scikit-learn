=================
Nearest Neighbors
=================

The principle behind nearest neighbor methods is to find the k
training points closest in euclidean distance to x0, and then classify
using a majority vote among the k neighbors.

Despite its simplicity, nearest neighbors has been successful in a
large number of classification problems, including handwritten digits
or satellite image scenes. It is ofter successful in situation where
the decission boundary is very irregular.

Classification
==============


.. autoclass:: scikits.learn.neighbors.Neighbors
   :members:

Regression
==========

Nearest neighbor regression is not (yet) implemented, yet it should be
straightforward using the BallTree class.

Examples
========

See :ref:`example_plot_neighbors.py` for an example of classification
using nearest neighbor.

BallTree
========

Behind the scenes, nearest neighbor search is done by the object
BallTree, which is a fast way to perform neighbor searches in data
sets of very high dimensionality.

.. autoclass:: scikits.learn.BallTree.BallTree
   :members:

.. autofunction:: scikits.learn.BallTree.knn_brute
