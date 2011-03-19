
.. _notes_neighbors:


.. currentmodule:: scikits.learn.neighbors

=====================================
scikits.learn.neighbors working notes
=====================================

barycenter
==========

Function :func:`barycenter` tries to find appropriate weights to
reconstruct the point x from a subset (y1, y2, ..., yn), where weights
sum to one.

This is just a simple case of Equality Constrained Least Squares
[#f1]_ with constrain dot(np.ones(n), x) = 1. In particular, the Q
matrix from the QR decomposition of B is the Householder reflection of
np.ones(n).


Purpose
-------

This method was added to ease some computations in the future manifold
module, namely in LLE. However, it is still to be shown that it is
useful and efficient in that context.


Performance
-----------

The algorithm has to iterate over n_samples, which is the main
bottleneck. It would be great to vectorize this loop. Also, the rank
updates could probably be moved outside the loop.

Also, least squares solution could be computed more efficiently by a
QR factorization, since probably we don't care about a minimum norm
solution for the undertermined case.

The paper 'An introduction to Locally Linear Embeddings', Saul &
Roweis solves the problem by the normal equation method over the
covariance matrix. However, it does not degrade grathefully when the
covariance is singular, requiring to explicitly add regularization.


Stability
---------

Should be good as it uses SVD to solve the LS problem. TODO: explicit
bounds.


API
---

The API is convenient to use from NeighborsBarycenter and
kneighbors_graph, but might not be very easy to use directly due to
the fact that Y must be a 3-D array.

It should be checked that it is usable in other contexts.


.. rubric:: Footnotes

.. [#f1] Section 12.1.4 ('Equality Constrained Least Squares'),
         'Matrix Computations' by Golub & Van Loan 
