
.. _manifold:

.. currentmodule:: scikits.learn.manifold

=================
Manifold learning
=================

.. rst-class:: quote

                 | Look for the bare necessities
                 | The simple bare necessities
                 | Forget about your worries and your strife
                 | I mean the bare necessities
                 | Old Mother Nature's recipes
                 | That bring the bare necessities of life
                 |
                 |             -- Baloo's song [The Jungle Book]



Manifold learning is an approach to nonlinear dimensionality reduction.
Algorithms for this task are based on the idea that the dimensionality of
many data sets is only artificially high.


.. figure:: ../auto_examples/manifold/images/plot_compare_methods.png
   :target: ../auto_examples/manifold/plot_compare_methods.html
   :align: center


Locally Linear Embedding
========================

Locally linear embedding (LLE) seeks a lower-dimensional projection of the data
which preserves distances within local neighborhoods.  It can be thought
of as a series of local *Principal Component* analyses which are optimized
to find the best nonlinear embedding.

Locally linear embedding can be performed with function
:func:`locally_linear_embedding` or its object-oriented counterpart
:class:`LocallyLinearEmbedding`.

``scikits.learn`` contains both standard LLE and several variants, which
can be accessed via the ``method`` keyword in either of the above routines.
They are:
* ``method = 'standard'``: Standard LLE.  This is the most computationally
  efficient of the LLE algorithms, but because of the necessity to
  regularize the local weight matrix, often leads to projections which
  distort the geometry of the embedded manifold.
* ``method = 'hessian'``: Hessian LLE, also known as Hessian Eigenmapping.
  This method addresses some of the distortion problems of standard LLE
  by constructing a more complicated *hessian estimator* in each local
  neighborhood. The size of the hessian estimator scales as
  :math:`O[d^2]`, where :math:`d` is the dimension of the embedded manifold.
  additionally, the number of nearest neighbors must be greater than
  :math:`d(d+3)/2`.  Because of this, hessian LLE is not suitable for
  some problems.
* ``method = 'modified'``: Modified LLE. This routine uses multiple weight 
  vectors within each local neighborhood to solve the regularization problem
  of standard LLE.  While this adds to the cost of computing the embedding,
  it is more efficient than the hessian LLE algorithm.
* ``method = 'ltsa'``: Local Tangent Space Alignment.  Though philosophically
  distinct from LLE, LTSA is included under LLE because it is algorithmically
  very similar to LLE and its variants.  LTSA seeks subspaces which are
  tangent to the data within local neighborhoods, and then performs a global
  analysis to align these spaces.  Its compuational cost is similar to that
  of Modified LLE.

These take as input a set of points in a high-dimensional space and return
those points embedded in a manifold of dimension specified by parameter
``out_dim``.

.. figure:: ../auto_examples/manifold/images/plot_lle_digits_3.png
   :target: ../auto_examples/manifold/plot_lle_digits.html
   :align: center


.. topic:: Examples:

    * See :ref:`example_manifold_plot_lle_digits.py` for an example of
      dimensionality reduction on handwritten digits.

    * See :ref:`example_manifold_plot_compare_methods.py` for an example of
      dimensionality reduction on a toy "S-curve" dataset.


Complexity
----------

For standard LLE, the complete algorithm scales using the 
`dense` eigensolver scales as
:math:`O(N log(N)) + O(D N K^3) + O(d N^2)`, where N is the number of samples,
D is the input dimension, d the output dimension and K the number of
neighbors. If the `lobcpg` or `arpack` solver is used, 
the last term can be reduced to sub-quadratic in N.

Isomap
======

Another approach to manifold learning is the Isomap algorithm, short for
Isometric Mapping.  Isomap seeks a lower-dimensional embedding which
maintains geodesic distances between all points.

Isomap can be performed with the function :func:`scikits.learn.manifold.isomap`
or the object :class:`scikits.learn.manifold.Isomap`.

Complexity
----------
The optimized isomap algorithm has a full cost of 
:math:`O[N log(N)] + O[N^2 (k + log(N))] + O[d N^2]` where the final
cost can be reduced to sub-quadratic through the use of ``ARPACK``.
The main bottleneck is a shortest-path graph search of all the points.
For this, `scikit-learn` uses either *Dijkstra's algorithm* with
*Fibonacci trees*, or the *Floyd-Warshall algorithm*, depending on the
parameters of the problem.


Tips on practical use
=====================

* Make sure the same scale is used over all features. Being this a
  nearest-neighbors method it will behave poorly otherwise.

* On certain problems, the `lobcpg` solver might converge slowly. Supply a
  generous value for `max_iter` if big oscillations are detected between runs.
