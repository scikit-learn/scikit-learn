
.. _manifold:

.. currentmodule:: scikits.learn.manifold

=================
Manifold learning
=================

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



Locally Linear Embedding
------------------------

Locally linear embedding can be performed with function
:func:`locally_linear_embedding` or its object-oriented counterpart
:class:`LocallyLinearEmbedding`.

These take as input a set of points in a high-dimensional space and return
those points embedded in a manifold of dimension specified by parameter
``out_dim``.

.. figure:: ../auto_examples/manifold/images/plot_lle_digits_1.png
   :target: ../auto_examples/manifold/plot_lle_digits.html
   :align: center


.. topic:: Examples:

    * See :ref:`example_manifold_plot_lle_digits.py` for an example of
      dimensionality reduction on handwritten digits.

    * See :ref:`example_manifold_plot_swissroll.py` for an example of
      locally linear embedding on the swiss roll.


Complexity
----------

The complete algorithm scales using the `dense` eigensolver scales as::

..math:: O(N log(N)) + O(D N K^3) + O(d N^2)

where N is the number of samples, D is the input dimension, d the output
dimension and K the number of neighbors. If the `lobcpg` solver is used, the
last term can be reduced to sub-quadratic in N.


Tips on practical use
---------------------

* Make sure the same scale is used over all features. Being this a
  nearest-neighbors method it will behave poorly otherwise.

* On certain problems, the `lobcpg` solver might converge slowly. Supply a
generous value for `max_iter` if big oscillations are detected between runs.
