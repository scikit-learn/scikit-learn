
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
                 |             -- Ballo's song [The Jungle Book]



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


Tips on practical use
---------------------

* Make sure the same scale is used over all features. Being this a
  nearest-neighbors method it will behave poorly otherwise.
