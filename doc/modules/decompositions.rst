
====================================================================
Decomposing signals in components (matrix factorization problems)
====================================================================

.. _PCA:

Principal component analysis (PCA)
====================================

.. currentmodule:: scikits.learn.pca

PCA is used to decompose a multivariate dataset in a set of successive
orthogonal components that explain a maximum amount of the variance. In
the scikit-learn, :class:`PCA` is implemented as a `transformer` object
that learns n components in its `fit` method, and can be used on new data
to project it on these components.

In addition, the :class:`ProbabilisticPCA` object provides a
probabilistic interpretation of the PCA that can give a likelihood of
data based on the amount of variance it explains. As such it implements a
`score` method that can be used in cross-validation.

Below is an example of the iris dataset, which is comprised of 4
features, projected on the 2 dimensions that explain most variance:

.. figure:: ../auto_examples/images/plot_pca.png
    :target: ../auto_examples/plot_pca.html
    :align: center
    :scale: 75%

.. topic:: Examples:

    * :ref:`example_plot_pca.py`

.. _ICA:

Independent component analysis (ICA)
=====================================

.. currentmodule:: scikits.learn.fastica

ICA finds components that are maximally independent. It is classically
used to separate mixed signals (a problem know as *blind source
separation*), as in the example below:

.. figure:: ../auto_examples/images/plot_ica_blind_source_separation.png
    :target: ../auto_examples/plot_ica_blind_source_separation.html
    :align: center
    :scale: 50%


.. topic:: Examples:

    * :ref:`example_plot_ica_blind_source_separation.py`
    * :ref:`example_plot_ica_vs_pca.py`

