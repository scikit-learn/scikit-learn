.. _cross_decomposition:

===================
Cross decomposition
===================

.. currentmodule:: sklearn.cross_decomposition

The cross decomposition module contains two main families of algorithms: the
partial least squares (PLS) and the canonical correlation analysis (CCA).

These families of algorithms are useful to find linear relations between two
multivariate datasets: the ``X`` and ``Y`` arguments of the ``fit`` method
are 2D arrays.

.. figure:: ../auto_examples/cross_decomposition/images/sphx_glr_plot_compare_cross_decomposition_001.png
   :target: ../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html
   :scale: 75%
   :align: center


Cross decomposition algorithms find the fundamental relations between two
matrices (X and Y). They are latent variable approaches to modeling the
covariance structures in these two spaces. They will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. PLS-regression is
particularly suited when the matrix of predictors has more variables than
observations, and when there is multicollinearity among X values. By contrast,
standard regression will fail in these cases.

Classes included in this module are :class:`PLSRegression`
:class:`PLSCanonical`, :class:`CCA` and :class:`PLSSVD`


.. topic:: Reference:

   * JA Wegelin
     `A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block case <https://www.stat.washington.edu/research/reports/2000/tr371.pdf>`_

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_cross_decomposition_plot_compare_cross_decomposition.py`



- PLS1: single output
- PLS2: multioutput

They're all transformers (they optionally use Y tho)
All can predict except for PLSSVD


PLSRegression:
   deflation_mode="regression", mode="A", norm_y_weights=False, algorith='nipals'
PLSCAnonincal
   deflation_mode="canonical", mode="A", norm_y_weights=True,
   algorithm=algorithm
CCA:
   deflation_mode="canonical", mode="B", norm_y_weights=True, algorithm="nipals"
PLSSVD
   mentioned slide 19
   http://vision.cse.psu.edu/seminars/talks/PLSpresentation.pdf, called
   Bookstein
   Also discussed in survey from Wegelin. says it can be used for predicting
   too

   Can only transform data, not predict


So norm_y_weights = (deflation == canonical)???

2 block vs.... ????? I think they're all 2-block

With a single output (make_regression, bigger or lower n_features):
- PLSRegression: OK, good score
- PLSCanonical: warning, bad score
- CCA: warning, good score
With multioutpt: no warnings. results OK


Look at slide 17 http://www.eigenvector.com/Docs/Wise_pls_properties.pdf for
tests about orthogonality and stuff

PLS1 i.e. PLSReg with single target is a regularization technique similar to
ridge.

They should all be equivalent when n_component == 1 (wegelin). (maybe not CCA
tho)