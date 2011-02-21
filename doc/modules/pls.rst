.. _pls:

======================
Partial Least Squares
======================

.. currentmodule:: scikits.learn.pls

PLS is a used to find the fundamental relations between two matrices
(X and Y), i.e. a latent variable approach to modeling the covariance
structures in these two spaces. A PLS model will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. PLS-regression is
particularly suited when the matrix of predictors has more variables
than observations, and when there is multicollinearity among X
values. By contrast, standard regression will fail in these cases.

Classes included in this module are :class:`PLSRegression`
:class:`PLSCanonical`, :class:`CCA` and :class:`PLSSVD`
