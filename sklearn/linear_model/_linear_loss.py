"""
Loss functions for linear models with raw_prediction = X @ coef
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from scipy import sparse

from sklearn.utils.extmath import squared_norm


def sandwich_dot(X, W):
    """Compute the sandwich product X.T @ diag(W) @ X."""
    # TODO: This "sandwich product" is the main computational bottleneck for solvers
    # that use the full hessian matrix. Here, thread parallelism would pay-off the
    # most.
    # While a dedicated Cython routine could exploit the symmetry, it is very hard to
    # beat BLAS GEMM, even thought the latter cannot exploit the symmetry, unless one
    # pays the price of taking square roots and implements
    #    sqrtWX = sqrt(W)[: None] * X
    #    return sqrtWX.T @ sqrtWX
    # which (might) detect the symmetry and use BLAS SYRK under the hood.
    n_samples = X.shape[0]
    if sparse.issparse(X):
        return (
            X.T @ sparse.dia_matrix((W, 0), shape=(n_samples, n_samples)) @ X
        ).toarray()
    else:
        # np.einsum may use less memory but the following, using BLAS matrix
        # multiplication (gemm), is by far faster.
        WX = W[:, None] * X
        return X.T @ WX


class LinearModelLoss:
    """General class for loss functions with raw_prediction = X @ coef + intercept.

    Note that raw_prediction is also known as linear predictor.

    The loss is the average of per sample losses and includes a term for L2
    regularization::

        loss = 1 / s_sum * sum_i s_i loss(y_i, X_i @ coef + intercept)
               + 1/2 * l2_reg_strength * ||coef||_2^2

    with sample weights s_i=1 if sample_weight=None and s_sum=sum_i s_i.

    Gradient and hessian, for simplicity without intercept, are::

        gradient = 1 / s_sum * X.T @ loss.gradient + l2_reg_strength * coef
        hessian = 1 / s_sum * X.T @ diag(loss.hessian) @ X
                  + l2_reg_strength * identity

    Conventions:
        if fit_intercept:
            n_dof =  n_features + 1
        else:
            n_dof = n_features

        if base_loss.is_multiclass:
            coef.shape = (n_classes, n_dof) or ravelled (n_classes * n_dof,)
        else:
            coef.shape = (n_dof,)

        The intercept term is at the end of the coef array:
        if base_loss.is_multiclass:
            if coef.shape (n_classes, n_dof):
                intercept = coef[:, -1]
            if coef.shape (n_classes * n_dof,)
                intercept = coef[n_features::n_dof] = coef[(n_dof-1)::n_dof]
            intercept.shape = (n_classes,)
        else:
            intercept = coef[-1]

        Shape of gradient follows shape of coef.
        gradient.shape = coef.shape

        But hessian (to make our lives simpler) are always 2-d:
        if base_loss.is_multiclass:
            hessian.shape = (n_classes * n_dof, n_classes * n_dof)
        else:
            hessian.shape = (n_dof, n_dof)

    Note: If coef has shape (n_classes * n_dof,), the 2d-array can be reconstructed as

        coef.reshape((n_classes, -1), order="F")

    The option order="F" makes coef[:, i] contiguous. This, in turn, makes the
    coefficients without intercept, coef[:, :-1], contiguous and speeds up
    matrix-vector computations.

    Note: If the average loss per sample is wanted instead of the sum of the loss per
    sample, one can simply use a rescaled sample_weight such that
    sum(sample_weight) = 1.

    Parameters
    ----------
    base_loss : instance of class BaseLoss from sklearn._loss.
    fit_intercept : bool
    """

    def __init__(self, base_loss, fit_intercept):
        self.base_loss = base_loss
        self.fit_intercept = fit_intercept

    def init_zero_coef(self, X, dtype=None):
        """Allocate coef of correct shape with zeros.

        Parameters:
        -----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        dtype : data-type, default=None
            Overrides the data type of coef. With dtype=None, coef will have the same
            dtype as X.

        Returns
        -------
        coef : ndarray of shape (n_dof,) or (n_classes, n_dof)
            Coefficients of a linear model.
        """
        n_features = X.shape[1]
        n_classes = self.base_loss.n_classes
        if self.fit_intercept:
            n_dof = n_features + 1
        else:
            n_dof = n_features
        if self.base_loss.is_multiclass:
            coef = np.zeros_like(X, shape=(n_classes, n_dof), dtype=dtype, order="F")
        else:
            coef = np.zeros_like(X, shape=n_dof, dtype=dtype)
        return coef

    def weight_intercept(self, coef):
        """Helper function to get coefficients and intercept.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").

        Returns
        -------
        weights : ndarray of shape (n_features,) or (n_classes, n_features)
            Coefficients without intercept term.
        intercept : float or ndarray of shape (n_classes,)
            Intercept terms.
        """
        if not self.base_loss.is_multiclass:
            if self.fit_intercept:
                intercept = coef[-1]
                weights = coef[:-1]
            else:
                intercept = 0.0
                weights = coef
        else:
            # reshape to (n_classes, n_dof)
            if coef.ndim == 1:
                weights = coef.reshape((self.base_loss.n_classes, -1), order="F")
            else:
                weights = coef
            if self.fit_intercept:
                intercept = weights[:, -1]
                weights = weights[:, :-1]
            else:
                intercept = 0.0

        return weights, intercept

    def weight_intercept_raw(self, coef, X):
        """Helper function to get coefficients, intercept and raw_prediction.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        Returns
        -------
        weights : ndarray of shape (n_features,) or (n_classes, n_features)
            Coefficients without intercept term.
        intercept : float or ndarray of shape (n_classes,)
            Intercept terms.
        raw_prediction : ndarray of shape (n_samples,) or \
            (n_samples, n_classes)
        """
        weights, intercept = self.weight_intercept(coef)

        if not self.base_loss.is_multiclass:
            raw_prediction = X @ weights + intercept
        else:
            # weights has shape (n_classes, n_dof)
            raw_prediction = X @ weights.T + intercept  # ndarray, likely C-contiguous

        return weights, intercept, raw_prediction

    def l2_penalty(self, weights, l2_reg_strength):
        """Compute L2 penalty term l2_reg_strength/2 *||w||_2^2."""
        norm2_w = weights @ weights if weights.ndim == 1 else squared_norm(weights)
        return 0.5 * l2_reg_strength * norm2_w

    def loss(
        self,
        coef,
        X,
        y,
        sample_weight=None,
        l2_reg_strength=0.0,
        n_threads=1,
        raw_prediction=None,
    ):
        """Compute the loss as weighted average over point-wise losses.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        y : contiguous array of shape (n_samples,)
            Observed, true target values.
        sample_weight : None or contiguous array of shape (n_samples,), default=None
            Sample weights.
        l2_reg_strength : float, default=0.0
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.
        raw_prediction : C-contiguous array of shape (n_samples,) or array of \
            shape (n_samples, n_classes)
            Raw prediction values (in link space). If provided, these are used. If
            None, then raw_prediction = X @ coef + intercept is calculated.

        Returns
        -------
        loss : float
            Weighted average of losses per sample, plus penalty.
        """
        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)

        loss = self.base_loss.loss(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=None,
            n_threads=n_threads,
        )
        loss = np.average(loss, weights=sample_weight)

        return loss + self.l2_penalty(weights, l2_reg_strength)

    def loss_gradient(
        self,
        coef,
        X,
        y,
        sample_weight=None,
        l2_reg_strength=0.0,
        n_threads=1,
        raw_prediction=None,
    ):
        """Computes the sum of loss and gradient w.r.t. coef.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        y : contiguous array of shape (n_samples,)
            Observed, true target values.
        sample_weight : None or contiguous array of shape (n_samples,), default=None
            Sample weights.
        l2_reg_strength : float, default=0.0
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.
        raw_prediction : C-contiguous array of shape (n_samples,) or array of \
            shape (n_samples, n_classes)
            Raw prediction values (in link space). If provided, these are used. If
            None, then raw_prediction = X @ coef + intercept is calculated.

        Returns
        -------
        loss : float
            Weighted average of losses per sample, plus penalty.

        gradient : ndarray of shape coef.shape
             The gradient of the loss.
        """
        (n_samples, n_features), n_classes = X.shape, self.base_loss.n_classes
        n_dof = n_features + int(self.fit_intercept)

        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)

        loss, grad_pointwise = self.base_loss.loss_gradient(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )
        sw_sum = n_samples if sample_weight is None else np.sum(sample_weight)
        loss = loss.sum() / sw_sum
        loss += self.l2_penalty(weights, l2_reg_strength)

        grad_pointwise /= sw_sum

        if not self.base_loss.is_multiclass:
            grad = np.empty_like(coef, dtype=weights.dtype)
            grad[:n_features] = X.T @ grad_pointwise + l2_reg_strength * weights
            if self.fit_intercept:
                grad[-1] = grad_pointwise.sum()
        else:
            grad = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
            # grad_pointwise.shape = (n_samples, n_classes)
            grad[:, :n_features] = grad_pointwise.T @ X + l2_reg_strength * weights
            if self.fit_intercept:
                grad[:, -1] = grad_pointwise.sum(axis=0)
            if coef.ndim == 1:
                grad = grad.ravel(order="F")

        return loss, grad

    def gradient(
        self,
        coef,
        X,
        y,
        sample_weight=None,
        l2_reg_strength=0.0,
        n_threads=1,
        raw_prediction=None,
    ):
        """Computes the gradient w.r.t. coef.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        y : contiguous array of shape (n_samples,)
            Observed, true target values.
        sample_weight : None or contiguous array of shape (n_samples,), default=None
            Sample weights.
        l2_reg_strength : float, default=0.0
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.
        raw_prediction : C-contiguous array of shape (n_samples,) or array of \
            shape (n_samples, n_classes)
            Raw prediction values (in link space). If provided, these are used. If
            None, then raw_prediction = X @ coef + intercept is calculated.

        Returns
        -------
        gradient : ndarray of shape coef.shape
             The gradient of the loss.
        """
        (n_samples, n_features), n_classes = X.shape, self.base_loss.n_classes
        n_dof = n_features + int(self.fit_intercept)

        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)

        grad_pointwise = self.base_loss.gradient(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )
        sw_sum = n_samples if sample_weight is None else np.sum(sample_weight)
        grad_pointwise /= sw_sum

        if not self.base_loss.is_multiclass:
            grad = np.empty_like(coef, dtype=weights.dtype)
            grad[:n_features] = X.T @ grad_pointwise + l2_reg_strength * weights
            if self.fit_intercept:
                grad[-1] = grad_pointwise.sum()
            return grad
        else:
            grad = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
            # gradient.shape = (n_samples, n_classes)
            grad[:, :n_features] = grad_pointwise.T @ X + l2_reg_strength * weights
            if self.fit_intercept:
                grad[:, -1] = grad_pointwise.sum(axis=0)
            if coef.ndim == 1:
                return grad.ravel(order="F")
            else:
                return grad

    def gradient_hessian(
        self,
        coef,
        X,
        y,
        sample_weight=None,
        l2_reg_strength=0.0,
        n_threads=1,
        gradient_out=None,
        hessian_out=None,
        raw_prediction=None,
    ):
        """Computes gradient and hessian w.r.t. coef.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        y : contiguous array of shape (n_samples,)
            Observed, true target values.
        sample_weight : None or contiguous array of shape (n_samples,), default=None
            Sample weights.
        l2_reg_strength : float, default=0.0
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.
        gradient_out : None or ndarray of shape coef.shape
            A location into which the gradient is stored. If None, a new array
            might be created.
        hessian_out : None or ndarray of shape (n_dof, n_dof) or \
            (n_classes * n_dof, n_classes * n_dof)
            A location into which the hessian is stored. If None, a new array
            might be created.
        raw_prediction : C-contiguous array of shape (n_samples,) or array of \
            shape (n_samples, n_classes)
            Raw prediction values (in link space). If provided, these are used. If
            None, then raw_prediction = X @ coef + intercept is calculated.

        Returns
        -------
        gradient : ndarray of shape coef.shape
             The gradient of the loss.

        hessian : ndarray of shape (n_dof, n_dof) or \
            (n_classes, n_dof, n_dof, n_classes)
            Hessian matrix.

        hessian_warning : bool
            True if pointwise hessian has more than 25% of its elements non-positive.
        """
        (n_samples, n_features), n_classes = X.shape, self.base_loss.n_classes
        n_dof = n_features + int(self.fit_intercept)
        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)
        sw_sum = n_samples if sample_weight is None else np.sum(sample_weight)

        # Allocate gradient.
        if gradient_out is None:
            grad = np.empty_like(coef, dtype=weights.dtype, order="F")
        elif gradient_out.shape != coef.shape:
            raise ValueError(
                f"gradient_out is required to have shape coef.shape = {coef.shape}; "
                f"got {gradient_out.shape}."
            )
        elif self.base_loss.is_multiclass and not gradient_out.flags.f_contiguous:
            raise ValueError("gradient_out must be F-contiguous.")
        else:
            grad = gradient_out
        # Allocate hessian.
        n = coef.size  # for multinomial this equals n_dof * n_classes
        if hessian_out is None:
            hess = np.empty((n, n), dtype=weights.dtype)
        elif hessian_out.shape != (n, n):
            raise ValueError(
                f"hessian_out is required to have shape ({n, n}); got "
                f"{hessian_out.shape=}."
            )
        elif self.base_loss.is_multiclass and (
            not hessian_out.flags.c_contiguous and not hessian_out.flags.f_contiguous
        ):
            raise ValueError("hessian_out must be contiguous.")
        else:
            hess = hessian_out

        if not self.base_loss.is_multiclass:
            grad_pointwise, hess_pointwise = self.base_loss.gradient_hessian(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad_pointwise /= sw_sum
            hess_pointwise /= sw_sum

            # For non-canonical link functions and far away from the optimum, the
            # pointwise hessian can be negative. We take care that 75% of the hessian
            # entries are positive.
            hessian_warning = (
                np.average(hess_pointwise <= 0, weights=sample_weight) > 0.25
            )
            hess_pointwise = np.abs(hess_pointwise)

            grad[:n_features] = X.T @ grad_pointwise + l2_reg_strength * weights
            if self.fit_intercept:
                grad[-1] = grad_pointwise.sum()

            if hessian_warning:
                # Exit early without computing the hessian.
                return grad, hess, hessian_warning

            hess[:n_features, :n_features] = sandwich_dot(X, hess_pointwise)

            if l2_reg_strength > 0:
                # The L2 penalty enters the Hessian on the diagonal only. To add those
                # terms, we use a flattened view of the array.
                order = "C" if hess.flags.c_contiguous else "F"
                hess.reshape(-1, order=order)[: (n_features * n_dof) : (n_dof + 1)] += (
                    l2_reg_strength
                )

            if self.fit_intercept:
                # With intercept included as added column to X, the hessian becomes
                # hess = (X, 1)' @ diag(h) @ (X, 1)
                #      = (X' @ diag(h) @ X, X' @ h)
                #        (           h @ X, sum(h))
                # The left upper part has already been filled, it remains to compute
                # the last row and the last column.
                Xh = X.T @ hess_pointwise
                hess[:-1, -1] = Xh
                hess[-1, :-1] = Xh
                hess[-1, -1] = hess_pointwise.sum()
        else:
            # Here we may safely assume HalfMultinomialLoss aka categorical
            # cross-entropy.
            # HalfMultinomialLoss computes only the diagonal part of the hessian, i.e.
            # diagonal in the classes. Here, we want the full hessian. Therefore, we
            # call gradient_proba.
            grad_pointwise, proba = self.base_loss.gradient_proba(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad_pointwise /= sw_sum
            grad = grad.reshape((n_classes, n_dof), order="F")
            grad[:, :n_features] = grad_pointwise.T @ X + l2_reg_strength * weights
            if self.fit_intercept:
                grad[:, -1] = grad_pointwise.sum(axis=0)
            if coef.ndim == 1:
                grad = grad.ravel(order="F")

            # The full hessian matrix, i.e. not only the diagonal part, dropping most
            # indices, is given by:
            #
            #   hess = X' @ h @ X
            #
            # Here, h is a priori a 4-dimensional matrix of shape
            # (n_samples, n_samples, n_classes, n_classes). It is diagonal its first
            # two dimensions (the ones with n_samples), i.e. it is
            # effectively a 3-dimensional matrix (n_samples, n_classes, n_classes).
            #
            #   h = diag(p) - p' p
            #
            # or with indices k and l for classes
            #
            #   h_kl = p_k * delta_kl - p_k * p_l
            #
            # with p_k the (predicted) probability for class k. Only the dimension in
            # n_samples multiplies with X.
            # For 3 classes and n_samples = 1, this looks like ("@" is a bit misused
            # here):
            #
            #   hess = X' @ (h00 h10 h20) @ X
            #               (h10 h11 h12)
            #               (h20 h12 h22)
            #        = (X' @ diag(h00) @ X, X' @ diag(h10), X' @ diag(h20))
            #          (X' @ diag(h10) @ X, X' @ diag(h11), X' @ diag(h12))
            #          (X' @ diag(h20) @ X, X' @ diag(h12), X' @ diag(h22))
            #
            # Now coef of shape (n_classes * n_dof) is contiguous in n_classes.
            # Therefore, we want the hessian to follow this convention, too, i.e.
            #     hess[:n_classes, :n_classes] = (x0' @ h00 @ x0, x0' @ h10 @ x0, ..)
            #                                    (x0' @ h10 @ x0, x0' @ h11 @ x0, ..)
            #                                    (x0' @ h20 @ x0, x0' @ h12 @ x0, ..)
            # is the first feature, x0, for all classes. In our implementation, we
            # still want to take advantage of BLAS "X.T @ X". Therefore, we have some
            # index/slicing battle to fight.
            if sample_weight is not None:
                sw = sample_weight / sw_sum
            else:
                sw = 1.0 / sw_sum

            for k in range(n_classes):
                # Diagonal terms (in classes) hess_kk.
                # Note that this also writes to some of the lower triangular part.
                h = proba[:, k] * (1 - proba[:, k]) * sw
                hess[
                    k : n_classes * n_features : n_classes,
                    k : n_classes * n_features : n_classes,
                ] = sandwich_dot(X, h)
                if self.fit_intercept:
                    # See above in the non multiclass case.
                    Xh = X.T @ h
                    hess[
                        k : n_classes * n_features : n_classes,
                        n_classes * n_features + k,
                    ] = Xh
                    hess[
                        n_classes * n_features + k,
                        k : n_classes * n_features : n_classes,
                    ] = Xh
                    hess[n_classes * n_features + k, n_classes * n_features + k] = (
                        h.sum()
                    )
                # Off diagonal terms (in classes) hess_kl.
                for l in range(k + 1, n_classes):
                    # Upper triangle (in classes).
                    h = -proba[:, k] * proba[:, l] * sw
                    hess[
                        k : n_classes * n_features : n_classes,
                        l : n_classes * n_features : n_classes,
                    ] = sandwich_dot(X, h)
                    if self.fit_intercept:
                        Xh = X.T @ h
                        hess[
                            k : n_classes * n_features : n_classes,
                            n_classes * n_features + l,
                        ] = Xh
                        hess[
                            n_classes * n_features + k,
                            l : n_classes * n_features : n_classes,
                        ] = Xh
                        hess[n_classes * n_features + k, n_classes * n_features + l] = (
                            h.sum()
                        )
                    # Fill lower triangle (in classes).
                    hess[l::n_classes, k::n_classes] = hess[k::n_classes, l::n_classes]

            if l2_reg_strength > 0:
                # See above in the non multiclass case.
                order = "C" if hess.flags.c_contiguous else "F"
                hess.reshape(-1, order=order)[
                    : (n_classes**2 * n_features * n_dof) : (n_classes * n_dof + 1)
                ] += l2_reg_strength

            # The pointwise hessian is always non-negative for the multinomial loss.
            hessian_warning = False

        return grad, hess, hessian_warning

    def gradient_hessian_product(
        self, coef, X, y, sample_weight=None, l2_reg_strength=0.0, n_threads=1
    ):
        """Computes gradient and hessp (hessian product function) w.r.t. coef.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        y : contiguous array of shape (n_samples,)
            Observed, true target values.
        sample_weight : None or contiguous array of shape (n_samples,), default=None
            Sample weights.
        l2_reg_strength : float, default=0.0
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.

        Returns
        -------
        gradient : ndarray of shape coef.shape
             The gradient of the loss.

        hessp : callable
            Function that takes in a vector input of shape of gradient and
            and returns matrix-vector product with hessian.
        """
        (n_samples, n_features), n_classes = X.shape, self.base_loss.n_classes
        n_dof = n_features + int(self.fit_intercept)
        weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        sw_sum = n_samples if sample_weight is None else np.sum(sample_weight)

        if not self.base_loss.is_multiclass:
            grad_pointwise, hess_pointwise = self.base_loss.gradient_hessian(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad_pointwise /= sw_sum
            hess_pointwise /= sw_sum
            grad = np.empty_like(coef, dtype=weights.dtype)
            grad[:n_features] = X.T @ grad_pointwise + l2_reg_strength * weights
            if self.fit_intercept:
                grad[-1] = grad_pointwise.sum()

            # Precompute as much as possible: hX, hX_sum and hessian_sum
            hessian_sum = hess_pointwise.sum()
            if sparse.issparse(X):
                hX = (
                    sparse.dia_matrix((hess_pointwise, 0), shape=(n_samples, n_samples))
                    @ X
                )
            else:
                hX = hess_pointwise[:, np.newaxis] * X

            if self.fit_intercept:
                # Calculate the double derivative with respect to intercept.
                # Note: In case hX is sparse, hX.sum is a matrix object.
                hX_sum = np.squeeze(np.asarray(hX.sum(axis=0)))
                # prevent squeezing to zero-dim array if n_features == 1
                hX_sum = np.atleast_1d(hX_sum)

            # With intercept included and l2_reg_strength = 0, hessp returns
            # res = (X, 1)' @ diag(h) @ (X, 1) @ s
            #     = (X, 1)' @ (hX @ s[:n_features], sum(h) * s[-1])
            # res[:n_features] = X' @ hX @ s[:n_features] + sum(h) * s[-1]
            # res[-1] = 1' @ hX @ s[:n_features] + sum(h) * s[-1]
            def hessp(s):
                ret = np.empty_like(s)
                if sparse.issparse(X):
                    ret[:n_features] = X.T @ (hX @ s[:n_features])
                else:
                    ret[:n_features] = np.linalg.multi_dot([X.T, hX, s[:n_features]])
                ret[:n_features] += l2_reg_strength * s[:n_features]

                if self.fit_intercept:
                    ret[:n_features] += s[-1] * hX_sum
                    ret[-1] = hX_sum @ s[:n_features] + hessian_sum * s[-1]
                return ret

        else:
            # Here we may safely assume HalfMultinomialLoss aka categorical
            # cross-entropy.
            # HalfMultinomialLoss computes only the diagonal part of the hessian, i.e.
            # diagonal in the classes. Here, we want the matrix-vector product of the
            # full hessian. Therefore, we call gradient_proba.
            grad_pointwise, proba = self.base_loss.gradient_proba(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad_pointwise /= sw_sum
            grad = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
            grad[:, :n_features] = grad_pointwise.T @ X + l2_reg_strength * weights
            if self.fit_intercept:
                grad[:, -1] = grad_pointwise.sum(axis=0)

            # Full hessian-vector product, i.e. not only the diagonal part of the
            # hessian. Derivation with some index battle for input vector s:
            #   - sample index i
            #   - feature indices j, m
            #   - class indices k, l
            #   - 1_{k=l} is one if k=l else 0
            #   - p_i_k is the (predicted) probability that sample i belongs to class k
            #     for all i: sum_k p_i_k = 1
            #   - s_l_m is input vector for class l and feature m
            #   - X' = X transposed
            #
            # Note: Hessian with dropping most indices is just:
            #       X' @ p_k (1(k=l) - p_l) @ X
            #
            # result_{k j} = sum_{i, l, m} Hessian_{i, k j, m l} * s_l_m
            #   = sum_{i, l, m} (X')_{ji} * p_i_k * (1_{k=l} - p_i_l)
            #                   * X_{im} s_l_m
            #   = sum_{i, m} (X')_{ji} * p_i_k
            #                * (X_{im} * s_k_m - sum_l p_i_l * X_{im} * s_l_m)
            #
            # See also https://github.com/scikit-learn/scikit-learn/pull/3646#discussion_r17461411
            def hessp(s):
                s = s.reshape((n_classes, -1), order="F")  # shape = (n_classes, n_dof)
                if self.fit_intercept:
                    s_intercept = s[:, -1]
                    s = s[:, :-1]  # shape = (n_classes, n_features)
                else:
                    s_intercept = 0
                tmp = X @ s.T + s_intercept  # X_{im} * s_k_m
                tmp += (-proba * tmp).sum(axis=1)[:, np.newaxis]  # - sum_l ..
                tmp *= proba  # * p_i_k
                if sample_weight is not None:
                    tmp *= sample_weight[:, np.newaxis]
                # hess_prod = empty_like(grad), but we ravel grad below and this
                # function is run after that.
                hess_prod = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
                hess_prod[:, :n_features] = (tmp.T @ X) / sw_sum + l2_reg_strength * s
                if self.fit_intercept:
                    hess_prod[:, -1] = tmp.sum(axis=0) / sw_sum
                if coef.ndim == 1:
                    return hess_prod.ravel(order="F")
                else:
                    return hess_prod

            if coef.ndim == 1:
                return grad.ravel(order="F"), hessp

        return grad, hessp
