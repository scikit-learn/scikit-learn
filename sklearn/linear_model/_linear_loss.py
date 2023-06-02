"""
Loss functions for linear models with raw_prediction = X @ coef
"""
import warnings

import numpy as np
from scipy import sparse
from ..utils.extmath import squared_norm


class LinearModelLoss:
    """General class for loss functions with raw_prediction = X @ coef + intercept.

    Note that raw_prediction is also known as linear predictor.

    The loss is the sum of per sample losses and includes a term for L2
    regularization::

        loss = sum_i s_i loss(y_i, X_i @ coef + intercept)
               + 1/2 * l2_reg_strength * ||coef||_2^2

    with sample weights s_i=1 if sample_weight=None. For the multiclass case, the
    squared L2 norm is replaced by the squared Frobenius norm, which equals the L2 norm
    of the flattend/ravelled coefficient array.

    Gradient and hessian, for simplicity without intercept, are::

        gradient = X.T @ loss.gradient + l2_reg_strength * coef
        hessian = X.T @ diag(loss.hessian) @ X + l2_reg_strength * identity

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
        """Compute the loss as sum over point-wise losses.

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
            Sum of losses per sample plus penalty.
        """
        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)

        loss = self.base_loss.loss(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )
        loss = loss.sum()

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
            Sum of losses per sample plus penalty.

        gradient : ndarray of shape coef.shape
             The gradient of the loss.
        """
        n_features, n_classes = X.shape[1], self.base_loss.n_classes
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
        loss = loss.sum()
        loss += self.l2_penalty(weights, l2_reg_strength)

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
        n_features, n_classes = X.shape[1], self.base_loss.n_classes
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
        hessian_out : None or ndarray
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

        hessian : ndarray
            Hessian matrix.

        hessian_warning : bool
            True if pointwise hessian has more than half of its elements non-positive.
        """
        n_samples, n_features = X.shape
        n_dof = n_features + int(self.fit_intercept)

        if raw_prediction is None:
            weights, intercept, raw_prediction = self.weight_intercept_raw(coef, X)
        else:
            weights, intercept = self.weight_intercept(coef)

        grad_pointwise, hess_pointwise = self.base_loss.gradient_hessian(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )

        # For non-canonical link functions and far away from the optimum, the pointwise
        # hessian can be negative. We take care that 75% of the hessian entries are
        # positive.
        hessian_warning = np.mean(hess_pointwise <= 0) > 0.25
        hess_pointwise = np.abs(hess_pointwise)

        if not self.base_loss.is_multiclass:
            # gradient
            if gradient_out is None:
                grad = np.empty_like(coef, dtype=weights.dtype)
            else:
                grad = gradient_out
            grad[:n_features] = X.T @ grad_pointwise + l2_reg_strength * weights
            if self.fit_intercept:
                grad[-1] = grad_pointwise.sum()

            # hessian
            if hessian_out is None:
                hess = np.empty(shape=(n_dof, n_dof), dtype=weights.dtype)
            else:
                hess = hessian_out

            if hessian_warning:
                # Exit early without computing the hessian.
                return grad, hess, hessian_warning

            # TODO: This "sandwich product", X' diag(W) X, is the main computational
            # bottleneck for solvers. A dedicated Cython routine might improve it
            # exploiting the symmetry (as opposed to, e.g., BLAS gemm).
            if sparse.issparse(X):
                hess[:n_features, :n_features] = (
                    X.T
                    @ sparse.dia_matrix(
                        (hess_pointwise, 0), shape=(n_samples, n_samples)
                    )
                    @ X
                ).toarray()
            else:
                # np.einsum may use less memory but the following, using BLAS matrix
                # multiplication (gemm), is by far faster.
                WX = hess_pointwise[:, None] * X
                hess[:n_features, :n_features] = np.dot(X.T, WX)

            if l2_reg_strength > 0:
                # The L2 penalty enters the Hessian on the diagonal only. To add those
                # terms, we use a flattened view on the array.
                hess.reshape(-1)[
                    : (n_features * n_dof) : (n_dof + 1)
                ] += l2_reg_strength

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
            raise NotImplementedError

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

        if not self.base_loss.is_multiclass:
            grad_pointwise, hess_pointwise = self.base_loss.gradient_hessian(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
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
            # See also https://github.com/scikit-learn/scikit-learn/pull/3646#discussion_r17461411  # noqa
            def hessp(s):
                # s.shape = (n_classes, n_features,) or (n_classes * n_features)
                # It could be the coefficients (including the intercept if present).
                # Returns an array of shape = (n_classes, n_dof) or (n_classes * n_dof)
                # depending on the shape of s.
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
                hess_prod[:, :n_features] = tmp.T @ X + l2_reg_strength * s
                if self.fit_intercept:
                    hess_prod[:, -1] = tmp.sum(axis=0)
                if coef.ndim == 1:
                    return hess_prod.ravel(order="F")
                else:
                    return hess_prod

            if coef.ndim == 1:
                return grad.ravel(order="F"), hessp

        return grad, hessp


class Multinomial_LDL_Decomposition:
    """A class for symbolic LDL' decomposition of multinomial hessian.

    The pointwise hessian of the multinomial loss with c = n_classes is given by
        h = diag(p) - p' p        or with indices
        h_ij = p_i * delta_ij - p_i * p_j    i and j are indices of classes.
    This holds for every single point (sample). The LDL decomposition
    of this 2-dim matrix is given in [1] for p_i > 0 and sum(p) <= 1 as (math indices)
        h = L D L'   with lower triangular L and diagonal D:: text

        q_0 = 1
        q_i = 1 - sum(p_k, k=1..c)
        L_ii = 1
        for i > j:
            L_ij = -p_i / q_j
        D_ii = p_i * q_i / q_{i-1}

    If the p_i sum to 1, then q_c = 0 and D_cc = 0, with c = n_classes.
    The inverse L^-1 is also lower triangular and given by:: text

        (L^-1)_ii = 1
        for i > j:
            (L^-1)_ij = f_i = p_i / q_{i-1} = -L_{i, i-1}

    Note that with Python indices (zero-based), we have: p_i = p[:, i-1],
    q_i = q[:, i-1] and we don't explicitly need q_0 = 1.

    The trick is to apply this LDL decomposition to all points (samples) at the same
    time. This class therefore provides methods to apply products with the matrices L
    and D to vectors/matrices. All objects have the 0-th dimension for n_samples and
    the 1st dimention for n_classes, i.e. shape (n_samples, n_classes).

    Parameters
    ----------
        proba : ndarray of shape (n_samples, n_classes)
            Array of predicted probabilities per class (and sample).

    Attributes
    ----------
        p : ndarray of shape (n_samples, n_classes)
            Array of predicted probabilities per class (and sample).

        q_inv : ndarray of shape (n_samples, n_classes)
            Helper array, inverse of q[:, j] = 1 - sum_{i=0}^j p[:, i], i.e. 1/q.

        sqrt_d : ndarray of shape (n_samples, n_classes)
            Square root of the diagonal matrix D, D_ii = p_i * q_i / q_{i-1}
            with q_{-1} = 1.

        proba_sum_to_1 : bool
            True if propabilities p sum to 1. This is most often expected to be true.

    References
    ----------
    .. [1] Kunio Tanabe and Masahiko Sagae. (1992) "An Exact Cholesky Decomposition and
           the Generalized Inverse of the Variance-Covariance Matrix of the Multinomial
           Distribution, with Applications"
           https://doi.org/10.1111/j.2517-6161.1992.tb01875.x
    """

    def __init__(self, *, proba, proba_sum_to_1=True):
        self.p = proba
        q = 1 - np.cumsum(self.p, axis=1)  # contiguity of p
        self.proba_sum_to_1 = proba_sum_to_1
        if self.p.dtype == np.float32:
            eps = 2 * np.finfo(np.float32).resolution
        else:
            eps = 2 * np.finfo(np.float64).resolution
        if np.any(q[:, -1] > eps):
            warnings.warn(
                "Probabilities proba are assumed to sum to 1, but they don't."
            )
        if self.proba_sum_to_1:
            # If np.sum(p, axis=1) = 1 then q[:, -1] = d[:, -1] = 0.
            q[:, -1] = 0.0
            # One might not need all classes for p to sum to 1, so we detect and
            # correct all values close to zero.
            q[q <= eps] = 0.0
            self.p[self.p <= eps] = 0.0
        d = self.p * q
        # From now on, q is always used in the denominator. We handle q == 0 by
        # setting q to 1 whenever q == 0 such that a division of q is a no-op in this
        # case.
        q[q == 0] = 1
        # And we use the inverse of q, 1/q.
        self.q_inv = 1 / q
        if self.proba_sum_to_1:
            # If q_{i - 1} = 0, then also q_i = 0.
            d[:, 1:-1] *= self.q_inv[:, :-2]  # d[:, -1] = 0 anyway.
        else:
            d[:, 1:] *= self.q_inv[:, :-1]
        self.sqrt_d = np.sqrt(d)

    def sqrt_D_Lt_matmul(self, x):
        """Compute sqrt(D) L' x from the multinomial LDL' decomposition.

        L' is the transpose of L, Lij is the i-th row and j-th column of L:
            Lij = -p[:, i] / q[:, j]

        For n_classes = 4, L' looks like:: text

            L' = (1 L10 L20 L30)
                 (0   1 L21 L31)
                 (0   0   1 L32)
                 (0   0   0   1)

        The carried out operation is matmul over n_classes (1st dimention) and
        element-wise multiplication over n_samples (0th dimension):

            x_ij = sum_k sqrt_d_{i, j} L'_{i, j, k} x_{i, k}

        Parameters
        ----------
        x : ndarray of shape (n_samples, n_classes)
            This array is overwritten with the result.

        Return
        ------
        x : ndarray of shape (n_samples, n_classes)
            Input array x, filled with the result.
        """
        n_classes = self.p.shape[1]
        for i in range(0, n_classes - 1):  # row i
            # L_ij = -p_i / q_j, we need transpose L'
            # for j in range(i + 1, n_classes):  # column j
            #     x[:, i] -= self.p[:, j] / self.q[:, i] * x[:, j]
            # The following is the same but faster.
            px = np.einsum(
                "ij,ij->i",
                self.p[:, i + 1 : n_classes],
                x[:, i + 1 : n_classes],
                order="A",
            )
            x[:, i] -= px * self.q_inv[:, i]
        x *= self.sqrt_d
        return x

    def L_sqrt_D_matmul(self, x):
        """Compute L sqrt(D) x from the multinomial LDL' decomposition.

        L is lower triangular and given by Lij = -p[:, i] / q[:, j]

        For n_classes = 4, L looks like:: text

            L = (1     0   0 0)
                (L10   1   0 0)
                (L20 L21   1 0)
                (L30 L31 L32 1)

        The carried out operation is matmul over n_classes (1st dimention) and
        element-wise multiplication over n_samples (0th dimension):

            x_ij = sum_k L'_{i, j, k} sqrt_d_{i, k}  x_{i, k}

        Parameters
        ----------
        x : ndarray of shape (n_samples, n_classes)
            This array is overwritten with the result.

        Return
        ------
        x : ndarray of shape (n_samples, n_classes)
            Input array x, filled with the result.
        """
        n_classes = self.p.shape[1]
        x *= self.sqrt_d
        for i in range(n_classes - 1, 0, -1):  # row i
            # L_ij = -p_i / q_j
            # for j in range(0, i):  # column j
            #     x[:, i] -= self.p[:, i] / self.q[:, j] * x[:, j]
            # The following is the same but faster.
            qx = np.einsum("ij,ij->i", self.q_inv[:, :i], x[:, :i], order="A")
            x[:, i] -= qx * self.p[:, i]
        return x

    def inverse_L_sqrt_D_matmul(self, x):
        """Compute 1/sqrt(D) L^(-1) x from the multinomial LDL' decomposition.

        L^(-1) is again lower triangular and given by:
            L^(-1)_ij = f[:, i] = p[:, i] / q[:, i-1]

        For n_classes = 4, L^(-1) looks like:: text

            L^(-1) = (1   0  0  0)
                     (f1  1  0  0)
                     (f2  f2 1  0)
                     (f3  f3 f3 1)

        Note that (L sqrt(D))^(-1) = 1/sqrt(D) L^(-1)

        Parameters
        ----------
        x : ndarray of shape (n_samples, n_classes)
            This array is overwritten with the result.

        Return
        ------
        x : ndarray of shape (n_samples, n_classes)
            Input array x, filled with the result.
        """
        n_classes = self.p.shape[1]
        for i in range(n_classes - 1, 0, -1):  # row i
            if i > 0:
                fj = self.p[:, i] * self.q_inv[:, i - 1]
            else:
                fj = self.p[:, i]
            # for j in range(0, i):  # column j
            #     x[:, i] += fj * x[:, j]
            # The following is the same but faster.
            x[:, i] += fj * np.sum(x[:, :i], axis=1)
        if self.proba_sum_to_1:
            # x[:, :-1] /= self.sqrt_d[:, :-1]
            mask = self.sqrt_d[:, :-1] == 0
            x[:, :-1] *= (~mask) / (self.sqrt_d[:, :-1] + mask)
            # Important Note:
            # Strictly speaking, the inverse of D does not exist.
            # We use 0 as the inverse of 0 and just set:
            # x[:, -1] = 0
            x[:, -1] = 0
        else:
            x /= self.sqrt_d
        return x
