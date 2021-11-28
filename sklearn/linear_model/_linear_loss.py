"""
Loss functions for linear models with raw_prediction = X @ coef
"""
import numpy as np
from scipy import sparse
from ..utils.extmath import squared_norm


class LinearLoss:
    """General class for loss functions with raw_prediction = X @ coef.

    The loss is the sum of per sample losses and includes an L2 term::

        loss = sum_i s_i loss(y_i, X_i @ coef + intercept) + 1/2 * alpha * ||coef||_2^2

    with sample weights s_i=1 if sample_weight=None.

    Gradient and hessian, for simplicity without intercept, are::

        gradient = X.T @ loss.gradient + alpha * coef
        hessian = X.T @ diag(loss.hessian) @ X + alpha * identity

    Conventions:
        if fit_intercept:
            n_dof =  n_features + 1
        else:
            n_dof = n_features

        if loss.is_multiclass:
            coef.shape = (n_classes * n_dof,)
            intercept.shape = (n_classes)
        else:
            coef.shape = (n_dof,)
            intercept.shape = (1) or it is a float

        The intercept term is at the end of the coef array:
        if loss.is_multiclass:
            coef[n_features::n_dof] = coef[(n_dof-1)::n_dof]
        else:
            coef[-1]

    Note: If the average loss per sample is wanted instead of the sum of the
    loss per sample, one can simply use a rescaled sample_weight such that
    sum(sample_weight) = 1.

    Parameters
    ----------
    _loss : instance of a loss function from sklearn._loss.
    fit_intercept : bool
    """

    def __init__(self, loss, fit_intercept):
        self._loss = loss
        self.fit_intercept = fit_intercept

    def _w_intercept_raw(self, coef, X):
        """Helper function to get coefficients, intercept and raw_prediction.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,) or (n_classes * n_dof,)
            Coefficients of a linear model.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        Returns
        -------
        w : ndarray of shape (n_features,) or (n_classes, n_features)
            Coefficients without intercept term.
        intercept : float or ndarray of shape (n_classes,)
            Intercept terms.
        raw_prediction : ndarray of shape (n_samples,) or \
            (n_samples, n_classes)
        """
        if not self._loss.is_multiclass:
            if self.fit_intercept:
                intercept = coef[-1]
                w = coef[:-1]
            else:
                intercept = 0.0
                w = coef
        else:
            # reshape to (n_classes, n_dof)
            w = coef.reshape(self._loss.n_classes, -1)
            if self.fit_intercept:
                intercept = w[:, -1]
                w = w[:, :-1]
            else:
                intercept = 0.0

        raw_prediction = X @ w.T + intercept
        return w, intercept, raw_prediction

    def loss(self, coef, X, y, sample_weight=None, alpha=0.0, n_threads=1):
        """Compute the loss as sum over point-wise losses.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,) or (n_classes * n_dof,)
            Coefficients of a linear model.
        y : C/F-contiguous array of shape (n_samples,)
            Observed, true target values.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        sample_weight : None or C/F-contiguous array of shape (n_samples,)
            Sample weights.
        alpha: float
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.

        Returns
        -------
        loss : float
            Sum of losses per sample plus penalty.
        """
        w, intercept, raw_prediction = self._w_intercept_raw(coef, X)

        loss = self._loss.loss(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )
        loss = loss.sum()

        if w.ndim == 1:
            return loss + 0.5 * alpha * (w @ w)
        else:
            return loss + 0.5 * alpha * squared_norm(w)

    def loss_gradient(self, coef, X, y, sample_weight=None, alpha=0.0, n_threads=1):
        """Computes the sum/average of loss and gradient.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,) or (n_classes * n_dof,)
            Coefficients of a linear model.
        y : C/F-contiguous array of shape (n_samples,)
            Observed, true target values.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        sample_weight : None or C/F-contiguous array of shape (n_samples,)
            Sample weights.
        alpha: float
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.

        Returns
        -------
        loss : float
            Sum of losses per sample plus penalty.

        gradient : ndarray of shape (n_dof,) or (n_classes * n_dof)
             The gradient of the loss as ravelled array.
        """
        n_features, n_classes = X.shape[1], self._loss.n_classes
        n_dof = n_features + self.fit_intercept
        w, intercept, raw_prediction = self._w_intercept_raw(coef, X)

        loss, gradient = self._loss.loss_gradient(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )
        loss = loss.sum()

        if not self._loss.is_multiclass:
            loss += 0.5 * alpha * (w @ w)
            grad = np.empty_like(coef, dtype=X.dtype)
            grad[:n_features] = X.T @ gradient + alpha * w
            if self.fit_intercept:
                grad[-1] = gradient.sum()
            return loss, grad
        else:
            loss += 0.5 * alpha * squared_norm(w)
            grad = np.empty((n_classes, n_dof), dtype=X.dtype)
            # gradient.shape = (n_samples, n_classes)
            grad[:, :n_features] = gradient.T @ X + alpha * w
            if self.fit_intercept:
                grad[:, -1] = gradient.sum(axis=0)
            return loss, grad.ravel()

    def gradient(self, coef, X, y, sample_weight=None, alpha=0.0, n_threads=1):
        """Computes the gradient.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,) or (n_classes * n_dof,)
            Coefficients of a linear model.
        y : C/F-contiguous array of shape (n_samples,)
            Observed, true target values.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        sample_weight : None or C/F-contiguous array of shape (n_samples,)
            Sample weights.
        alpha: float
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.

        Returns
        -------
        gradient : ndarray of shape (n_dof,) or (n_classes * n_dof)
             The gradient of the loss as ravelled array.
        """
        n_features, n_classes = X.shape[1], self._loss.n_classes
        n_dof = n_features + self.fit_intercept
        w, intercept, raw_prediction = self._w_intercept_raw(coef, X)

        gradient = self._loss.gradient(
            y_true=y,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )

        if not self._loss.is_multiclass:
            grad = np.empty_like(coef, dtype=X.dtype)
            grad[:n_features] = X.T @ gradient + alpha * w
            if self.fit_intercept:
                grad[-1] = gradient.sum()
            return grad
        else:
            grad = np.empty((n_classes, n_dof), dtype=X.dtype)
            # gradient.shape = (n_samples, n_classes)
            grad[:, :n_features] = gradient.T @ X + alpha * w
            if self.fit_intercept:
                grad[:, -1] = gradient.sum(axis=0)
            return grad.ravel()

    def gradient_hessp(self, coef, X, y, sample_weight=None, alpha=0.0, n_threads=1):
        """Computes gradient and hessp (hessian product function).

        Parameters
        ----------
        coef : ndarray of shape (n_dof,) or (n_classes * n_dof,)
            Coefficients of a linear model.
        y : C/F-contiguous array of shape (n_samples,)
            Observed, true target values.
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.
        sample_weight : None or C/F-contiguous array of shape (n_samples,)
            Sample weights.
        alpha: float
            L2 regularization strength
        n_threads : int, default=1
            Number of OpenMP threads to use.

        Returns
        -------
        gradient : ndarray of shape (n_dof,) or (n_classes * n_dof)
             The gradient of the loss as ravelled array.

        hessp : callable
            Function that takes in a vector input of shape of gradient and
            and returns matrix-vector product with hessian.
        """
        (n_samples, n_features), n_classes = X.shape, self._loss.n_classes
        w, intercept, raw_prediction = self._w_intercept_raw(coef, X)

        if not self._loss.is_multiclass:
            gradient, hessian = self._loss.gradient_hessian(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad = np.empty_like(coef, dtype=X.dtype)
            grad[:n_features] = X.T @ gradient + alpha * w
            if self.fit_intercept:
                grad[-1] = gradient.sum()

            # Precompute as much as possible: hX, hh_intercept and hsum
            hsum = hessian.sum()
            if sparse.issparse(X):
                hX = sparse.dia_matrix((hessian, 0), shape=(n_samples, n_samples)) @ X
            else:
                hX = hessian[:, np.newaxis] * X

            if self.fit_intercept:
                # Calculate the double derivative with respect to intercept.
                # Note: In case hX is sparse, hX.sum is a matrix object.
                hh_intercept = np.squeeze(np.array(hX.sum(axis=0)))

            # With intercept included and alpha = 0, hessp returns
            # res = (X, 1)' @ diag(h) @ (X, 1) @ s
            #     = (X, 1)' @ (hX @ s[:n_features], sum(h) * s[-1])
            # res[:n_features] = X' @ hX @ s[:n_features] + sum(h) * s[-1]
            # res[:-1] = 1' @ hX @ s[:n_features] + sum(h) * s[-1]
            def hessp(s):
                ret = np.empty_like(s)
                ret[:n_features] = X.T @ (hX @ s[:n_features])
                ret[:n_features] += alpha * s[:n_features]

                if self.fit_intercept:
                    ret[:n_features] += s[-1] * hh_intercept
                    ret[-1] = hh_intercept @ s[:n_features] + hsum * s[-1]
                return ret

        else:
            # Here we may safely assume HalfMultinomialLoss aka categorical
            # cross-entropy.
            # HalfMultinomialLoss computes only the diagonal part of the hessian, i.e.
            # diagonal in the classes. Here, we want the matrix-vector product of the
            # full hessian. Therefore, we call gradient_proba.
            gradient, proba = self._loss.gradient_proba(
                y_true=y,
                raw_prediction=raw_prediction,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            grad = np.empty_like(coef.reshape(n_classes, -1), dtype=X.dtype)
            grad[:, :n_features] = gradient.T @ X + alpha * w
            if self.fit_intercept:
                grad[:, -1] = gradient.sum(axis=0)

            # Full hessian-vector product, i.e. not only the diagonal part of the
            # hessian. Derivation with some index battle for inupt vector s:
            #   - sample index i
            #   - feature indices j, m
            #   - class indices k, l
            #   - 1_{k=l} is one if k=l else 0
            #   - p_i_k is class probability of sample i and class k
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
                s = s.reshape(n_classes, -1)  # shape = (n_classes, n_dof)
                if self.fit_intercept:
                    s_intercept = s[:, -1]
                    s = s[:, :-1]
                else:
                    s_intercept = 0
                tmp = X @ s.T + s_intercept  # X_{im} * s_k_m
                tmp += (-proba * tmp).sum(axis=1)[:, np.newaxis]  # - sum_l ..
                tmp *= proba  # * p_i_k
                if sample_weight is not None:
                    tmp *= sample_weight[:, np.newaxis]
                hessProd = np.empty_like(grad)
                hessProd[:, :n_features] = tmp.T @ X + alpha * s
                if self.fit_intercept:
                    hessProd[:, -1] = tmp.sum(axis=0)
                return hessProd.ravel()

        return grad.ravel(), hessp
