"""
Logistic Regression
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numbers
import warnings
from numbers import Integral, Real

import numpy as np
from scipy import optimize

from sklearn._loss.loss import HalfBinomialLoss, HalfMultinomialLoss
from sklearn.base import _fit_context
from sklearn.linear_model._base import (
    BaseEstimator,
    LinearClassifierMixin,
    SparseCoefMixin,
)
from sklearn.linear_model._glm.glm import NewtonCholeskySolver
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.linear_model._sag import sag_solver
from sklearn.metrics import get_scorer, get_scorer_names
from sklearn.model_selection import check_cv
from sklearn.preprocessing import LabelEncoder
from sklearn.svm._base import _fit_liblinear
from sklearn.utils import (
    Bunch,
    check_array,
    check_consistent_length,
    check_random_state,
    compute_class_weight,
)
from sklearn.utils._param_validation import Hidden, Interval, StrOptions
from sklearn.utils.extmath import row_norms, softmax
from sklearn.utils.fixes import _get_additional_lbfgs_options_dict
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    _routing_enabled,
    process_routing,
)
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.optimize import _check_optimize_result, _newton_cg
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.validation import (
    _check_method_params,
    _check_sample_weight,
    check_is_fitted,
    validate_data,
)

_LOGISTIC_SOLVER_CONVERGENCE_MSG = (
    "Please also refer to the documentation for alternative solver options:\n"
    "    https://scikit-learn.org/stable/modules/linear_model.html"
    "#logistic-regression"
)


def _check_solver(solver, penalty, dual):
    if solver not in ["liblinear", "saga"] and penalty not in ("l2", None):
        raise ValueError(
            f"Solver {solver} supports only 'l2' or None penalties, got {penalty} "
            "penalty."
        )
    if solver != "liblinear" and dual:
        raise ValueError(f"Solver {solver} supports only dual=False, got dual={dual}")

    if penalty == "elasticnet" and solver != "saga":
        raise ValueError(
            f"Only 'saga' solver supports elasticnet penalty, got solver={solver}."
        )

    if solver == "liblinear" and penalty is None:
        # TODO(1.10): update message to remove "as well as penalty=None".
        raise ValueError(
            "C=np.inf as well as penalty=None is not supported for the liblinear solver"
        )

    return solver


def _logistic_regression_path(
    X,
    y,
    *,
    classes,
    Cs=10,
    fit_intercept=True,
    max_iter=100,
    tol=1e-4,
    verbose=0,
    solver="lbfgs",
    coef=None,
    class_weight=None,
    dual=False,
    penalty="l2",
    intercept_scaling=1.0,
    random_state=None,
    check_input=True,
    max_squared_sum=None,
    sample_weight=None,
    l1_ratio=None,
    n_threads=1,
):
    """Compute a Logistic Regression model for a list of regularization
    parameters.

    This is an implementation that uses the result of the previous model
    to speed up computations along the set of solutions, making it faster
    than sequentially calling LogisticRegression for the different parameters.
    Note that there will be no speedup with liblinear solver, since it does
    not handle warm-starting.

    Read more in the :ref:`User Guide <logistic_regression>`.

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input data.

    y : array-like of shape (n_samples,) or (n_samples, n_targets)
        Input data, target values.

    classes : ndarray
        A list of class labels known to the classifier.

    Cs : int or array-like of shape (n_cs,), default=10
        List of values for the regularization parameter or integer specifying
        the number of regularization parameters that should be used. In this
        case, the parameters will be chosen in a logarithmic scale between
        1e-4 and 1e4.

    fit_intercept : bool, default=True
        Whether to fit an intercept for the model. In this case the shape of
        the returned array is (n_cs, n_features + 1).

    max_iter : int, default=100
        Maximum number of iterations for the solver.

    tol : float, default=1e-4
        Stopping criterion. For the newton-cg and lbfgs solvers, the iteration
        will stop when ``max{|g_i | i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient.

    verbose : int, default=0
        For the liblinear and lbfgs solvers set verbose to any positive
        number for verbosity.

    solver : {'lbfgs', 'liblinear', 'newton-cg', 'newton-cholesky', 'sag', 'saga'}, \
            default='lbfgs'
        Numerical solver to use.

    coef : array-like of shape (n_classes, features + int(fit_intercept)) or \
            (1, n_features + int(fit_intercept)) or \
            (n_features + int(fit_intercept)), default=None
        Initialization value for coefficients of logistic regression.
        Useless for liblinear solver.

    class_weight : dict or 'balanced', default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

    dual : bool, default=False
        Dual or primal formulation. Dual formulation is only implemented for
        l2 penalty with liblinear solver. Prefer dual=False when
        n_samples > n_features.

    penalty : {'l1', 'l2', 'elasticnet'}, default='l2'
        Used to specify the norm used in the penalization. The 'newton-cg',
        'sag' and 'lbfgs' solvers support only l2 penalties. 'elasticnet' is
        only supported by the 'saga' solver.

    intercept_scaling : float, default=1.
        Useful only when the solver `liblinear` is used
        and `self.fit_intercept` is set to `True`. In this case, `x` becomes
        `[x, self.intercept_scaling]`,
        i.e. a "synthetic" feature with constant value equal to
        `intercept_scaling` is appended to the instance vector.
        The intercept becomes
        ``intercept_scaling * synthetic_feature_weight``.

        .. note::
            The synthetic feature weight is subject to L1 or L2
            regularization as all other features.
            To lessen the effect of regularization on synthetic feature weight
            (and therefore on the intercept) `intercept_scaling` has to be increased.

    random_state : int, RandomState instance, default=None
        Used when ``solver`` == 'sag', 'saga' or 'liblinear' to shuffle the
        data. See :term:`Glossary <random_state>` for details.

    check_input : bool, default=True
        If False, the input arrays X and y will not be checked.

    max_squared_sum : float, default=None
        Maximum squared sum of X over samples. Used only in SAG solver.
        If None, it will be computed, going through all the samples.
        The value should be precomputed to speed up cross validation.

    sample_weight : array-like of shape (n_samples,), default=None
        Array of weights that are assigned to individual samples.
        If not provided, then each sample is given unit weight.

    l1_ratio : float, default=None
        The Elastic-Net mixing parameter, with ``0 <= l1_ratio <= 1``. Only
        used if ``penalty='elasticnet'``. Setting ``l1_ratio=0`` is equivalent
        to using ``penalty='l2'``, while setting ``l1_ratio=1`` is equivalent
        to using ``penalty='l1'``. For ``0 < l1_ratio <1``, the penalty is a
        combination of L1 and L2.

    n_threads : int, default=1
       Number of OpenMP threads to use.

    Returns
    -------
    coefs : ndarray of shape (n_cs, n_classes, n_features + int(fit_intercept)) or \
            (n_cs, n_features + int(fit_intercept))
        List of coefficients for the Logistic Regression model. If fit_intercept is set
        to True, then the last dimension will be n_features + 1, where the last item
        represents the intercept.
        For binary problems the second dimension in n_classes is dropped, i.e. the shape
        will be `(n_cs, n_features + int(fit_intercept))`.

    Cs : ndarray
        Grid of Cs used for cross-validation.

    n_iter : array of shape (n_cs,)
        Actual number of iteration for each C in Cs.

    Notes
    -----
    You might get slightly different results with the solver liblinear than
    with the others since this uses LIBLINEAR which penalizes the intercept.

    .. versionchanged:: 0.19
        The "copy" parameter was removed.
    """
    if isinstance(Cs, numbers.Integral):
        Cs = np.logspace(-4, 4, Cs)

    solver = _check_solver(solver, penalty, dual)

    # Preprocessing.
    if check_input:
        X = check_array(
            X,
            accept_sparse="csr",
            dtype=np.float64,
            accept_large_sparse=solver not in ["liblinear", "sag", "saga"],
        )
        y = check_array(y, ensure_2d=False, dtype=None)
        check_consistent_length(X, y)

    if sample_weight is not None or class_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype, copy=True)

    n_samples, n_features = X.shape
    n_classes = len(classes)
    is_binary = n_classes == 2

    if solver == "liblinear" and not is_binary:
        raise ValueError(
            "The 'liblinear' solver does not support multiclass classification"
            " (n_classes >= 3). Either use another solver or wrap the "
            "estimator in a OneVsRestClassifier to keep applying a "
            "one-versus-rest scheme."
        )

    random_state = check_random_state(random_state)

    le = LabelEncoder().fit(classes)
    if class_weight is not None:
        class_weight_ = compute_class_weight(
            class_weight, classes=classes, y=y, sample_weight=sample_weight
        )
        sample_weight *= class_weight_[le.transform(y)]

    if is_binary:
        w0 = np.zeros(n_features + int(fit_intercept), dtype=X.dtype)
        mask = y == classes[1]
        y_bin = np.ones(y.shape, dtype=X.dtype)
        if solver == "liblinear":
            y_bin[~mask] = -1.0
        else:
            # HalfBinomialLoss, used for those solvers, represents y in [0, 1] instead
            # of in [-1, 1].
            y_bin[~mask] = 0.0
    else:
        # All solvers capable of a multinomial need LabelEncoder, not LabelBinarizer,
        # i.e. y as a 1d-array of integers. LabelEncoder also saves memory
        # compared to LabelBinarizer, especially when n_classes is large.
        Y_multi = le.transform(y).astype(X.dtype, copy=False)
        # It is important that w0 is F-contiguous.
        w0 = np.zeros(
            (classes.size, n_features + int(fit_intercept)), order="F", dtype=X.dtype
        )

    # IMPORTANT NOTE:
    # All solvers relying on LinearModelLoss need to scale the penalty with n_samples
    # or the sum of sample weights because the implemented logistic regression
    # objective here is (unfortunately)
    #     C * sum(pointwise_loss) + penalty
    # instead of (as LinearModelLoss does)
    #     mean(pointwise_loss) + 1/C * penalty
    if solver in ["lbfgs", "newton-cg", "newton-cholesky"]:
        # This needs to be calculated after sample_weight is multiplied by
        # class_weight. It is even tested that passing class_weight is equivalent to
        # passing sample_weights according to class_weight.
        sw_sum = n_samples if sample_weight is None else np.sum(sample_weight)

    if coef is not None:
        if is_binary:
            if coef.ndim == 1 and coef.shape[0] == n_features + int(fit_intercept):
                w0[:] = coef
            elif (
                coef.ndim == 2
                and coef.shape[0] == 1
                and coef.shape[1] == n_features + int(fit_intercept)
            ):
                w0[:] = coef[0]
            else:
                msg = (
                    f"Initialization coef is of shape {coef.shape}, expected shape "
                    f"{w0.shape} or (1, {w0.shape[0]})"
                )
                raise ValueError(msg)
        else:
            if (
                coef.ndim == 2
                and coef.shape[0] == n_classes
                and coef.shape[1] == n_features + int(fit_intercept)
            ):
                w0[:, : coef.shape[1]] = coef
            else:
                msg = (
                    f"Initialization coef is of shape {coef.shape}, expected shape "
                    f"{w0.shape}"
                )
                raise ValueError(msg)

    if is_binary:
        target = y_bin
        loss = LinearModelLoss(
            base_loss=HalfBinomialLoss(), fit_intercept=fit_intercept
        )
        if solver == "lbfgs":
            func = loss.loss_gradient
        elif solver == "newton-cg":
            func = loss.loss
            grad = loss.gradient
            hess = loss.gradient_hessian_product  # hess = [gradient, hessp]
        warm_start_sag = {"coef": np.expand_dims(w0, axis=1)}
    else:  # multinomial
        loss = LinearModelLoss(
            base_loss=HalfMultinomialLoss(n_classes=classes.size),
            fit_intercept=fit_intercept,
        )
        target = Y_multi
        if solver in ["lbfgs", "newton-cg", "newton-cholesky"]:
            # scipy.optimize.minimize and newton-cg accept only ravelled parameters,
            # i.e. 1d-arrays. LinearModelLoss expects classes to be contiguous and
            # reconstructs the 2d-array via w0.reshape((n_classes, -1), order="F").
            # As w0 is F-contiguous, ravel(order="F") also avoids a copy.
            w0 = w0.ravel(order="F")
        if solver == "lbfgs":
            func = loss.loss_gradient
        elif solver == "newton-cg":
            func = loss.loss
            grad = loss.gradient
            hess = loss.gradient_hessian_product  # hess = [gradient, hessp]
        warm_start_sag = {"coef": w0.T}

    coefs = list()
    n_iter = np.zeros(len(Cs), dtype=np.int32)
    for i, C in enumerate(Cs):
        if solver == "lbfgs":
            l2_reg_strength = 1.0 / (C * sw_sum)
            iprint = [-1, 50, 1, 100, 101][
                np.searchsorted(np.array([0, 1, 2, 3]), verbose)
            ]
            opt_res = optimize.minimize(
                func,
                w0,
                method="L-BFGS-B",
                jac=True,
                args=(X, target, sample_weight, l2_reg_strength, n_threads),
                options={
                    "maxiter": max_iter,
                    "maxls": 50,  # default is 20
                    "gtol": tol,
                    "ftol": 64 * np.finfo(float).eps,
                    **_get_additional_lbfgs_options_dict("iprint", iprint),
                },
            )
            n_iter_i = _check_optimize_result(
                solver,
                opt_res,
                max_iter,
                extra_warning_msg=_LOGISTIC_SOLVER_CONVERGENCE_MSG,
            )
            w0, loss = opt_res.x, opt_res.fun
        elif solver == "newton-cg":
            l2_reg_strength = 1.0 / (C * sw_sum)
            args = (X, target, sample_weight, l2_reg_strength, n_threads)
            w0, n_iter_i = _newton_cg(
                grad_hess=hess,
                func=func,
                grad=grad,
                x0=w0,
                args=args,
                maxiter=max_iter,
                tol=tol,
                verbose=verbose,
            )
        elif solver == "newton-cholesky":
            l2_reg_strength = 1.0 / (C * sw_sum)
            sol = NewtonCholeskySolver(
                coef=w0,
                linear_loss=loss,
                l2_reg_strength=l2_reg_strength,
                tol=tol,
                max_iter=max_iter,
                n_threads=n_threads,
                verbose=verbose,
            )
            w0 = sol.solve(X=X, y=target, sample_weight=sample_weight)
            n_iter_i = sol.iteration
        elif solver == "liblinear":
            coef_, intercept_, n_iter_i = _fit_liblinear(
                X,
                target,
                C,
                fit_intercept,
                intercept_scaling,
                None,
                penalty,
                dual,
                verbose,
                max_iter,
                tol,
                random_state,
                sample_weight=sample_weight,
            )
            if fit_intercept:
                w0 = np.concatenate([coef_.ravel(), intercept_])
            else:
                w0 = coef_.ravel()
            # n_iter_i is an array for each class. However, `target` is always encoded
            # in {-1, 1}, so we only take the first element of n_iter_i.
            n_iter_i = n_iter_i.item()

        elif solver in ["sag", "saga"]:
            if is_binary:
                loss = "log"
            else:
                target = target.astype(X.dtype, copy=False)
                loss = "multinomial"
            # alpha is for L2-norm, beta is for L1-norm
            if penalty == "l1":
                alpha = 0.0
                beta = 1.0 / C
            elif penalty == "l2":
                alpha = 1.0 / C
                beta = 0.0
            else:  # Elastic-Net penalty
                alpha = (1.0 / C) * (1 - l1_ratio)
                beta = (1.0 / C) * l1_ratio

            w0, n_iter_i, warm_start_sag = sag_solver(
                X,
                target,
                sample_weight,
                loss,
                alpha,
                beta,
                max_iter,
                tol,
                verbose,
                random_state,
                False,
                max_squared_sum,
                warm_start_sag,
                is_saga=(solver == "saga"),
            )

        else:
            msg = (
                "solver must be one of {'lbfgs', 'liblinear', 'newton-cg', "
                "'newton-cholesky', 'sag', 'saga'}, "
                f"got '{solver}' instead."
            )
            raise ValueError(msg)

        if is_binary:
            coefs.append(w0.copy())
        else:
            if solver in ["lbfgs", "newton-cg", "newton-cholesky"]:
                multi_w0 = np.reshape(w0, (n_classes, -1), order="F")
            else:
                multi_w0 = w0
            coefs.append(multi_w0.copy())

        n_iter[i] = n_iter_i

    return np.array(coefs), np.array(Cs), n_iter


# helper function for LogisticCV
def _log_reg_scoring_path(
    X,
    y,
    train,
    test,
    *,
    classes,
    Cs,
    scoring,
    fit_intercept,
    max_iter,
    tol,
    class_weight,
    verbose,
    solver,
    penalty,
    dual,
    intercept_scaling,
    random_state,
    max_squared_sum,
    sample_weight,
    l1_ratio,
    score_params,
):
    """Computes scores across logistic_regression_path

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Training data.

    y : array-like of shape (n_samples,) or (n_samples, n_targets)
        Target labels.

    train : list of indices
        The indices of the train set.

    test : list of indices
        The indices of the test set.

    classes : ndarray
        A list of class labels known to the classifier.

    Cs : int or list of floats
        Each of the values in Cs describes the inverse of
        regularization strength. If Cs is as an int, then a grid of Cs
        values are chosen in a logarithmic scale between 1e-4 and 1e4.

    scoring : str, callable or None
        The scoring method to use for cross-validation. Options:

        - str: see :ref:`scoring_string_names` for options.
        - callable: a scorer callable object (e.g., function) with signature
          ``scorer(estimator, X, y)``. See :ref:`scoring_callable` for details.
        - `None`: :ref:`accuracy <accuracy_score>` is used.

    fit_intercept : bool
        If False, then the bias term is set to zero. Else the last
        term of each coef_ gives us the intercept.

    max_iter : int
        Maximum number of iterations for the solver.

    tol : float
        Tolerance for stopping criteria.

    class_weight : dict or 'balanced'
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

    verbose : int
        For the liblinear and lbfgs solvers set verbose to any positive
        number for verbosity.

    solver : {'lbfgs', 'liblinear', 'newton-cg', 'newton-cholesky', 'sag', 'saga'}
        Decides which solver to use.

    penalty : {'l1', 'l2', 'elasticnet'}
        Used to specify the norm used in the penalization. The 'newton-cg',
        'sag' and 'lbfgs' solvers support only l2 penalties. 'elasticnet' is
        only supported by the 'saga' solver.

    dual : bool
        Dual or primal formulation. Dual formulation is only implemented for
        l2 penalty with liblinear solver. Prefer dual=False when
        n_samples > n_features.

    intercept_scaling : float
        Useful only when the solver `liblinear` is used
        and `self.fit_intercept` is set to `True`. In this case, `x` becomes
        `[x, self.intercept_scaling]`,
        i.e. a "synthetic" feature with constant value equal to
        `intercept_scaling` is appended to the instance vector.
        The intercept becomes
        ``intercept_scaling * synthetic_feature_weight``.

        .. note::
            The synthetic feature weight is subject to L1 or L2
            regularization as all other features.
            To lessen the effect of regularization on synthetic feature weight
            (and therefore on the intercept) `intercept_scaling` has to be increased.

    random_state : int, RandomState instance
        Used when ``solver`` == 'sag', 'saga' or 'liblinear' to shuffle the
        data. See :term:`Glossary <random_state>` for details.

    max_squared_sum : float
        Maximum squared sum of X over samples. Used only in SAG solver.
        If None, it will be computed, going through all the samples.
        The value should be precomputed to speed up cross validation.

    sample_weight : array-like of shape (n_samples,)
        Array of weights that are assigned to individual samples.
        If not provided, then each sample is given unit weight.

    l1_ratio : float
        The Elastic-Net mixing parameter, with ``0 <= l1_ratio <= 1``. Only
        used if ``penalty='elasticnet'``. Setting ``l1_ratio=0`` is equivalent
        to using ``penalty='l2'``, while setting ``l1_ratio=1`` is equivalent
        to using ``penalty='l1'``. For ``0 < l1_ratio <1``, the penalty is a
        combination of L1 and L2.

    score_params : dict
        Parameters to pass to the `score` method of the underlying scorer.

    Returns
    -------
    coefs : ndarray of shape (n_cs, n_classes, n_features + int(fit_intercept)) or \
            (n_cs, n_features + int(fit_intercept))
        List of coefficients for the Logistic Regression model. If fit_intercept is set
        to True, then the last dimension will be n_features + 1, where the last item
        represents the intercept.
        For binary problems the second dimension in n_classes is dropped, i.e. the shape
        will be `(n_cs, n_features + int(fit_intercept))`.

    Cs : ndarray of shape (n_cs,)
        Grid of Cs used for cross-validation.

    scores : ndarray of shape (n_cs,)
        Scores obtained for each Cs.

    n_iter : ndarray of shape (n_cs,)
        Actual number of iteration for each C in Cs.
    """
    X_train = X[train]
    X_test = X[test]
    y_train = y[train]
    y_test = y[test]

    sw_train, sw_test = None, None
    if sample_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, X)
        sw_train = sample_weight[train]
        sw_test = sample_weight[test]

    # Note: We pass classes for the whole dataset to avoid inconsistencies,
    # i.e. different number of classes in different folds. This way, if a class
    # is not present in a fold, _logistic_regression_path will still return
    # coefficients associated to this class.
    coefs, Cs, n_iter = _logistic_regression_path(
        X_train,
        y_train,
        classes=classes,
        Cs=Cs,
        l1_ratio=l1_ratio,
        fit_intercept=fit_intercept,
        solver=solver,
        max_iter=max_iter,
        class_weight=class_weight,
        tol=tol,
        verbose=verbose,
        dual=dual,
        penalty=penalty,
        intercept_scaling=intercept_scaling,
        random_state=random_state,
        check_input=False,
        max_squared_sum=max_squared_sum,
        sample_weight=sw_train,
    )

    log_reg = LogisticRegression(solver=solver)

    # The score method of Logistic Regression has a classes_ attribute.
    log_reg.classes_ = classes

    scores = list()

    scoring = get_scorer(scoring)
    for w in coefs:
        if fit_intercept:
            log_reg.coef_ = w[..., :-1]
            log_reg.intercept_ = w[..., -1]
        else:
            log_reg.coef_ = w
            log_reg.intercept_ = 0.0

        if scoring is None:
            scores.append(log_reg.score(X_test, y_test, sample_weight=sw_test))
        else:
            score_params = score_params or {}
            score_params = _check_method_params(X=X, params=score_params, indices=test)
            # FIXME: If scoring = "neg_brier_score" and if not all class labels
            # are present in y_test, the following fails. Maybe we can pass
            # "labels=classes" to the call of scoring.
            scores.append(scoring(log_reg, X_test, y_test, **score_params))
    return coefs, Cs, np.array(scores), n_iter


class LogisticRegression(LinearClassifierMixin, SparseCoefMixin, BaseEstimator):
    """
    Logistic Regression (aka logit, MaxEnt) classifier.

    This class implements regularized logistic regression using a set of available
    solvers. **Note that regularization is applied by default**. It can handle both
    dense and sparse input `X`. Use C-ordered arrays or CSR matrices containing 64-bit
    floats for optimal performance; any other input format will be converted (and
    copied).

    The solvers 'lbfgs', 'newton-cg', 'newton-cholesky' and 'sag' support only L2
    regularization with primal formulation, or no regularization. The 'liblinear'
    solver supports both L1 and L2 regularization (but not both, i.e. elastic-net),
    with a dual formulation only for the L2 penalty. The Elastic-Net (combination of L1
    and L2) regularization is only supported by the 'saga' solver.

    For :term:`multiclass` problems (whenever `n_classes >= 3`), all solvers except
    'liblinear' optimize the (penalized) multinomial loss. 'liblinear' only handles
    binary classification but can be extended to handle multiclass by using
    :class:`~sklearn.multiclass.OneVsRestClassifier`.

    Read more in the :ref:`User Guide <logistic_regression>`.

    Parameters
    ----------
    penalty : {'l1', 'l2', 'elasticnet', None}, default='l2'
        Specify the norm of the penalty:

        - `None`: no penalty is added;
        - `'l2'`: add a L2 penalty term and it is the default choice;
        - `'l1'`: add a L1 penalty term;
        - `'elasticnet'`: both L1 and L2 penalty terms are added.

        .. warning::
           Some penalties may not work with some solvers. See the parameter
           `solver` below, to know the compatibility between the penalty and
           solver.

        .. versionadded:: 0.19
           l1 penalty with SAGA solver (allowing 'multinomial' + L1)

        .. deprecated:: 1.8
           `penalty` was deprecated in version 1.8 and will be removed in 1.10.
           Use `l1_ratio` instead. `l1_ratio=0` for `penalty='l2'`, `l1_ratio=1` for
           `penalty='l1'` and `l1_ratio` set to any float between 0 and 1 for
           `'penalty='elasticnet'`.

    C : float, default=1.0
        Inverse of regularization strength; must be a positive float.
        Like in support vector machines, smaller values specify stronger
        regularization. `C=np.inf` results in unpenalized logistic regression.
        For a visual example on the effect of tuning the `C` parameter
        with an L1 penalty, see:
        :ref:`sphx_glr_auto_examples_linear_model_plot_logistic_path.py`.

    l1_ratio : float, default=0.0
        The Elastic-Net mixing parameter, with `0 <= l1_ratio <= 1`. Setting
        `l1_ratio=1` gives a pure L1-penalty, setting `l1_ratio=0` a pure L2-penalty.
        Any value between 0 and 1 gives an Elastic-Net penalty of the form
        `l1_ratio * L1 + (1 - l1_ratio) * L2`.

        .. warning::
           Certain values of `l1_ratio`, i.e. some penalties, may not work with some
           solvers. See the parameter `solver` below, to know the compatibility between
           the penalty and solver.

        .. versionchanged:: 1.8
            Default value changed from None to 0.0.

        .. deprecated:: 1.8
            `None` is deprecated and will be removed in version 1.10. Always use
            `l1_ratio` to specify the penalty type.

    dual : bool, default=False
        Dual (constrained) or primal (regularized, see also
        :ref:`this equation <regularized-logistic-loss>`) formulation. Dual formulation
        is only implemented for l2 penalty with liblinear solver. Prefer `dual=False`
        when n_samples > n_features.

    tol : float, default=1e-4
        Tolerance for stopping criteria.

    fit_intercept : bool, default=True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the decision function.

    intercept_scaling : float, default=1
        Useful only when the solver `liblinear` is used
        and `self.fit_intercept` is set to `True`. In this case, `x` becomes
        `[x, self.intercept_scaling]`,
        i.e. a "synthetic" feature with constant value equal to
        `intercept_scaling` is appended to the instance vector.
        The intercept becomes
        ``intercept_scaling * synthetic_feature_weight``.

        .. note::
            The synthetic feature weight is subject to L1 or L2
            regularization as all other features.
            To lessen the effect of regularization on synthetic feature weight
            (and therefore on the intercept) `intercept_scaling` has to be increased.

    class_weight : dict or 'balanced', default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

        .. versionadded:: 0.17
           *class_weight='balanced'*

    random_state : int, RandomState instance, default=None
        Used when ``solver`` == 'sag', 'saga' or 'liblinear' to shuffle the
        data. See :term:`Glossary <random_state>` for details.

    solver : {'lbfgs', 'liblinear', 'newton-cg', 'newton-cholesky', 'sag', 'saga'}, \
            default='lbfgs'

        Algorithm to use in the optimization problem. Default is 'lbfgs'.
        To choose a solver, you might want to consider the following aspects:

        - 'lbfgs' is a good default solver because it works reasonably well for a wide
          class of problems.
        - For :term:`multiclass` problems (`n_classes >= 3`), all solvers except
          'liblinear' minimize the full multinomial loss, 'liblinear' will raise an
          error.
        - 'newton-cholesky' is a good choice for
          `n_samples` >> `n_features * n_classes`, especially with one-hot encoded
          categorical features with rare categories. Be aware that the memory usage
          of this solver has a quadratic dependency on `n_features * n_classes`
          because it explicitly computes the full Hessian matrix.
        - For small datasets, 'liblinear' is a good choice, whereas 'sag'
          and 'saga' are faster for large ones;
        - 'liblinear' can only handle binary classification by default. To apply a
          one-versus-rest scheme for the multiclass setting one can wrap it with the
          :class:`~sklearn.multiclass.OneVsRestClassifier`.

        .. warning::
           The choice of the algorithm depends on the penalty chosen (`l1_ratio=0`
           for L2-penalty, `l1_ratio=1` for L1-penalty and `0 < l1_ratio < 1` for
           Elastic-Net) and on (multinomial) multiclass support:

           ================= ======================== ======================
           solver            l1_ratio                 multinomial multiclass
           ================= ======================== ======================
           'lbfgs'           l1_ratio=0               yes
           'liblinear'       l1_ratio=1 or l1_ratio=0 no
           'newton-cg'       l1_ratio=0               yes
           'newton-cholesky' l1_ratio=0               yes
           'sag'             l1_ratio=0               yes
           'saga'            0<=l1_ratio<=1           yes
           ================= ======================== ======================

        .. note::
           'sag' and 'saga' fast convergence is only guaranteed on features
           with approximately the same scale. You can preprocess the data with
           a scaler from :mod:`sklearn.preprocessing`.

        .. seealso::
           Refer to the :ref:`User Guide <Logistic_regression>` for more
           information regarding :class:`LogisticRegression` and more specifically the
           :ref:`Table <logistic_regression_solvers>`
           summarizing solver/penalty supports.

        .. versionadded:: 0.17
           Stochastic Average Gradient (SAG) descent solver. Multinomial support in
           version 0.18.
        .. versionadded:: 0.19
           SAGA solver.
        .. versionchanged:: 0.22
           The default solver changed from 'liblinear' to 'lbfgs' in 0.22.
        .. versionadded:: 1.2
           newton-cholesky solver. Multinomial support in version 1.6.

    max_iter : int, default=100
        Maximum number of iterations taken for the solvers to converge.

    verbose : int, default=0
        For the liblinear and lbfgs solvers set verbose to any positive
        number for verbosity.

    warm_start : bool, default=False
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        Useless for liblinear solver. See :term:`the Glossary <warm_start>`.

        .. versionadded:: 0.17
           *warm_start* to support *lbfgs*, *newton-cg*, *sag*, *saga* solvers.

    n_jobs : int, default=None
        Does not have any effect.

        .. deprecated:: 1.8
           `n_jobs` is deprecated in version 1.8 and will be removed in 1.10.

    Attributes
    ----------

    classes_ : ndarray of shape (n_classes, )
        A list of class labels known to the classifier.

    coef_ : ndarray of shape (1, n_features) or (n_classes, n_features)
        Coefficient of the features in the decision function.

        `coef_` is of shape (1, n_features) when the given problem is binary.

    intercept_ : ndarray of shape (1,) or (n_classes,)
        Intercept (a.k.a. bias) added to the decision function.

        If `fit_intercept` is set to False, the intercept is set to zero.
        `intercept_` is of shape (1,) when the given problem is binary.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_iter_ : ndarray of shape (1, )
        Actual number of iterations for all classes.

        .. versionchanged:: 0.20

            In SciPy <= 1.0.0 the number of lbfgs iterations may exceed
            ``max_iter``. ``n_iter_`` will now report at most ``max_iter``.

    See Also
    --------
    SGDClassifier : Incrementally trained logistic regression (when given
        the parameter ``loss="log_loss"``).
    LogisticRegressionCV : Logistic regression with built-in cross validation.

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    Predict output may not match that of standalone liblinear in certain
    cases. See :ref:`differences from liblinear <liblinear_differences>`
    in the narrative documentation.

    References
    ----------

    L-BFGS-B -- Software for Large-scale Bound-constrained Optimization
        Ciyou Zhu, Richard Byrd, Jorge Nocedal and Jose Luis Morales.
        http://users.iems.northwestern.edu/~nocedal/lbfgsb.html

    LIBLINEAR -- A Library for Large Linear Classification
        https://www.csie.ntu.edu.tw/~cjlin/liblinear/

    SAG -- Mark Schmidt, Nicolas Le Roux, and Francis Bach
        Minimizing Finite Sums with the Stochastic Average Gradient
        https://hal.inria.fr/hal-00860051/document

    SAGA -- Defazio, A., Bach F. & Lacoste-Julien S. (2014).
            :arxiv:`"SAGA: A Fast Incremental Gradient Method With Support
            for Non-Strongly Convex Composite Objectives" <1407.0202>`

    Hsiang-Fu Yu, Fang-Lan Huang, Chih-Jen Lin (2011). Dual coordinate descent
        methods for logistic regression and maximum entropy models.
        Machine Learning 85(1-2):41-75.
        https://www.csie.ntu.edu.tw/~cjlin/papers/maxent_dual.pdf

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.linear_model import LogisticRegression
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = LogisticRegression(random_state=0).fit(X, y)
    >>> clf.predict(X[:2, :])
    array([0, 0])
    >>> clf.predict_proba(X[:2, :])
    array([[9.82e-01, 1.82e-02, 1.44e-08],
           [9.72e-01, 2.82e-02, 3.02e-08]])
    >>> clf.score(X, y)
    0.97

    For a comparison of the LogisticRegression with other classifiers see:
    :ref:`sphx_glr_auto_examples_classification_plot_classification_probability.py`.
    """

    _parameter_constraints: dict = {
        "penalty": [
            StrOptions({"l1", "l2", "elasticnet"}),
            None,
            Hidden(StrOptions({"deprecated"})),
        ],
        "C": [Interval(Real, 0, None, closed="right")],
        "l1_ratio": [Interval(Real, 0, 1, closed="both"), None],
        "dual": ["boolean"],
        "tol": [Interval(Real, 0, None, closed="left")],
        "fit_intercept": ["boolean"],
        "intercept_scaling": [Interval(Real, 0, None, closed="neither")],
        "class_weight": [dict, StrOptions({"balanced"}), None],
        "random_state": ["random_state"],
        "solver": [
            StrOptions(
                {"lbfgs", "liblinear", "newton-cg", "newton-cholesky", "sag", "saga"}
            )
        ],
        "max_iter": [Interval(Integral, 0, None, closed="left")],
        "verbose": ["verbose"],
        "warm_start": ["boolean"],
        "n_jobs": [None, Integral],
    }

    def __init__(
        self,
        penalty="deprecated",
        *,
        C=1.0,
        l1_ratio=0.0,
        dual=False,
        tol=1e-4,
        fit_intercept=True,
        intercept_scaling=1,
        class_weight=None,
        random_state=None,
        solver="lbfgs",
        max_iter=100,
        verbose=0,
        warm_start=False,
        n_jobs=None,
    ):
        self.penalty = penalty
        self.C = C
        self.l1_ratio = l1_ratio
        self.dual = dual
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.class_weight = class_weight
        self.random_state = random_state
        self.solver = solver
        self.max_iter = max_iter
        self.verbose = verbose
        self.warm_start = warm_start
        self.n_jobs = n_jobs

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """
        Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like of shape (n_samples,) default=None
            Array of weights that are assigned to individual samples.
            If not provided, then each sample is given unit weight.

            .. versionadded:: 0.17
               *sample_weight* support to LogisticRegression.

        Returns
        -------
        self
            Fitted estimator.

        Notes
        -----
        The SAGA solver supports both float64 and float32 bit arrays.
        """
        if self.penalty == "deprecated":
            if self.l1_ratio == 0 or self.l1_ratio is None:
                penalty = "l2"
                if self.l1_ratio is None:
                    warnings.warn(
                        (
                            "'l1_ratio=None' was deprecated in version 1.8 and will "
                            "trigger an error in 1.10. Use 0<=l1_ratio<=1 instead."
                        ),
                        FutureWarning,
                    )
            elif self.l1_ratio == 1:
                penalty = "l1"
            else:
                penalty = "elasticnet"
            if self.C == np.inf:
                penalty = None
        else:
            penalty = self.penalty
            warnings.warn(
                (
                    "'penalty' was deprecated in version 1.8 and will be removed in"
                    " 1.10. To avoid this warning, leave 'penalty' set to its default"
                    " value and use 'l1_ratio' or 'C' instead."
                    " Use l1_ratio=0 instead of penalty='l2',"
                    " l1_ratio=1 instead of penalty='l1', and "
                    "C=np.inf instead of penalty=None."
                ),
                FutureWarning,
            )

        solver = _check_solver(self.solver, penalty, self.dual)

        if penalty != "elasticnet" and (
            self.l1_ratio is not None and 0 < self.l1_ratio < 1
        ):
            warnings.warn(
                "l1_ratio parameter is only used when penalty is "
                "'elasticnet'. Got "
                "(penalty={})".format(penalty)
            )
        if (self.penalty == "l2" and self.l1_ratio != 0) or (
            self.penalty == "l1" and self.l1_ratio != 1
        ):
            warnings.warn(
                f"Inconsistent values: penalty={self.penalty} with "
                f"l1_ratio={self.l1_ratio}. penalty is deprecated. Please use "
                f"l1_ratio only."
            )
        if penalty == "elasticnet" and self.l1_ratio is None:
            raise ValueError("l1_ratio must be specified when penalty is elasticnet.")

        if penalty is None:
            if self.C != 1.0:  # default values
                warnings.warn(
                    "Setting penalty=None will ignore the C and l1_ratio parameters"
                )
                # Note that check for l1_ratio is done right above
            C_ = np.inf
            penalty = "l2"
        else:
            C_ = self.C

        msg = (
            "'n_jobs' has no effect since 1.8 and will be removed in 1.10. "
            f"You provided 'n_jobs={self.n_jobs}', please leave it unspecified."
        )
        if self.n_jobs is not None:
            warnings.warn(msg, category=FutureWarning)

        if solver == "lbfgs":
            _dtype = np.float64
        else:
            _dtype = [np.float64, np.float32]

        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse="csr",
            dtype=_dtype,
            order="C",
            accept_large_sparse=solver not in ["liblinear", "sag", "saga"],
        )
        n_features = X.shape[1]
        check_classification_targets(y)
        self.classes_ = np.unique(y)
        n_classes = len(self.classes_)
        is_binary = n_classes == 2

        if solver == "liblinear":
            if not is_binary:
                raise ValueError(
                    "The 'liblinear' solver does not support multiclass classification"
                    " (n_classes >= 3). Either use another solver or wrap the "
                    "estimator in a OneVsRestClassifier to keep applying a "
                    "one-versus-rest scheme."
                )
            if np.max(X) > 1e30:
                raise ValueError(
                    "Using the 'liblinear' solver while X contains a maximum "
                    "value > 1e30 results in a frozen fit. Please choose another "
                    "solver or rescale the input X."
                )
            self.coef_, self.intercept_, self.n_iter_ = _fit_liblinear(
                X,
                y,
                self.C,
                self.fit_intercept,
                self.intercept_scaling,
                self.class_weight,
                penalty,
                self.dual,
                self.verbose,
                self.max_iter,
                self.tol,
                self.random_state,
                sample_weight=sample_weight,
            )
            return self

        if solver in ["sag", "saga"]:
            max_squared_sum = row_norms(X, squared=True).max()
        else:
            max_squared_sum = None

        if n_classes < 2:
            raise ValueError(
                "This solver needs samples of at least 2 classes"
                " in the data, but the data contains only one"
                " class: %r" % self.classes_[0]
            )

        if self.warm_start:
            warm_start_coef = getattr(self, "coef_", None)
        else:
            warm_start_coef = None
        if warm_start_coef is not None and self.fit_intercept:
            warm_start_coef = np.append(
                warm_start_coef, self.intercept_[:, np.newaxis], axis=1
            )

        # TODO: enable multi-threading if benchmarks show a positive effect,
        # see https://github.com/scikit-learn/scikit-learn/issues/32162
        n_threads = 1

        coefs, _, n_iter = _logistic_regression_path(
            X,
            y,
            classes=self.classes_,
            Cs=[C_],
            l1_ratio=self.l1_ratio,
            fit_intercept=self.fit_intercept,
            tol=self.tol,
            verbose=self.verbose,
            solver=solver,
            max_iter=self.max_iter,
            class_weight=self.class_weight,
            check_input=False,
            random_state=self.random_state,
            coef=warm_start_coef,
            penalty=penalty,
            max_squared_sum=max_squared_sum,
            sample_weight=sample_weight,
            n_threads=n_threads,
        )

        self.n_iter_ = np.asarray(n_iter, dtype=np.int32)

        self.coef_ = coefs[0]
        if self.fit_intercept:
            if is_binary:
                self.intercept_ = self.coef_[-1:]
                self.coef_ = self.coef_[:-1][None, :]
            else:
                self.intercept_ = self.coef_[:, -1]
                self.coef_ = self.coef_[:, :-1]
        else:
            if is_binary:
                self.intercept_ = np.zeros(1, dtype=X.dtype)
                self.coef_ = self.coef_[None, :]
            else:
                self.intercept_ = np.zeros(n_classes, dtype=X.dtype)

        return self

    def predict_proba(self, X):
        """
        Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        For a multiclass / multinomial problem the softmax function is used to find
        the predicted probability of each class.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Vector to be scored, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        T : array-like of shape (n_samples, n_classes)
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in ``self.classes_``.
        """
        check_is_fitted(self)

        is_binary = self.classes_.size <= 2
        if is_binary:
            return super()._predict_proba_lr(X)
        else:
            decision_2d = self.decision_function(X)
            return softmax(decision_2d, copy=False)

    def predict_log_proba(self, X):
        """
        Predict logarithm of probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Vector to be scored, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        T : array-like of shape (n_samples, n_classes)
            Returns the log-probability of the sample for each class in the
            model, where classes are ordered as they are in ``self.classes_``.
        """
        return np.log(self.predict_proba(X))

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = True
        if self.solver == "liblinear":
            tags.classifier_tags.multi_class = False

        return tags


class LogisticRegressionCV(LogisticRegression, LinearClassifierMixin, BaseEstimator):
    """Logistic Regression CV (aka logit, MaxEnt) classifier.

    See glossary entry for :term:`cross-validation estimator`.

    This class implements regularized logistic regression with implicit cross
    validation for the penalty parameters `C` and `l1_ratio`, see
    :class:`LogisticRegression`, using a set of available solvers.

    The solvers 'lbfgs', 'newton-cg', 'newton-cholesky' and 'sag' support only L2
    regularization with primal formulation. The 'liblinear'
    solver supports both L1 and L2 regularization (but not both, i.e. elastic-net),
    with a dual formulation only for the L2 penalty. The Elastic-Net (combination of L1
    and L2) regularization is only supported by the 'saga' solver.

    For the grid of `Cs` values and `l1_ratios` values, the best hyperparameter
    is selected by the cross-validator
    :class:`~sklearn.model_selection.StratifiedKFold`, but it can be changed
    using the :term:`cv` parameter. All solvers except 'liblinear' can warm-start the
    coefficients (see :term:`Glossary<warm_start>`).

    Read more in the :ref:`User Guide <logistic_regression>`.

    Parameters
    ----------
    Cs : int or list of floats, default=10
        Each of the values in Cs describes the inverse of regularization
        strength. If Cs is as an int, then a grid of Cs values are chosen
        in a logarithmic scale between 1e-4 and 1e4.
        Like in support vector machines, smaller values specify stronger
        regularization.

    l1_ratios : array-like of shape (n_l1_ratios), default=None
        Floats between 0 and 1 passed as Elastic-Net mixing parameter (scaling between
        L1 and L2 penalties). For `l1_ratio = 0` the penalty is an L2 penalty. For
        `l1_ratio = 1` it is an L1 penalty. For `0 < l1_ratio < 1`, the penalty is a
        combination of L1 and L2.
        All the values of the given array-like are tested by cross-validation and the
        one giving the best prediction score is used.

        .. warning::
           Certain values of `l1_ratios`, i.e. some penalties, may not work with some
           solvers. See the parameter `solver` below, to know the compatibility between
           the penalty and solver.

        .. deprecated:: 1.8
            `l1_ratios=None` is deprecated in 1.8 and will raise an error
            in version 1.10. Default value will change from `None` to `(0.0,)`
            in version 1.10.

    fit_intercept : bool, default=True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the decision function.

    cv : int or cross-validation generator, default=None
        The default cross-validation generator used is Stratified K-Folds.
        If an integer is provided, it specifies the number of folds, `n_folds`, used.
        See the module :mod:`sklearn.model_selection` module for the
        list of possible cross-validation objects.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    dual : bool, default=False
        Dual (constrained) or primal (regularized, see also
        :ref:`this equation <regularized-logistic-loss>`) formulation. Dual formulation
        is only implemented for l2 penalty with liblinear solver. Prefer dual=False when
        n_samples > n_features.

    penalty : {'l1', 'l2', 'elasticnet'}, default='l2'
        Specify the norm of the penalty:

        - `'l2'`: add a L2 penalty term (used by default);
        - `'l1'`: add a L1 penalty term;
        - `'elasticnet'`: both L1 and L2 penalty terms are added.

        .. warning::
           Some penalties may not work with some solvers. See the parameter
           `solver` below, to know the compatibility between the penalty and
           solver.

        .. deprecated:: 1.8
           `penalty` was deprecated in version 1.8 and will be removed in 1.10.
           Use `l1_ratio` instead. `l1_ratio=0` for `penalty='l2'`, `l1_ratio=1` for
           `penalty='l1'` and `l1_ratio` set to any float between 0 and 1 for
           `'penalty='elasticnet'`.

    scoring : str or callable, default=None
        The scoring method to use for cross-validation. Options:

        - str: see :ref:`scoring_string_names` for options.
        - callable: a scorer callable object (e.g., function) with signature
          ``scorer(estimator, X, y)``. See :ref:`scoring_callable` for details.
        - `None`: :ref:`accuracy <accuracy_score>` is used.

    solver : {'lbfgs', 'liblinear', 'newton-cg', 'newton-cholesky', 'sag', 'saga'}, \
            default='lbfgs'

        Algorithm to use in the optimization problem. Default is 'lbfgs'.
        To choose a solver, you might want to consider the following aspects:

        - 'lbfgs' is a good default solver because it works reasonably well for a wide
          class of problems.
        - For :term:`multiclass` problems (`n_classes >= 3`), all solvers except
          'liblinear' minimize the full multinomial loss, 'liblinear' will raise an
          error.
        - 'newton-cholesky' is a good choice for
          `n_samples` >> `n_features * n_classes`, especially with one-hot encoded
          categorical features with rare categories. Be aware that the memory usage
          of this solver has a quadratic dependency on `n_features * n_classes`
          because it explicitly computes the full Hessian matrix.
        - For small datasets, 'liblinear' is a good choice, whereas 'sag'
          and 'saga' are faster for large ones;
        - 'liblinear' might be slower in :class:`LogisticRegressionCV`
          because it does not handle warm-starting.
        - 'liblinear' can only handle binary classification by default. To apply a
          one-versus-rest scheme for the multiclass setting one can wrap it with the
          :class:`~sklearn.multiclass.OneVsRestClassifier`.

        .. warning::
           The choice of the algorithm depends on the penalty (`l1_ratio=0` for
           L2-penalty, `l1_ratio=1` for L1-penalty and `0 < l1_ratio < 1` for
           Elastic-Net) chosen and on (multinomial) multiclass support:

           ================= ======================== ======================
           solver            l1_ratio                 multinomial multiclass
           ================= ======================== ======================
           'lbfgs'           l1_ratio=0               yes
           'liblinear'       l1_ratio=1 or l1_ratio=0 no
           'newton-cg'       l1_ratio=0               yes
           'newton-cholesky' l1_ratio=0               yes
           'sag'             l1_ratio=0               yes
           'saga'            0<=l1_ratio<=1           yes
           ================= ======================== ======================

        .. note::
           'sag' and 'saga' fast convergence is only guaranteed on features
           with approximately the same scale. You can preprocess the data with
           a scaler from :mod:`sklearn.preprocessing`.

        .. versionadded:: 0.17
           Stochastic Average Gradient (SAG) descent solver. Multinomial support in
           version 0.18.
        .. versionadded:: 0.19
           SAGA solver.
        .. versionadded:: 1.2
           newton-cholesky solver. Multinomial support in version 1.6.

    tol : float, default=1e-4
        Tolerance for stopping criteria.

    max_iter : int, default=100
        Maximum number of iterations of the optimization algorithm.

    class_weight : dict or 'balanced', default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

        .. versionadded:: 0.17
           class_weight == 'balanced'

    n_jobs : int, default=None
        Number of CPU cores used during the cross-validation loop.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : int, default=0
        For the 'liblinear', 'sag' and 'lbfgs' solvers set verbose to any
        positive number for verbosity.

    refit : bool, default=True
        If set to True, the scores are averaged across all folds, and the
        coefs and the C that corresponds to the best score is taken, and a
        final refit is done using these parameters.
        Otherwise the coefs, intercepts and C that correspond to the
        best scores across folds are averaged.

    intercept_scaling : float, default=1
        Useful only when the solver `liblinear` is used
        and `self.fit_intercept` is set to `True`. In this case, `x` becomes
        `[x, self.intercept_scaling]`,
        i.e. a "synthetic" feature with constant value equal to
        `intercept_scaling` is appended to the instance vector.
        The intercept becomes
        ``intercept_scaling * synthetic_feature_weight``.

        .. note::
            The synthetic feature weight is subject to L1 or L2
            regularization as all other features.
            To lessen the effect of regularization on synthetic feature weight
            (and therefore on the intercept) `intercept_scaling` has to be increased.

    random_state : int, RandomState instance, default=None
        Used when `solver='sag'`, 'saga' or 'liblinear' to shuffle the data.
        Note that this only applies to the solver and not the cross-validation
        generator. See :term:`Glossary <random_state>` for details.

    use_legacy_attributes : bool, default=True
        If True, use legacy values for attributes:

        - `C_` is an ndarray of shape (n_classes,) with the same value repeated
        - `l1_ratio_` is an ndarray of shape (n_classes,) with the same value repeated
        - `coefs_paths_` is a dict with class labels as keys and ndarrays as values
        - `scores_` is a dict with class labels as keys and ndarrays as values
        - `n_iter_` is an ndarray of shape (1, n_folds, n_cs) or similar

        If False, use new values for attributes:

        - `C_` is a float
        - `l1_ratio_` is a float
        - `coefs_paths_` is an ndarray of shape
          (n_folds, n_l1_ratios, n_cs, n_classes, n_features)
          For binary problems (n_classes=2), the 2nd last dimension is 1.
        - `scores_` is an ndarray of shape (n_folds, n_l1_ratios, n_cs)
        - `n_iter_` is an ndarray of shape (n_folds, n_l1_ratios, n_cs)

        .. versionchanged:: 1.10
           The default will change from True to False in version 1.10.
        .. deprecated:: 1.10
           `use_legacy_attributes` will be deprecated in version 1.10 and be removed in
           1.12.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes, )
        A list of class labels known to the classifier.

    coef_ : ndarray of shape (1, n_features) or (n_classes, n_features)
        Coefficient of the features in the decision function.

        `coef_` is of shape (1, n_features) when the given problem
        is binary.

    intercept_ : ndarray of shape (1,) or (n_classes,)
        Intercept (a.k.a. bias) added to the decision function.

        If `fit_intercept` is set to False, the intercept is set to zero.
        `intercept_` is of shape (1,) when the problem is binary.

    Cs_ : ndarray of shape (n_cs)
        Array of C i.e. inverse of regularization parameter values used
        for cross-validation.

    l1_ratios_ : ndarray of shape (n_l1_ratios)
        Array of l1_ratios used for cross-validation. If l1_ratios=None is used
        (i.e. penalty is not 'elasticnet'), this is set to ``[None]``

    coefs_paths_ : dict of ndarray of shape (n_folds, n_cs, n_dof) or \
            (n_folds, n_cs, n_l1_ratios, n_dof)
        A dict with classes as the keys, and the path of coefficients obtained
        during cross-validating across each fold (`n_folds`) and then across each Cs
        (`n_cs`).
        The size of the coefficients is the number of degrees of freedom (`n_dof`),
        i.e. without intercept `n_dof=n_features` and with intercept
        `n_dof=n_features+1`.
        If `penalty='elasticnet'`, there is an additional dimension for the number of
        l1_ratio values (`n_l1_ratios`), which gives a shape of
        ``(n_folds, n_cs, n_l1_ratios_, n_dof)``.
        See also parameter `use_legacy_attributes`.

    scores_ : dict
        A dict with classes as the keys, and the values as the
        grid of scores obtained during cross-validating each fold.
        The same score is repeated across all classes. Each dict value
        has shape ``(n_folds, n_cs)`` or ``(n_folds, n_cs, n_l1_ratios)`` if
        ``penalty='elasticnet'``.
        See also parameter `use_legacy_attributes`.

    C_ : ndarray of shape (n_classes,) or (1,)
        The value of C that maps to the best score, repeated n_classes times.
        If refit is set to False, the best C is the average of the
        C's that correspond to the best score for each fold.
        `C_` is of shape (1,) when the problem is binary.
        See also parameter `use_legacy_attributes`.

    l1_ratio_ : ndarray of shape (n_classes,) or (n_classes - 1,)
        The value of l1_ratio that maps to the best score, repeated n_classes times.
        If refit is set to False, the best l1_ratio is the average of the
        l1_ratio's that correspond to the best score for each fold.
        `l1_ratio_` is of shape (1,) when the problem is binary.
        See also parameter `use_legacy_attributes`.

    n_iter_ : ndarray of shape (1, n_folds, n_cs) or (1, n_folds, n_cs, n_l1_ratios)
        Actual number of iterations for all classes, folds and Cs.
        If `penalty='elasticnet'`, the shape is `(1, n_folds, n_cs, n_l1_ratios)`.
        See also parameter `use_legacy_attributes`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    LogisticRegression : Logistic regression without tuning the
        hyperparameter `C`.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.linear_model import LogisticRegressionCV
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = LogisticRegressionCV(
    ...     cv=5, random_state=0, use_legacy_attributes=False, l1_ratios=(0,)
    ... ).fit(X, y)
    >>> clf.predict(X[:2, :])
    array([0, 0])
    >>> clf.predict_proba(X[:2, :]).shape
    (2, 3)
    >>> clf.score(X, y)
    0.98...
    """

    _parameter_constraints: dict = {**LogisticRegression._parameter_constraints}

    for param in ["C", "warm_start", "l1_ratio"]:
        _parameter_constraints.pop(param)

    _parameter_constraints.update(
        {
            "Cs": [Interval(Integral, 1, None, closed="left"), "array-like"],
            "l1_ratios": ["array-like", None, Hidden(StrOptions({"warn"}))],
            "cv": ["cv_object"],
            "scoring": [StrOptions(set(get_scorer_names())), callable, None],
            "refit": ["boolean"],
            "penalty": [
                StrOptions({"l1", "l2", "elasticnet"}),
                Hidden(StrOptions({"deprecated"})),
            ],
            "use_legacy_attributes": ["boolean", Hidden(StrOptions({"warn"}))],
        }
    )

    def __init__(
        self,
        *,
        Cs=10,
        l1_ratios="warn",
        fit_intercept=True,
        cv=None,
        dual=False,
        penalty="deprecated",
        scoring=None,
        solver="lbfgs",
        tol=1e-4,
        max_iter=100,
        class_weight=None,
        n_jobs=None,
        verbose=0,
        refit=True,
        intercept_scaling=1.0,
        random_state=None,
        use_legacy_attributes="warn",
    ):
        self.Cs = Cs
        self.l1_ratios = l1_ratios
        self.fit_intercept = fit_intercept
        self.cv = cv
        self.dual = dual
        self.penalty = penalty
        self.scoring = scoring
        self.tol = tol
        self.max_iter = max_iter
        self.class_weight = class_weight
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.solver = solver
        self.refit = refit
        self.intercept_scaling = intercept_scaling
        self.random_state = random_state
        self.use_legacy_attributes = use_legacy_attributes

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None, **params):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like of shape (n_samples,) default=None
            Array of weights that are assigned to individual samples.
            If not provided, then each sample is given unit weight.

        **params : dict
            Parameters to pass to the underlying splitter and scorer.

            .. versionadded:: 1.4

        Returns
        -------
        self : object
            Fitted LogisticRegressionCV estimator.
        """
        _raise_for_params(params, self, "fit")

        if isinstance(self.l1_ratios, str) and self.l1_ratios == "warn":
            l1_ratios = None
            warnings.warn(
                (
                    "The default value for l1_ratios will change from None to (0.0,) "
                    "in version 1.10. From version 1.10 onwards, only array-like "
                    "with values in [0, 1] will be allowed, None will be forbidden. "
                    "To avoid this warning, explicitly set a value, "
                    "e.g. l1_ratios=(0,)."
                ),
                FutureWarning,
            )
        else:
            l1_ratios = self.l1_ratios

        if self.penalty == "deprecated":
            if self.l1_ratios is None:
                warnings.warn(
                    (
                        "'l1_ratios=None' was deprecated in version 1.8 and will "
                        "trigger an error in 1.10. Use an array-like with values"
                        "in [0, 1] instead."
                    ),
                    FutureWarning,
                )
            if np.all(np.asarray(l1_ratios) == 0) or l1_ratios is None:
                penalty = "l2"
            elif np.all(np.asarray(l1_ratios) == 1):
                penalty = "l1"
            else:
                penalty = "elasticnet"
        else:
            penalty = self.penalty
            warnings.warn(
                (
                    "'penalty' was deprecated in version 1.8 and will be removed in"
                    " 1.10. To avoid this warning, leave 'penalty' set to its default"
                    " value and use 'l1_ratios' instead."
                    " Use l1_ratios=(0,) instead of penalty='l2' "
                    " and l1_ratios=(1,) instead of penalty='l1'."
                ),
                FutureWarning,
            )

        if self.use_legacy_attributes == "warn":
            warnings.warn(
                f"The fitted attributes of {self.__class__.__name__} will be "
                "simplified in scikit-learn 1.10 to remove redundancy. Set"
                "`use_legacy_attributes=False` to enable the new behavior now, or "
                "set it to `True` to silence this warning during the transition period "
                "while keeping the deprecated behavior for the time being. The default "
                "value of use_legacy_attributes will change from True to False in "
                f"scikit-learn 1.10. See the docstring of {self.__class__.__name__} "
                "for more details.",
                FutureWarning,
            )
            use_legacy_attributes = True
        else:
            use_legacy_attributes = self.use_legacy_attributes

        solver = _check_solver(self.solver, penalty, self.dual)

        if penalty == "elasticnet":
            if (
                l1_ratios is None
                or len(l1_ratios) == 0
                or any(
                    (
                        not isinstance(l1_ratio, numbers.Number)
                        or l1_ratio < 0
                        or l1_ratio > 1
                    )
                    for l1_ratio in l1_ratios
                )
            ):
                raise ValueError(
                    "l1_ratios must be an array-like of numbers between "
                    "0 and 1; got (l1_ratios=%r)" % l1_ratios
                )
            l1_ratios_ = l1_ratios
        else:
            if l1_ratios is not None and self.penalty != "deprecated":
                warnings.warn(
                    "l1_ratios parameter is only used when penalty "
                    "is 'elasticnet'. Got (penalty={})".format(penalty)
                )

            if l1_ratios is None:
                l1_ratios_ = [None]
            else:
                l1_ratios_ = l1_ratios

        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse="csr",
            dtype=np.float64,
            order="C",
            accept_large_sparse=solver not in ["liblinear", "sag", "saga"],
        )
        n_features = X.shape[1]
        check_classification_targets(y)

        class_weight = self.class_weight

        # Encode for string labels
        label_encoder = LabelEncoder().fit(y)

        # The original class labels
        classes_only_pos_if_binary = self.classes_ = label_encoder.classes_
        n_classes = len(self.classes_)
        is_binary = n_classes == 2

        if n_classes < 2:
            raise ValueError(
                "This solver needs samples of at least 2 classes"
                " in the data, but the data contains only one"
                f" class: {self.classes_[0]}."
            )

        if solver in ["sag", "saga"]:
            max_squared_sum = row_norms(X, squared=True).max()
        else:
            max_squared_sum = None

        if _routing_enabled():
            routed_params = process_routing(
                self,
                "fit",
                sample_weight=sample_weight,
                **params,
            )
        else:
            routed_params = Bunch()
            routed_params.splitter = Bunch(split={})
            routed_params.scorer = Bunch(score=params)
            if sample_weight is not None:
                routed_params.scorer.score["sample_weight"] = sample_weight

        # init cross-validation generator
        cv = check_cv(self.cv, y, classifier=True)
        folds = list(cv.split(X, y, **routed_params.splitter.split))

        if isinstance(class_weight, dict):
            if not (set(class_weight.keys()) <= set(self.classes_)):
                msg = (
                    "The given class_weight dict must have the class labels as keys; "
                    f"classes={self.classes_} but key={class_weight.keys()}"
                )
                raise ValueError(msg)
        elif class_weight == "balanced":
            # compute the class weights for the entire dataset y
            class_weight = compute_class_weight(
                class_weight,
                classes=self.classes_,
                y=y,
                sample_weight=sample_weight,
            )
            class_weight = dict(zip(self.classes_, class_weight))

        if is_binary:
            n_classes = 1
            classes_only_pos_if_binary = classes_only_pos_if_binary[1:]

        path_func = delayed(_log_reg_scoring_path)

        # The SAG solver releases the GIL so it's more efficient to use
        # threads for this solver.
        if self.solver in ["sag", "saga"]:
            prefer = "threads"
        else:
            prefer = "processes"

        fold_coefs_ = Parallel(n_jobs=self.n_jobs, verbose=self.verbose, prefer=prefer)(
            path_func(
                X,
                y,
                train,
                test,
                classes=self.classes_,
                Cs=self.Cs,
                fit_intercept=self.fit_intercept,
                penalty=penalty,
                dual=self.dual,
                solver=solver,
                tol=self.tol,
                max_iter=self.max_iter,
                verbose=self.verbose,
                class_weight=class_weight,
                scoring=self.scoring,
                intercept_scaling=self.intercept_scaling,
                random_state=self.random_state,
                max_squared_sum=max_squared_sum,
                sample_weight=sample_weight,
                l1_ratio=l1_ratio,
                score_params=routed_params.scorer.score,
            )
            for train, test in folds
            for l1_ratio in l1_ratios_
        )

        # fold_coefs_ is a list and would have shape (n_folds * n_l1_ratios, ..)
        # After reshaping,
        # - coefs_paths is of shape (n_classes, n_folds, n_Cs, n_l1_ratios, n_features)
        # - scores is of shape (n_classes, n_folds, n_Cs, n_l1_ratios)
        # - n_iter is of shape (1, n_folds, n_Cs, n_l1_ratios)
        coefs_paths, Cs, scores, n_iter_ = zip(*fold_coefs_)
        self.Cs_ = Cs[0]  # the same for all folds and l1_ratios
        if is_binary:
            coefs_paths = np.reshape(
                coefs_paths, (len(folds), len(l1_ratios_), len(self.Cs_), -1)
            )
            # coefs_paths.shape = (n_folds, n_l1_ratios, n_Cs, n_features)
            coefs_paths = np.swapaxes(coefs_paths, 1, 2)[None, ...]
        else:
            coefs_paths = np.reshape(
                coefs_paths, (len(folds), len(l1_ratios_), len(self.Cs_), n_classes, -1)
            )
            # coefs_paths.shape = (n_folds, n_l1_ratios, n_Cs, n_classes, n_features)
            coefs_paths = np.moveaxis(coefs_paths, (0, 1, 3), (1, 3, 0))
        # n_iter_.shape = (n_folds, n_l1_ratios, n_Cs)
        n_iter_ = np.reshape(n_iter_, (len(folds), len(l1_ratios_), len(self.Cs_)))
        self.n_iter_ = np.swapaxes(n_iter_, 1, 2)[None, ...]
        # scores.shape = (n_folds, n_l1_ratios, n_Cs)
        scores = np.reshape(scores, (len(folds), len(l1_ratios_), len(self.Cs_)))
        scores = np.swapaxes(scores, 1, 2)[None, ...]
        # repeat same scores across all classes
        scores = np.tile(scores, (n_classes, 1, 1, 1))
        self.scores_ = dict(zip(classes_only_pos_if_binary, scores))
        self.coefs_paths_ = dict(zip(classes_only_pos_if_binary, coefs_paths))

        self.C_ = list()
        self.l1_ratio_ = list()
        self.coef_ = np.empty((n_classes, n_features))
        self.intercept_ = np.zeros(n_classes)

        # All scores are the same across classes
        scores = self.scores_[classes_only_pos_if_binary[0]]

        if self.refit:
            # best_index over folds
            scores_sum = scores.sum(axis=0)  # shape (n_cs, n_l1_ratios)
            best_index = np.unravel_index(np.argmax(scores_sum), scores_sum.shape)

            C_ = self.Cs_[best_index[0]]
            self.C_.append(C_)

            l1_ratio_ = l1_ratios_[best_index[1]]
            self.l1_ratio_.append(l1_ratio_)

            if is_binary:
                coef_init = np.mean(coefs_paths[0, :, *best_index, :], axis=0)
            else:
                coef_init = np.mean(coefs_paths[:, :, *best_index, :], axis=1)

            # Note that y is label encoded
            w, _, _ = _logistic_regression_path(
                X,
                y,
                classes=self.classes_,
                Cs=[C_],
                solver=solver,
                fit_intercept=self.fit_intercept,
                coef=coef_init,
                max_iter=self.max_iter,
                tol=self.tol,
                penalty=penalty,
                class_weight=class_weight,
                verbose=max(0, self.verbose - 1),
                random_state=self.random_state,
                check_input=False,
                max_squared_sum=max_squared_sum,
                sample_weight=sample_weight,
                l1_ratio=l1_ratio_,
            )
            w = w[0]

        else:
            # Take the best scores across every fold and the average of
            # all coefficients corresponding to the best scores.
            n_folds, n_cs, n_l1_ratios = scores.shape
            scores = scores.reshape(n_folds, -1)  # (n_folds, n_cs * n_l1_ratios)
            best_indices = np.argmax(scores, axis=1)  # (n_folds,)
            best_indices = np.unravel_index(best_indices, (n_cs, n_l1_ratios))
            best_indices = list(zip(*best_indices))  # (n_folds, 2)
            # each row of best_indices has the 2 indices for Cs and l1_ratios
            if is_binary:
                w = np.mean(
                    [coefs_paths[0, i, *best_indices[i], :] for i in range(len(folds))],
                    axis=0,
                )
            else:
                w = np.mean(
                    [
                        coefs_paths[:, i, best_indices[i][0], best_indices[i][1], :]
                        for i in range(len(folds))
                    ],
                    axis=0,
                )

            best_indices = np.asarray(best_indices)
            best_indices_C = best_indices[:, 0]
            self.C_.append(np.mean(self.Cs_[best_indices_C]))

            if penalty == "elasticnet":
                best_indices_l1 = best_indices[:, 1]
                self.l1_ratio_.append(np.mean(l1_ratios_[best_indices_l1]))
            else:
                self.l1_ratio_.append(None)

        if is_binary:
            self.coef_ = w[:, :n_features] if w.ndim == 2 else w[:n_features][None, :]
            if self.fit_intercept:
                self.intercept_[0] = w[0, -1] if w.ndim == 2 else w[-1]
        else:
            self.C_ = np.tile(self.C_, n_classes)
            self.l1_ratio_ = np.tile(self.l1_ratio_, n_classes)
            self.coef_ = w[:, :n_features]
            if self.fit_intercept:
                self.intercept_ = w[:, -1]

        self.C_ = np.asarray(self.C_)
        self.l1_ratio_ = np.asarray(self.l1_ratio_)
        self.l1_ratios_ = np.asarray(l1_ratios_)
        if l1_ratios is None:
            # if elasticnet was not used, remove the l1_ratios dimension of some
            # attributes
            for cls, coefs_path in self.coefs_paths_.items():
                self.coefs_paths_[cls] = coefs_path[:, :, 0, :]
            for cls, score in self.scores_.items():
                self.scores_[cls] = score[:, :, 0]
            self.n_iter_ = self.n_iter_[:, :, :, 0]

        if not use_legacy_attributes:
            n_folds = len(folds)
            n_cs = self.Cs_.size
            n_dof = X.shape[1] + int(self.fit_intercept)
            self.C_ = float(self.C_[0])
            newpaths = np.concatenate(list(self.coefs_paths_.values()))
            newscores = self.scores_[
                classes_only_pos_if_binary[0]
            ]  # same for all classes
            newniter = self.n_iter_[0]
            if l1_ratios is None:
                if n_classes <= 2:
                    newpaths = newpaths.reshape(1, n_folds, n_cs, 1, n_dof)
                else:
                    newpaths = newpaths.reshape(n_classes, n_folds, n_cs, 1, n_dof)
                newscores = newscores.reshape(n_folds, n_cs, 1)
                newniter = newniter.reshape(n_folds, n_cs, 1)
                if self.penalty == "l1":
                    self.l1_ratio_ = 1.0
                else:
                    self.l1_ratio_ = 0.0
            else:
                n_l1_ratios = len(self.l1_ratios_)
                self.l1_ratio_ = float(self.l1_ratio_[0])
                if n_classes <= 2:
                    newpaths = newpaths.reshape(1, n_folds, n_cs, n_l1_ratios, n_dof)
                else:
                    newpaths = newpaths.reshape(
                        n_classes, n_folds, n_cs, n_l1_ratios, n_dof
                    )
            # newpaths.shape = (n_classes, n_folds, n_cs, n_l1_ratios, n_dof)
            # self.coefs_paths_.shape should be
            # (n_folds, n_l1_ratios, n_cs, n_classes, n_dof)
            self.coefs_paths_ = np.moveaxis(newpaths, (0, 1, 3), (3, 0, 1))
            # newscores.shape = (n_folds, n_cs, n_l1_ratios)
            # self.scores_.shape should be (n_folds, n_l1_ratios, n_cs)
            self.scores_ = np.moveaxis(newscores, (1, 2), (2, 1))
            self.n_iter_ = np.moveaxis(newniter, (1, 2), (2, 1))

        return self

    def score(self, X, y, sample_weight=None, **score_params):
        """Score using the `scoring` option on the given test data and labels.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test samples.

        y : array-like of shape (n_samples,)
            True labels for X.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        **score_params : dict
            Parameters to pass to the `score` method of the underlying scorer.

            .. versionadded:: 1.4

        Returns
        -------
        score : float
            Score of self.predict(X) w.r.t. y.
        """
        _raise_for_params(score_params, self, "score")

        scoring = self._get_scorer()
        if _routing_enabled():
            routed_params = process_routing(
                self,
                "score",
                sample_weight=sample_weight,
                **score_params,
            )
        else:
            routed_params = Bunch()
            routed_params.scorer = Bunch(score={})
            if sample_weight is not None:
                routed_params.scorer.score["sample_weight"] = sample_weight

        return scoring(
            self,
            X,
            y,
            **routed_params.scorer.score,
        )

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.4

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """

        router = (
            MetadataRouter(owner=self)
            .add_self_request(self)
            .add(
                splitter=self.cv,
                method_mapping=MethodMapping().add(caller="fit", callee="split"),
            )
            .add(
                scorer=self._get_scorer(),
                method_mapping=MethodMapping()
                .add(caller="score", callee="score")
                .add(caller="fit", callee="score"),
            )
        )
        return router

    def _get_scorer(self):
        """Get the scorer based on the scoring method specified.
        The default scoring method is `accuracy`.
        """
        scoring = self.scoring or "accuracy"
        return get_scorer(scoring)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = True
        return tags
