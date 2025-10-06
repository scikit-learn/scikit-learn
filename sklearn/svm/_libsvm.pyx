"""
Binding for libsvm_skl
----------------------

These are the bindings for libsvm_skl, which is a fork of libsvm[1]
that adds to libsvm some capabilities, like index of support vectors
and efficient representation of dense matrices.

These are low-level routines, but can be used for flexibility or
performance reasons. See sklearn.svm for a higher-level API.

Low-level memory management is done in libsvm_helper.c. If we happen
to run out of memory a MemoryError will be raised. In practice this is
not very helpful since high chances are malloc fails inside svm.cpp,
where no sort of memory checks are done.

[1] https://www.csie.ntu.edu.tw/~cjlin/libsvm/

Notes
-----
The signature mode='c' is somewhat superficial, since we already
check that arrays are C-contiguous in svm.py

Authors
-------
2010: Fabian Pedregosa <fabian.pedregosa@inria.fr>
      Gael Varoquaux <gael.varoquaux@normalesup.org>
"""

import  numpy as np
from libc.stdlib cimport free
from sklearn.utils._cython_blas cimport _dot
from sklearn.utils._typedefs cimport float64_t, int32_t, intp_t

include "_libsvm.pxi"

cdef extern from *:
    ctypedef struct svm_parameter:
        pass


################################################################################
# Internal variables
LIBSVM_KERNEL_TYPES = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']


################################################################################
# Wrapper functions

def fit(
    const float64_t[:, ::1] X,
    const float64_t[::1] Y,
    int svm_type=0,
    kernel='rbf',
    int degree=3,
    double gamma=0.1,
    double coef0=0.0,
    double tol=1e-3,
    double C=1.0,
    double nu=0.5,
    double epsilon=0.1,
    const float64_t[::1] class_weight=np.empty(0),
    const float64_t[::1] sample_weight=np.empty(0),
    int shrinking=1,
    int probability=0,
    double cache_size=100.,
    int max_iter=-1,
    int random_seed=0,
):
    """
    Train the model using libsvm (low-level method)

    Parameters
    ----------
    X : array-like, dtype=float64 of shape (n_samples, n_features)

    Y : array, dtype=float64 of shape (n_samples,)
        target vector

    svm_type : {0, 1, 2, 3, 4}, default=0
        Type of SVM: C_SVC, NuSVC, OneClassSVM, EpsilonSVR or NuSVR
        respectively.

    kernel : {'linear', 'rbf', 'poly', 'sigmoid', 'precomputed'}, default="rbf"
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int32, default=3
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial).

    gamma : float64, default=0.1
        Gamma parameter in rbf, poly and sigmoid kernels. Ignored by other
        kernels.

    coef0 : float64, default=0
        Independent parameter in poly/sigmoid kernel.

    tol : float64, default=1e-3
        Numeric stopping criterion (WRITEME).

    C : float64, default=1
        C parameter in C-Support Vector Classification.

    nu : float64, default=0.5
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].

    epsilon : double, default=0.1
        Epsilon parameter in the epsilon-insensitive loss function.

    class_weight : array, dtype=float64, shape (n_classes,), \
            default=np.empty(0)
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one.

    sample_weight : array, dtype=float64, shape (n_samples,), \
            default=np.empty(0)
        Weights assigned to each sample.

    shrinking : int, default=1
        Whether to use the shrinking heuristic.

    probability : int, default=0
        Whether to enable probability estimates.

    cache_size : float64, default=100
        Cache size for gram matrix columns (in megabytes).

    max_iter : int (-1 for no limit), default=-1
        Stop solver after this many iterations regardless of accuracy
        (XXX Currently there is no API to know whether this kicked in.)

    random_seed : int, default=0
        Seed for the random number generator used for probability estimates.

    Returns
    -------
    support : array of shape (n_support,)
        Index of support vectors.

    support_vectors : array of shape (n_support, n_features)
        Support vectors (equivalent to X[support]). Will return an
        empty array in the case of precomputed kernel.

    n_class_SV : array of shape (n_class,)
        Number of support vectors in each class.

    sv_coef : array of shape (n_class-1, n_support)
        Coefficients of support vectors in decision function.

    intercept : array of shape (n_class*(n_class-1)/2,)
        Intercept in decision function.

    probA, probB : array of shape (n_class*(n_class-1)/2,)
        Probability estimates, empty array for probability=False.

    n_iter : ndarray of shape (max(1, (n_class * (n_class - 1) // 2)),)
        Number of iterations run by the optimization routine to fit the model.
    """

    cdef svm_parameter param
    cdef svm_problem problem
    cdef svm_model *model
    cdef const char *error_msg
    cdef intp_t SV_len

    if len(sample_weight) == 0:
        sample_weight = np.ones(X.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == X.shape[0], (
            f"sample_weight and X have incompatible shapes: sample_weight has "
            f"{sample_weight.shape[0]} samples while X has {X.shape[0]}"
        )

    kernel_index = LIBSVM_KERNEL_TYPES.index(kernel)
    set_problem(
        &problem,
        <char*> &X[0, 0],
        <char*> &Y[0],
        <char*> &sample_weight[0],
        <intp_t*> X.shape,
        kernel_index,
    )
    if problem.x == NULL:
        raise MemoryError("Seems we've run out of memory")
    cdef int32_t[::1] class_weight_label = np.arange(
        class_weight.shape[0], dtype=np.int32
    )
    set_parameter(
        &param,
        svm_type,
        kernel_index,
        degree,
        gamma,
        coef0,
        nu,
        cache_size,
        C,
        tol,
        epsilon,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char*> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char*> &class_weight[0] if class_weight.size > 0 else NULL,
        max_iter,
        random_seed,
    )

    error_msg = svm_check_parameter(&problem, &param)
    if error_msg:
        # for SVR: epsilon is called p in libsvm
        error_repl = error_msg.decode('utf-8').replace("p < 0", "epsilon < 0")
        raise ValueError(error_repl)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    # this does the real work
    cdef int fit_status = 0
    with nogil:
        model = svm_train(&problem, &param, &fit_status, &blas_functions)

    # from here until the end, we just copy the data returned by
    # svm_train
    SV_len = get_l(model)
    n_class = get_nr(model)

    cdef int[::1] n_iter = np.empty(max(1, n_class * (n_class - 1) // 2), dtype=np.intc)
    copy_n_iter(<char*> &n_iter[0], model)

    cdef float64_t[:, ::1] sv_coef = np.empty((n_class-1, SV_len), dtype=np.float64)
    copy_sv_coef(<char*> &sv_coef[0, 0] if sv_coef.size > 0 else NULL, model)

    # the intercept is just model.rho but with sign changed
    cdef float64_t[::1] intercept = np.empty(
        int((n_class*(n_class-1))/2), dtype=np.float64
    )
    copy_intercept(<char*> &intercept[0], model, <intp_t*> intercept.shape)

    cdef int32_t[::1] support = np.empty(SV_len, dtype=np.int32)
    copy_support(<char*> &support[0] if support.size > 0 else NULL, model)

    # copy model.SV
    cdef float64_t[:, ::1] support_vectors
    if kernel_index == 4:
        # precomputed kernel
        support_vectors = np.empty((0, 0), dtype=np.float64)
    else:
        support_vectors = np.empty((SV_len, X.shape[1]), dtype=np.float64)
        copy_SV(
            <char*> &support_vectors[0, 0] if support_vectors.size > 0 else NULL,
            model,
            <intp_t*> support_vectors.shape,
        )

    cdef int32_t[::1] n_class_SV
    if svm_type == 0 or svm_type == 1:
        n_class_SV = np.empty(n_class, dtype=np.int32)
        copy_nSV(<char*> &n_class_SV[0] if n_class_SV.size > 0 else NULL, model)
    else:
        # OneClass and SVR are considered to have 2 classes
        n_class_SV = np.array([SV_len, SV_len], dtype=np.int32)

    cdef float64_t[::1] probA
    cdef float64_t[::1] probB
    if probability != 0:
        if svm_type < 2:  # SVC and NuSVC
            probA = np.empty(int(n_class*(n_class-1)/2), dtype=np.float64)
            probB = np.empty(int(n_class*(n_class-1)/2), dtype=np.float64)
            copy_probB(<char*> &probB[0], model, <intp_t*> probB.shape)
        else:
            probA = np.empty(1, dtype=np.float64)
            probB = np.empty(0, dtype=np.float64)
        copy_probA(<char*> &probA[0], model, <intp_t*> probA.shape)
    else:
        probA = np.empty(0, dtype=np.float64)
        probB = np.empty(0, dtype=np.float64)

    svm_free_and_destroy_model(&model)
    free(problem.x)

    return (
        support.base,
        support_vectors.base,
        n_class_SV.base,
        sv_coef.base,
        intercept.base,
        probA.base,
        probB.base,
        fit_status,
        n_iter.base,
    )


cdef void set_predict_params(
    svm_parameter *param,
    int svm_type,
    kernel,
    int degree,
    double gamma,
    double coef0,
    double cache_size,
    int probability,
    int nr_weight,
    char *weight_label,
    char *weight,
) except *:
    """Fill param with prediction time-only parameters."""

    # training-time only parameters
    cdef double C = 0.0
    cdef double epsilon = 0.1
    cdef int max_iter = 0
    cdef double nu = 0.5
    cdef int shrinking = 0
    cdef double tol = 0.1
    cdef int random_seed = -1

    kernel_index = LIBSVM_KERNEL_TYPES.index(kernel)

    set_parameter(
        param,
        svm_type,
        kernel_index,
        degree,
        gamma,
        coef0,
        nu,
        cache_size,
        C,
        tol,
        epsilon,
        shrinking,
        probability,
        nr_weight,
        weight_label,
        weight,
        max_iter,
        random_seed,
    )


def predict(
    const float64_t[:, ::1] X,
    const int32_t[::1] support,
    const float64_t[:, ::1] SV,
    const int32_t[::1] nSV,
    const float64_t[:, ::1] sv_coef,
    const float64_t[::1] intercept,
    const float64_t[::1] probA=np.empty(0),
    const float64_t[::1] probB=np.empty(0),
    int svm_type=0,
    kernel='rbf',
    int degree=3,
    double gamma=0.1,
    double coef0=0.0,
    const float64_t[::1] class_weight=np.empty(0),
    const float64_t[::1] sample_weight=np.empty(0),
    double cache_size=100.0,
):
    """
    Predict target values of X given a model (low-level method)

    Parameters
    ----------
    X : array-like, dtype=float of shape (n_samples, n_features)

    support : array of shape (n_support,)
        Index of support vectors in training set.

    SV : array of shape (n_support, n_features)
        Support vectors.

    nSV : array of shape (n_class,)
        Number of support vectors in each class.

    sv_coef : array of shape (n_class-1, n_support)
        Coefficients of support vectors in decision function.

    intercept : array of shape (n_class*(n_class-1)/2)
        Intercept in decision function.

    probA, probB : array of shape (n_class*(n_class-1)/2,)
        Probability estimates.

    svm_type : {0, 1, 2, 3, 4}, default=0
        Type of SVM: C_SVC, NuSVC, OneClassSVM, EpsilonSVR or NuSVR
        respectively.

    kernel : {'linear', 'rbf', 'poly', 'sigmoid', 'precomputed'}, default="rbf"
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int32, default=3
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial).

    gamma : float64, default=0.1
        Gamma parameter in rbf, poly and sigmoid kernels. Ignored by other
        kernels.

    coef0 : float64, default=0.0
        Independent parameter in poly/sigmoid kernel.

    Returns
    -------
    dec_values : array
        Predicted values.
    """
    cdef float64_t[::1] dec_values
    cdef svm_parameter param
    cdef svm_model *model
    cdef int rv

    cdef int32_t[::1] class_weight_label = np.arange(
        class_weight.shape[0], dtype=np.int32
    )

    set_predict_params(
        &param,
        svm_type,
        kernel,
        degree,
        gamma,
        coef0,
        cache_size,
        0,
        <int>class_weight.shape[0],
        <char*> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char*> &class_weight[0] if class_weight.size > 0 else NULL,
    )
    model = set_model(
        &param,
        <int> nSV.shape[0],
        <char*> &SV[0, 0] if SV.size > 0 else NULL,
        <intp_t*> SV.shape,
        <char*> &support[0] if support.size > 0 else NULL,
        <intp_t*> support.shape,
        <intp_t*> sv_coef.strides,
        <char*> &sv_coef[0, 0] if sv_coef.size > 0 else NULL,
        <char*> &intercept[0],
        <char*> &nSV[0],
        <char*> &probA[0] if probA.size > 0 else NULL,
        <char*> &probB[0] if probB.size > 0 else NULL,
    )
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    # TODO: use check_model
    try:
        dec_values = np.empty(X.shape[0])
        with nogil:
            rv = copy_predict(
                <char*> &X[0, 0],
                model,
                <intp_t*> X.shape,
                <char*> &dec_values[0],
                &blas_functions,
            )
        if rv < 0:
            raise MemoryError("We've run out of memory")
    finally:
        free_model(model)

    return dec_values.base


def predict_proba(
    const float64_t[:, ::1] X,
    const int32_t[::1] support,
    const float64_t[:, ::1] SV,
    const int32_t[::1] nSV,
    float64_t[:, ::1] sv_coef,
    float64_t[::1] intercept,
    float64_t[::1] probA=np.empty(0),
    float64_t[::1] probB=np.empty(0),
    int svm_type=0,
    kernel='rbf',
    int degree=3,
    double gamma=0.1,
    double coef0=0.0,
    float64_t[::1] class_weight=np.empty(0),
    float64_t[::1] sample_weight=np.empty(0),
    double cache_size=100.0,
):
    """
    Predict probabilities

    svm_model stores all parameters needed to predict a given value.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

    See sklearn.svm.predict for a complete list of parameters.

    Parameters
    ----------
    X : array-like, dtype=float of shape (n_samples, n_features)

    support : array of shape (n_support,)
        Index of support vectors in training set.

    SV : array of shape (n_support, n_features)
        Support vectors.

    nSV : array of shape (n_class,)
        Number of support vectors in each class.

    sv_coef : array of shape (n_class-1, n_support)
        Coefficients of support vectors in decision function.

    intercept : array of shape (n_class*(n_class-1)/2,)
        Intercept in decision function.

    probA, probB : array of shape (n_class*(n_class-1)/2,)
        Probability estimates.

    svm_type : {0, 1, 2, 3, 4}, default=0
        Type of SVM: C_SVC, NuSVC, OneClassSVM, EpsilonSVR or NuSVR
        respectively.

    kernel : {'linear', 'rbf', 'poly', 'sigmoid', 'precomputed'}, default="rbf"
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int32, default=3
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial).

    gamma : float64, default=0.1
        Gamma parameter in rbf, poly and sigmoid kernels. Ignored by other
        kernels.

    coef0 : float64, default=0.0
        Independent parameter in poly/sigmoid kernel.

    Returns
    -------
    dec_values : array
        Predicted values.
    """
    cdef float64_t[:, ::1] dec_values
    cdef svm_parameter param
    cdef svm_model *model
    cdef int32_t[::1] class_weight_label = np.arange(
        class_weight.shape[0], dtype=np.int32
    )
    cdef int rv

    set_predict_params(
        &param,
        svm_type,
        kernel,
        degree,
        gamma,
        coef0,
        cache_size,
        1,
        <int> class_weight.shape[0],
        <char*> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char*> &class_weight[0] if class_weight.size > 0 else NULL,
    )
    model = set_model(
        &param,
        <int> nSV.shape[0],
        <char*> &SV[0, 0] if SV.size > 0 else NULL,
        <intp_t*> SV.shape,
        <char*> &support[0],
        <intp_t*> support.shape,
        <intp_t*> sv_coef.strides,
        <char*> &sv_coef[0, 0],
        <char*> &intercept[0],
        <char*> &nSV[0],
        <char*> &probA[0] if probA.size > 0 else NULL,
        <char*> &probB[0] if probB.size > 0 else NULL,
    )

    cdef intp_t n_class = get_nr(model)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    try:
        dec_values = np.empty((X.shape[0], n_class), dtype=np.float64)
        with nogil:
            rv = copy_predict_proba(
                <char*> &X[0, 0],
                model,
                <intp_t*> X.shape,
                <char*> &dec_values[0, 0],
                &blas_functions,
            )
        if rv < 0:
            raise MemoryError("We've run out of memory")
    finally:
        free_model(model)

    return dec_values.base


def decision_function(
    const float64_t[:, ::1] X,
    const int32_t[::1] support,
    const float64_t[:, ::1] SV,
    const int32_t[::1] nSV,
    const float64_t[:, ::1] sv_coef,
    const float64_t[::1] intercept,
    const float64_t[::1] probA=np.empty(0),
    const float64_t[::1] probB=np.empty(0),
    int svm_type=0,
    kernel='rbf',
    int degree=3,
    double gamma=0.1,
    double coef0=0.0,
    const float64_t[::1] class_weight=np.empty(0),
    const float64_t[::1] sample_weight=np.empty(0),
    double cache_size=100.0,
):
    """
    Predict margin (libsvm name for this is predict_values)

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

    Parameters
    ----------
    X : array-like, dtype=float, size=[n_samples, n_features]

    support : array, shape=[n_support]
        Index of support vectors in training set.

    SV : array, shape=[n_support, n_features]
        Support vectors.

    nSV : array, shape=[n_class]
        Number of support vectors in each class.

    sv_coef : array, shape=[n_class-1, n_support]
        Coefficients of support vectors in decision function.

    intercept : array, shape=[n_class*(n_class-1)/2]
        Intercept in decision function.

    probA, probB : array, shape=[n_class*(n_class-1)/2]
        Probability estimates.

    svm_type : {0, 1, 2, 3, 4}, optional
        Type of SVM: C_SVC, NuSVC, OneClassSVM, EpsilonSVR or NuSVR
        respectively. 0 by default.

    kernel : {'linear', 'rbf', 'poly', 'sigmoid', 'precomputed'}, optional
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed. 'rbf' by default.

    degree : int32, optional
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial), 3 by default.

    gamma : float64, optional
        Gamma parameter in rbf, poly and sigmoid kernels. Ignored by other
        kernels. 0.1 by default.

    coef0 : float64, optional
        Independent parameter in poly/sigmoid kernel. 0 by default.

    Returns
    -------
    dec_values : array
        Predicted values.
    """
    cdef float64_t[:, ::1] dec_values
    cdef svm_parameter param
    cdef svm_model *model
    cdef intp_t n_class

    cdef int32_t[::1] class_weight_label = np.arange(
        class_weight.shape[0], dtype=np.int32
    )

    cdef int rv

    set_predict_params(
        &param,
        svm_type,
        kernel,
        degree,
        gamma,
        coef0,
        cache_size,
        0,
        <int> class_weight.shape[0],
        <char*> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char*> &class_weight[0] if class_weight.size > 0 else NULL,
    )

    model = set_model(
        &param,
        <int> nSV.shape[0],
        <char*> &SV[0, 0] if SV.size > 0 else NULL,
        <intp_t*> SV.shape,
        <char*> &support[0],
        <intp_t*> support.shape,
        <intp_t*> sv_coef.strides,
        <char*> &sv_coef[0, 0],
        <char*> &intercept[0],
        <char*> &nSV[0],
        <char*> &probA[0] if probA.size > 0 else NULL,
        <char*> &probB[0] if probB.size > 0 else NULL,
    )

    if svm_type > 1:
        n_class = 1
    else:
        n_class = get_nr(model)
        n_class = n_class * (n_class - 1) // 2
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    try:
        dec_values = np.empty((X.shape[0], n_class), dtype=np.float64)
        with nogil:
            rv = copy_predict_values(
                <char*> &X[0, 0],
                model,
                <intp_t*> X.shape,
                <char*> &dec_values[0, 0],
                n_class,
                &blas_functions,
            )
        if rv < 0:
            raise MemoryError("We've run out of memory")
    finally:
        free_model(model)

    return dec_values.base


def cross_validation(
    const float64_t[:, ::1] X,
    const float64_t[::1] Y,
    int n_fold,
    int svm_type=0,
    kernel='rbf',
    int degree=3,
    double gamma=0.1,
    double coef0=0.0,
    double tol=1e-3,
    double C=1.0,
    double nu=0.5,
    double epsilon=0.1,
    float64_t[::1] class_weight=np.empty(0),
    float64_t[::1] sample_weight=np.empty(0),
    int shrinking=0,
    int probability=0,
    double cache_size=100.0,
    int max_iter=-1,
    int random_seed=0,
):
    """
    Binding of the cross-validation routine (low-level routine)

    Parameters
    ----------

    X : array-like, dtype=float of shape (n_samples, n_features)

    Y : array, dtype=float of shape (n_samples,)
        target vector

    n_fold : int32
        Number of folds for cross validation.

    svm_type : {0, 1, 2, 3, 4}, default=0
        Type of SVM: C_SVC, NuSVC, OneClassSVM, EpsilonSVR or NuSVR
        respectively.

    kernel : {'linear', 'rbf', 'poly', 'sigmoid', 'precomputed'}, default='rbf'
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int32, default=3
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial).

    gamma : float64, default=0.1
        Gamma parameter in rbf, poly and sigmoid kernels. Ignored by other
        kernels.

    coef0 : float64, default=0.0
        Independent parameter in poly/sigmoid kernel.

    tol : float64, default=1e-3
        Numeric stopping criterion (WRITEME).

    C : float64, default=1
        C parameter in C-Support Vector Classification.

    nu : float64, default=0.5
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].

    epsilon : double, default=0.1
        Epsilon parameter in the epsilon-insensitive loss function.

    class_weight : array, dtype=float64, shape (n_classes,), \
            default=np.empty(0)
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one.

    sample_weight : array, dtype=float64, shape (n_samples,), \
            default=np.empty(0)
        Weights assigned to each sample.

    shrinking : int, default=1
        Whether to use the shrinking heuristic.

    probability : int, default=0
        Whether to enable probability estimates.

    cache_size : float64, default=100
        Cache size for gram matrix columns (in megabytes).

    max_iter : int (-1 for no limit), default=-1
        Stop solver after this many iterations regardless of accuracy
        (XXX Currently there is no API to know whether this kicked in.)

    random_seed : int, default=0
        Seed for the random number generator used for probability estimates.

    Returns
    -------
    target : array, float

    """

    cdef svm_parameter param
    cdef svm_problem problem
    cdef const char *error_msg

    if len(sample_weight) == 0:
        sample_weight = np.ones(X.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == X.shape[0], (
            f"sample_weight and X have incompatible shapes: sample_weight has "
            f"{sample_weight.shape[0]} samples while X has {X.shape[0]}"
        )

    if X.shape[0] < n_fold:
        raise ValueError("Number of samples is less than number of folds")

    # set problem
    kernel_index = LIBSVM_KERNEL_TYPES.index(kernel)
    set_problem(
        &problem,
        <char*> &X[0, 0],
        <char*> &Y[0],
        <char*> &sample_weight[0] if sample_weight.size > 0 else NULL,
        <intp_t*> X.shape,
        kernel_index,
    )
    if problem.x == NULL:
        raise MemoryError("Seems we've run out of memory")
    cdef int32_t[::1] class_weight_label = np.arange(
        class_weight.shape[0], dtype=np.int32
    )

    # set parameters
    set_parameter(
        &param,
        svm_type,
        kernel_index,
        degree,
        gamma,
        coef0,
        nu,
        cache_size,
        C,
        tol,
        tol,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char*> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char*> &class_weight[0] if class_weight.size > 0 else NULL,
        max_iter,
        random_seed,
    )

    error_msg = svm_check_parameter(&problem, &param)
    if error_msg:
        raise ValueError(error_msg)

    cdef float64_t[::1] target
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    try:
        target = np.empty((X.shape[0]), dtype=np.float64)
        with nogil:
            svm_cross_validation(
                &problem,
                &param,
                n_fold,
                <double *> &target[0],
                &blas_functions,
            )
    finally:
        free(problem.x)

    return target.base


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
