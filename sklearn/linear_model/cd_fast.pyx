# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
#
# License: BSD Style.

cimport numpy as np
import numpy as np
import numpy.linalg as linalg
cimport cython
from cpython cimport bool
import warnings

cdef extern from "math.h":
    double fabs(double f)
    double sqrt(double f)

cdef inline double fmax(double x, double y):
    if x > y:
        return x
    return y

cdef inline double fsign(double f):
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0

cdef extern from "cblas.h":
    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY)
    #enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
    #enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
    #                      AtlasConj=114};
    void dgemv "cblas_dgemv"(int order, int transpose, int M, int N,
                             double alpha, double *A, int lda, double *X,
                             int incX, double beta, double *Y, int incY)
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)


ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_std(unsigned int n_samples,
               unsigned int n_features,
               np.ndarray[DOUBLE, ndim=1] X_data,
               np.ndarray[INTEGER, ndim=1] X_indices,
               np.ndarray[INTEGER, ndim=1] X_indptr,
               np.ndarray[DOUBLE, ndim=1] X_mean=None):
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int nnz_ii
    cdef double X_sum_ii
    cdef double X_mean_ii
    cdef double diff
    cdef double X_std_ii

    cdef np.ndarray[DOUBLE, ndim = 1] X_std = np.zeros(n_features, np.float64)

    if X_mean is None:
        X_mean = np.zeros(n_features, np.float64)

        for ii in xrange(n_features):
            # Computes the mean
            X_sum_ii = 0.0
            for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                X_sum_ii += X_data[jj]
            X_mean[ii] = X_sum_ii / n_samples

    for ii in xrange(n_features):
        X_mean_ii = X_mean[ii]
        X_sum_ii = 0.0
        nnz_ii = 0
        for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
            diff = X_data[jj] - X_mean_ii
            X_sum_ii += diff * diff
            nnz_ii += 1

        X_std[ii] = (X_sum_ii + (n_samples - nnz_ii) * X_mean_ii * X_mean_ii)
    return np.sqrt(X_std)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double calculate_gap(np.ndarray[DOUBLE, ndim=1] w,
                            double l1_reg, double l2_reg,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y, bint positive):
    cdef double gap

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R

    # not efficient
    R = y - np.dot(X, w)

    XtA = np.dot(X.T, R) - l2_reg * w
    if positive:
        dual_norm_XtA = np.max(XtA)
    else:
        dual_norm_XtA = linalg.norm(XtA, np.inf)

    # TODO: use squared L2 norm directly
    R_norm = linalg.norm(R)
    w_norm = linalg.norm(w, 2)
    if (dual_norm_XtA > l1_reg):
        const = l1_reg / dual_norm_XtA
        A_norm = R_norm * const
        gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
    else:
        const = 1.0
        gap = R_norm ** 2

    gap += l1_reg * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * l2_reg * (1 + const ** 2) * (w_norm ** 2)
    return gap


def create_mapping(length, nz_index):
    cdef np.ndarray[INTEGER, ndim=1] map_to_ac = np.arange(length, dtype=np.int32)

    for (counter, nz) in enumerate(nz_index):
        tmp = map_to_ac[counter]
        map_to_ac[counter] = map_to_ac[nz]
        map_to_ac[nz] = tmp

    map_back = np.argsort(map_to_ac)
    return map_to_ac, map_back


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double l1_reg, double l2_reg,
                            np.ndarray[DOUBLE, ndim=2, mode='fortran'] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol, bint positive=False,
                            int memory_limit=250000):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)
    cdef double tmp
    cdef double tmp_gradient
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int m_pos
    cdef unsigned int org_pos
    cdef unsigned int n_iter
    cdef bint use_cache = False
    cdef bint initialize_cache = False
    cdef bint search_missing_feature = False
    cdef bint over_all = True
    cdef is_cached
    cdef int n_active_features = 0
    cdef int n_cached_features

    if l1_reg == 0:
        warnings.warn("Coordinate descent with l2_reg=0 may lead to unexpected"
            " results and is discouraged.")

    if memory_limit < n_features * n_features:
        warnings.warn("Allowed memory is not sufficient "
            "for caching some values need to be recalculated.")

    tol = tol * linalg.norm(y) ** 2

    cdef np.ndarray[DOUBLE, ndim=1] Xy = np.dot(X.T, y)
    cdef np.ndarray[DOUBLE, ndim=1] gradient = np.zeros(n_features, dtype=np.float64)
    cdef np.ndarray[INTEGER, ndim=1] nz_index = np.arange(n_features, dtype=np.int32)
    cdef np.ndarray[INTEGER, ndim=1] map_back = np.arange(n_features, dtype=np.int32)
    cdef np.ndarray[INTEGER, ndim=1] map_to_ac = np.arange(n_features, dtype=np.int32)
    cdef np.ndarray[INTEGER, ndim=1] active_set = np.arange(n_features, dtype=np.int32)
    cdef np.ndarray[INTEGER, ndim=1] iter_range = np.arange(n_features, dtype=np.int32)
    cdef np.ndarray[DOUBLE, ndim=2, mode='c'] feature_inner_product
    cdef np.ndarray[DOUBLE, ndim=1] tmp_feature_inner_product = np.zeros(n_features, dtype=np.float64)

    cdef int row_major = 101
    cdef int col_major = 102
    cdef int no_trans = 111
    cdef int trans = 112

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0

        n_active_features = len(active_set)

        # black magic conditions
        if n_iter > 2 and n_active_features > 2:
            active_set = np.nonzero(w)[0]
            iter_range = active_set
            n_active_features = len(active_set)

        # check if memory is now sufficient for caching
        if not use_cache:
            if n_active_features ** 2 <= memory_limit:
                nz_index = active_set
                iter_range = np.arange(n_active_features, dtype=np.int32)
                feature_inner_product = \
                    np.zeros(shape=(n_features, n_active_features),
                             dtype=np.float64, order='C')
                # resize
                map_to_ac, map_back = create_mapping(n_features, nz_index)
                w = w[map_to_ac]
                gradient = gradient[map_to_ac]
                n_cached_features = n_active_features
                use_cache = True
                initialize_cache = True
                over_all = False

        if search_missing_feature:
            iter_range = np.arange(n_features, dtype=np.int32)

        for ii in iter_range:  # Loop over coordinates

            if norm_cols_X[map_back[ii]] == 0.0:
                continue

            if use_cache:
                m_pos = map_to_ac[ii]
            else:
                m_pos = ii

            if over_all:
                org_pos = ii
            else:
                org_pos = map_to_ac[ii]

            w_ii = w[ii]  # Store previous value

            # if feature is not located at the beginning of the array it's
            # not cached
            is_cached = m_pos < n_cached_features

            # initial calculation
            if initialize_cache:
                #tmp_feature_inner_product = np.dot(X[:, nz_index[ii]], X)
                dgemv(col_major, trans, n_samples, n_features,
                          1, &X[0,0], n_samples, &X[0, org_pos], 
                          1, 0, &tmp_feature_inner_product[0], 1)
                feature_inner_product[:, ii] = tmp_feature_inner_product[map_to_ac]
                #gradient[ii] = Xy[nz_index[ii]] - \
                #        np.dot(feature_inner_product[:, ii], w)
                gradient[ii] = Xy[org_pos] - \
                        ddot(n_active_features, &feature_inner_product[0,ii], n_active_features, &w[0], 1)

            if not use_cache:
                #tmp_feature_inner_product = np.dot(X[:, ii], X)
                #gradient[ii] = Xy[ii] - \
                #         np.dot(tmp_feature_inner_product, w)
                dgemv(col_major, trans, n_samples, n_features,
                          1, &X[0,0], n_samples, &X[0,org_pos], 
                          1, 0, &tmp_feature_inner_product[0], 1)

                gradient[m_pos] = Xy[org_pos] - \
                        ddot(n_features, &tmp_feature_inner_product[0],1 , &w[0], 1)

            if search_missing_feature:
                if not is_cached:
                    dgemv(col_major, trans, n_samples, n_features,
                              1, &X[0,0], n_samples, &X[0,org_pos], 
                              1, 0, &tmp_feature_inner_product[0], 1)

                    gradient[m_pos] = Xy[org_pos] - \
                            ddot(n_features, &tmp_feature_inner_product[0],1 , &w[0], 1)

            tmp = gradient[m_pos] + w_ii * norm_cols_X[org_pos]

            if positive and tmp < 0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - l1_reg, 0) \
                    / (norm_cols_X[org_pos] + l2_reg)

            # update gradients, if w changed
            if w_ii != w[ii]:
                if search_missing_feature:
                    # add feature to the active-set that improved the objective 
                    #  but has not bin in the set before
                    pass
                if use_cache:
                    # gradient -= feature_inner_product[ii, :] * \
                    #                                    (w[ii] - w_ii)
                    daxpy(n_cached_features, -(w[ii] - w_ii),
                          &feature_inner_product[ii, 0], 1, &gradient[0], 1)

            # update the maximum absolute coefficient update
            d_w_ii = fabs(w[ii] - w_ii)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            if fabs(w[ii]) > w_max:
                w_max = fabs(w[ii])

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion
            if use_cache:
                gap = calculate_gap(w[map_back], l1_reg, l2_reg, X, y, positive)
            else:
                gap = calculate_gap(w, l1_reg, l2_reg, X, y, positive)

            if gap < tol:
                # return if we reached desired tolerance
                break
            else:
                print "dual gap check failed"
                search_missing_feature = True
                over_all = True

        initialize_cache = False

    if use_cache:
        w = w[map_back]
    return w, gap, tol


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=1] X_data,
                            np.ndarray[INTEGER, ndim=1] X_indices,
                            np.ndarray[INTEGER, ndim=1] X_indptr,
                            np.ndarray[DOUBLE, ndim=1] y,
                            np.ndarray[DOUBLE, ndim=1] X_mean,
                            int max_iter, double tol, bint positive=False):
    """Cython version of the coordinate descent algorithm for Elastic-Net

    We minimize:

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # compute norms of the columns of X
    cdef unsigned int ii
    cdef np.ndarray[DOUBLE, ndim = 1] norm_cols_X = np.zeros(n_features,
                                                           np.float64)
    for ii in xrange(n_features):
        norm_cols_X[ii] = ((X_data[X_indptr[ii]:X_indptr[ii + 1]] - \
        X_mean[ii]) ** 2).sum() + \
        (n_samples - X_indptr[ii + 1] + X_indptr[ii]) * X_mean[ii] ** 2

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim = 1] R

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double X_mean_ii
    cdef double R_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef bint center = (X_mean != 0).any()

    # initialize the residuals
    R = y.copy()

    for ii in xrange(n_features):
        # sparse X column / dense w dot product
        for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
            R[X_indices[jj]] -= X_data[jj] * w[ii]
        if center:
            R += X_mean[ii] * w[ii]

    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(max_iter):

        w_max = 0.0
        d_w_max = 0.0

        for ii in xrange(n_features):  # Loop over coordinates

            if norm_cols_X[ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value
            X_mean_ii = X_mean[ii]

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    R[X_indices[jj]] += X_data[jj] * w_ii
                if center:
                    for jj in xrange(n_samples):
                        R[jj] -= X_mean_ii * w[ii]

            # tmp = (X[:,ii] * R).sum()
            tmp = 0.0
            for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                tmp += R[X_indices[jj]] * X_data[jj]

            if center:
                R_sum = 0.0
                for jj in xrange(n_samples):
                    R_sum += R[jj]
                tmp -= R_sum * X_mean_ii

            if positive and tmp < 0.0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    R[X_indices[jj]] -= X_data[jj] * w[ii]

                if center:
                    for jj in xrange(n_samples):
                        R[jj] += X_mean_ii * w[ii]

            # update the maximum absolute coefficient update
            d_w_ii = fabs(w[ii] - w_ii)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            if w[ii] > w_max:
                w_max = w[ii]

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            # sparse X.T / dense R dot product
            X_T_R = np.zeros(n_features)
            for ii in xrange(n_features):
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    X_T_R[ii] += X_data[jj] * R[X_indices[jj]]
                X_T_R[ii] -= X_mean[ii] * R.sum()

            XtA = X_T_R - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            # TODO: use squared L2 norm directly
            R_norm = linalg.norm(R)
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
            else:
                const = 1.0
                gap = R_norm ** 2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_gram(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] Q,
                            np.ndarray[DOUBLE, ndim=1] q,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol, bool positive=False):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 w^T Q w - q^T w + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                      ----
        2                                        2

        which amount to the Elastic-Net problem when:
        Q = X^T X (Gram matrix)
        q = X^T y
    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = Q.shape[0]

    # initial value "Q w" which will be kept of up to date in the iterations
    cdef np.ndarray[DOUBLE, ndim=1] H = np.dot(Q, w)

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int n_iter

    cdef double y_norm2 = linalg.norm(y) ** 2
    tol = tol * y_norm2

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features):  # Loop over coordinates
            if Q[ii, ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value

            if w_ii != 0.0:
                # H -= w_ii * Q[ii]
                daxpy(n_features, -w_ii,
                      <DOUBLE*>(Q.data + ii * n_features * sizeof(DOUBLE)), 1,
                      <DOUBLE*>H.data, 1)

            tmp = q[ii] - H[ii]

            if positive and tmp < 0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (Q[ii, ii] + beta)

            if w[ii] != 0.0:
                # H +=  w[ii] * Q[ii] # Update H = X.T X w
                daxpy(n_features, w[ii],
                      <DOUBLE*>(Q.data + ii * n_features * sizeof(DOUBLE)), 1,
                      <DOUBLE*>H.data, 1)

            # update the maximum absolute coefficient update
            d_w_ii = fabs(w[ii] - w_ii)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            if fabs(w[ii]) > w_max:
                w_max = fabs(w[ii])

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            q_dot_w = np.dot(w, q)

            XtA = q - H - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            R_norm2 = y_norm2 + np.sum(w * H) - 2.0 * q_dot_w
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            gap += alpha * linalg.norm(w, 1) \
                   - const * y_norm2 \
                   + const * q_dot_w + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol

# ------------------ old code, to be removed later -------------------------


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_old(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol, bool positive=False):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int n_iter

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    R = y - np.dot(X, w)

    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features):  # Loop over coordinates
            if norm_cols_X[ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                daxpy(n_samples, w_ii,
                      <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                      <DOUBLE*>R.data, 1)

            # tmp = (X[:,ii]*R).sum()
            tmp = ddot(n_samples,
                       <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                       <DOUBLE*>R.data, 1)

            if positive and tmp < 0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                daxpy(n_samples, -w[ii],
                      <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                      <DOUBLE*>R.data, 1)

            # update the maximum absolute coefficient update
            d_w_ii = fabs(w[ii] - w_ii)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            if fabs(w[ii]) > w_max:
                w_max = fabs(w[ii])

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            XtA = np.dot(X.T, R) - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            # TODO: use squared L2 norm directly
            R_norm = linalg.norm(R)
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
            else:
                const = 1.0
                gap = R_norm ** 2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol