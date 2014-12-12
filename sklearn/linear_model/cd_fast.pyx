# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
#         Manoj Kumar <manojkumarsivaraj334@gmail.com>
#         Olivier Fercoq <olivier.fercoq@telecom-paristech.fr>
#
# Licence: BSD 3 clause

from libc.math cimport fabs, sqrt
cimport numpy as np
import numpy as np
import numpy.linalg as linalg

cimport cython
from cpython cimport bool
import warnings

ctypedef np.float64_t DOUBLE
ctypedef np.uint32_t UINT32_t

np.import_array()


# The following two functions are shamelessly copied from the tree code.

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline UINT32_t rand_int(UINT32_t end, UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


cdef inline double fmax(double x, double y) nogil:
    if x > y:
        return x
    return y

cdef inline double fmin(double x, double y) nogil:
    if x < y:
        return x
    return y

cdef inline double fsign(double f) nogil:
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0


cdef double abs_max(int n, double* a) nogil:
    """np.max(np.abs(a))"""
    cdef int i
    cdef double m = fabs(a[0])
    cdef double d
    for i in range(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m


cdef double max(int n, double* a) nogil:
    """np.max(a)"""
    cdef int i
    cdef double m = a[0]
    cdef double d
    for i in range(1, n):
        d = a[i]
        if d > m:
            m = d
    return m


cdef double diff_abs_max(int n, double* a, double* b) nogil:
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef double m = fabs(a[0] - b[0])
    cdef double d
    for i in range(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m


cdef extern from "cblas.h":
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY) nogil
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY
                             ) nogil
    double dasum "cblas_dasum"(int N, double *X, int incX) nogil
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                double *X, int incX, double *Y, int incY, double *A, int lda) nogil
    void dgemv "cblas_dgemv"(CBLAS_ORDER Order,
                      CBLAS_TRANSPOSE TransA, int M, int N,
                      double alpha, double *A, int lda,
                      double *X, int incX, double beta,
                      double *Y, int incY) nogil
    double dnrm2 "cblas_dnrm2"(int N, double *X, int incX) nogil
    void dcopy "cblas_dcopy"(int N, double *X, int incX, double *Y, int incY) nogil
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX) nogil


#####################################################################
# Dense matrices


# Function to compute the duality gap
cdef double enet_duality_gap( # Data
                        unsigned int n_samples,
                        unsigned int n_features,
                        unsigned int n_tasks,
                        DOUBLE * X_data,
                        DOUBLE * y_data,
                        DOUBLE * R_data,
                        DOUBLE * w_data,
                        # Variables intended to be modified
                        DOUBLE * XtA_data, 
                        double * dual_norm_XtA,
                        # Parameters
                        double alpha,
                        double beta,
                        bint positive,
                        ) nogil:
    
    cdef double R_norm2
    cdef double w_norm2
    cdef double ytA
    cdef double y_norm2
    cdef double l1_norm
    cdef double const
    cdef double gap
    cdef unsigned int ii
    
    # XtA = np.dot(X.T, R) - beta * w
    for ii in range(n_features):
        XtA_data[ii] = ddot(
            n_samples,
            X_data + ii * n_samples,
            1, R_data, 1) - beta * w_data[ii]
            
    if positive:
        dual_norm_XtA[0] = max(n_features, XtA_data)
    else:
        dual_norm_XtA[0] = abs_max(n_features, XtA_data)

    # R_norm2 = np.dot(R, R)
    R_norm2 = ddot(n_samples, R_data, 1, R_data, 1)

    # w_norm2 = np.dot(w, w)
    w_norm2 = ddot(n_features, w_data, 1, w_data, 1)

    # ytA = np.dot(y.T, R)
    ytA = ddot(n_samples, R_data, 1, y_data, n_tasks)

    # Determining the dual feasible point which is
    # proportional to R and closest to y / alpha
    if dual_norm_XtA[0] <= 0:
        if R_norm2 == 0:
            const = 1.
        else:
            const = ytA / R_norm2
    elif positive:
        const = alpha * fmin(ytA / (alpha * R_norm2),
                             1. / dual_norm_XtA[0])
    else:
        const = alpha * fmax(-1. / dual_norm_XtA[0],
                         fmin(ytA / (alpha * R_norm2),
                              1. / dual_norm_XtA[0]) )
       
    l1_norm = dasum(n_features, w_data, 1)

    gap = 0.5 * (1 + const ** 2) * (R_norm2 + beta * w_norm2) \
          + alpha * l1_norm - const * ytA
    
    return gap


# Function to compute the primal value
cdef double enet_primal_value( # Data
                        unsigned int n_samples,
                        unsigned int n_features,
                        DOUBLE * R_data,
                        DOUBLE * w_data,
                        # Parameters
                        double alpha,
                        double beta,
                        ) nogil:
    
    cdef double R_norm2
    cdef double w_norm2 = 0
    cdef double l1_norm
    cdef double fval
    
    # R_norm2 = np.dot(R, R)
    R_norm2 = ddot(n_samples, R_data, 1,
                   R_data, 1)

    # w_norm2 = np.dot(w, w)
    if beta>0:
        w_norm2 = ddot(n_features, w_data, 1,
                       w_data, 1)
    
    l1_norm = dasum(n_features, w_data, 1)

    fval = 0.5 * R_norm2 + alpha * l1_norm + 0.5 * beta * w_norm2 

    return fval

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_accelerated_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol,
                            object rng, bint random=1,
                            bint positive=0,
                            ):
    """Cython version of the accelerated coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)
    
    # additional variables
    cdef np.ndarray[DOUBLE, ndim=1] we = np.empty(n_features) # exploration variable
    cdef np.ndarray[DOUBLE, ndim=1] wc = np.zeros(n_features) # correction variable
    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R = np.empty(n_samples)  # normal variable
    cdef np.ndarray[DOUBLE, ndim=1] Re = np.empty(n_samples) # exploration variable
    cdef np.ndarray[DOUBLE, ndim=1] Rc = np.zeros(n_samples) # correction variable

    cdef np.ndarray[DOUBLE, ndim=1] XtA = np.empty(n_features)
    cdef np.ndarray[DOUBLE, ndim=1] Xty = np.empty(n_features)
    
    cdef double tmp
    cdef double we_ii
    cdef double wc_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA

    # initializing momentum parameter
    cdef double theta0 = 1. / fmax(1., n_features)
    cdef double theta = theta0
    
    cdef double test_restart = 0
    cdef unsigned int nnz
    cdef unsigned int nnze
    
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed
   
    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    if not random:
        warnings.warn("Accelerated coordinate descent with cyclic sampling"
            " may be slower than with random sampling and is discouraged.")     
    

    with nogil:

        # R = y - np.dot(X, w)
        for jj in range(n_samples):
            R[jj] = y[jj]
        for ii in range(n_features):
            # R -=  w[ii] * X[:,ii] # Update residual
            daxpy(n_samples, -w[ii],
                  <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                  1, <DOUBLE*>R.data, 1)
        # Re = R,  Rc = - np.dot(X, wc) = 0 (already initialized)
        for jj in range(n_samples):
            Re[jj] = R[jj]
        for ii in range(n_features):
            we[ii] = w[ii]

        ### scaling tolerance
        # tol *= np.dot(y, y)
        tol *= fmax(1e-15,ddot(n_samples, <DOUBLE*>y.data, n_tasks,
                    <DOUBLE*>y.data, n_tasks))

        ### main loop
        for n_iter in range(max_iter):
                  
            ###################
            # Accelerated Coordinate descent
            # cf. Fercoq, O., & Richtárik, P. (2013). Accelerated, parallel
            #  and proximal coordinate descent. arXiv preprint arXiv:1312.5799.

            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                we_ii = we[ii]  # Store previous value
                wc_ii = wc[ii]

                # partial derivative at wd = (theta ** 2 * wc + we) and recall at we
                # tmp = (X[:,ii]*Rd).sum() + norm_cols_X[ii] * theta/theta0 * we_ii
                tmp = (theta ** 2) * ddot(n_samples,
                            <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                            1, <DOUBLE*>Rc.data, 1) \
                         + ddot(n_samples,
                            <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                            1, <DOUBLE*>Re.data, 1) \
                         + norm_cols_X[ii] * theta / theta0 * we_ii

                if positive and tmp < 0:
                    we[ii] = 0.0
                else:
                    we[ii] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                             / (norm_cols_X[ii]*theta/theta0 + beta))

                if we[ii] != we_ii:
                    # Re -=  (we[ii]-we_ii) * X[:,ii] # Update exploration residual
                    daxpy(n_samples, -we[ii] + we_ii,
                          <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                          1, <DOUBLE*>Re.data, 1)
                    
                    wc[ii] -= (we[ii] - we_ii) * (1-theta/theta0)/(theta**2)
                    # Rc -=  (wc[ii]-wc_ii) * X[:,ii] # Update correction residual
                    daxpy(n_samples, -wc[ii] + wc_ii,
                          <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                          1, <DOUBLE*>Rc.data, 1)

                # update the test for restarting when minimizing a strongly convex functions
                eta = 1. / n_features 
                # The formula for eta is inspired from
                #  Glasmachers, T., & Dogan, Ü. (2014). Coordinate Descent with Online
                #  Adaptation of Coordinate Frequencies. arXiv preprint arXiv:1401.3737.
                test_restart = (1-eta) * test_restart \
                               - eta * (we[ii]-we_ii)*( we[ii] - (theta**2 * wc_ii+we_ii) )
                
                # update the maximum absolute coefficient update
                d_w_ii = fabs((theta ** 2) * (wc[ii] - wc_ii) + (we[ii] - we_ii))
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                # update theta
                theta = 0.5 * (sqrt((theta ** 4) + 4 * (theta ** 2)) - theta ** 2)


            #######################
            # Tricks to make accelerated coordinate descent efficient for the lasso
            
            tmp = (theta ** 2)/(1-theta) 
            for ii in range(n_features):
                w[ii] = tmp * wc[ii] + we[ii]
                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])
            for jj in range(n_samples):
                R[jj] = tmp * Rc[jj] + Re[jj]

            nnz=0
            nnze=0
            for ii in range(n_features):
                if w[ii] != 0:
                    nnz += 1
                if we[ii] != 0:
                    nnze +=1
            
            # if test_restart > 0: restart
            # (improves behaviour in presence of strong convexity)
            # cf. O’Donoghue, B., & Candes, E. (2012). Adaptive restart for accelerated
            #  gradient schemes. Foundations of Computational Mathematics, 1-18.
            if test_restart>0:
                for ii in range(n_features):
                    wc[ii] = 0
                    we[ii] = w[ii]
                for jj in range(n_samples):
                    Rc[jj] = 0
                    Re[jj] = R[jj]
                    
                theta = theta0
                                
            # if nnze < nnz and F(we) < F(w) then w <- we
            # (helps promoting sparsity to w)
            elif (nnze < nnz) \
                     and ( enet_primal_value(n_samples, n_features, <DOUBLE*>Re.data,
                                        <DOUBLE*>we.data, alpha, beta) <
                           enet_primal_value(n_samples, n_features, <DOUBLE*>R.data,
                                        <DOUBLE*>w.data, alpha, beta)):
                for ii in range(n_features):
                    wc[ii] = 0
                    w[ii] = we[ii]
                for jj in range(n_samples):
                    Rc[jj] = 0
                    R[jj] = Re[jj]
                d_w_max = 10000 * w_max * d_w_tol  # we've made a big change

            #######################
            # Termination criterion

            if (w_max == 0.0
                    or d_w_max / w_max < d_w_tol
                    or n_iter == max_iter - 1):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion

                gap = enet_duality_gap(n_samples, n_features, n_tasks,
                                  <DOUBLE*>X.data, <DOUBLE*>y.data,
                                  <DOUBLE*>R.data, <DOUBLE*>w.data,
                                  <DOUBLE*>XtA.data, &dual_norm_XtA,
                                  alpha, beta, positive)
                               
                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return w, gap, tol, n_iter + 1



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol,
                            object rng, bint random=0, bint positive=0):
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

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R = np.empty(n_samples)

    cdef np.ndarray[DOUBLE, ndim=1] XtA = np.empty(n_features)
    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA
    cdef unsigned int ii
    cdef unsigned int i
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    if True: #with nogil:

        # R = y - np.dot(X, w)
        for i in range(n_samples):
            R[i] = y[i] - ddot(n_features,
                               <DOUBLE*>(X.data + i * sizeof(DOUBLE)),
                               n_samples, <DOUBLE*>w.data, 1)

        # tol *= np.dot(y, y)
        tol *= fmax(1e-15, ddot(n_samples, <DOUBLE*>y.data, n_tasks,
                    <DOUBLE*>y.data, n_tasks))

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                # tmp = (X[:,ii]*R).sum() + norm_cols_X[ii] * w_ii
                tmp = ddot(n_samples,
                           <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                           1, <DOUBLE*>R.data, 1) + norm_cols_X[ii] * w_ii

                # soft thresholding
                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                             / (norm_cols_X[ii] + beta))

                if w[ii] != w_ii:
                    # R -=  (w[ii]-w_ii) * X[:,ii] # Update residual
                    daxpy(n_samples, -w[ii] + w_ii,
                          <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                          1, <DOUBLE*>R.data, 1)

                # update the maximum absolute coefficient update
                d_w_ii = fabs(w[ii] - w_ii)
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])

            if (w_max == 0.0
                    or d_w_max / w_max < d_w_tol
                    or n_iter == max_iter - 1):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion

                gap = enet_duality_gap(n_samples, n_features, n_tasks,
                                  <DOUBLE*>X.data, <DOUBLE*>y.data,
                                  <DOUBLE*>R.data, <DOUBLE*>w.data,
                                  <DOUBLE*>XtA.data, &dual_norm_XtA,
                                  alpha, beta, positive)
                
                if gap < tol:
                    # return if we reached desired tolerance
                    break
                              
    return w, gap, tol, n_iter + 1


#####################################################################
# Sparse matrices

# Function to compute the duality gap
cdef double sparse_enet_duality_gap( # Data
                        unsigned int n_samples,
                        unsigned int n_features,
                        unsigned int n_tasks,
                        DOUBLE * X_data,
                        int * X_indptr,
                        int * X_indices,
                        DOUBLE * X_mean_data,
                        DOUBLE * y_data,
                        DOUBLE * R_data,
                        DOUBLE * w_data,
                        # Variables intended to be modified
                        DOUBLE * XtA_data, 
                        double * dual_norm_XtA,
                        # Parameters
                        double alpha,
                        double beta,
                        bint positive,
                        ) nogil:
    
    cdef double R_norm2
    cdef double R_sum
    cdef double w_norm2 = 0
    cdef double ytA
    cdef double y_norm2
    cdef double l1_norm
    cdef double const
    cdef double gap
    cdef unsigned int ii

    # sparse X.T / dense R dot product
    for ii in range(n_features):
        XtA_data[ii] = 0.0
        for jj in range(X_indptr[ii], X_indptr[ii + 1]):
            XtA_data[ii] += X_data[jj] * R_data[X_indices[jj]]
        R_sum = 0.0
        for jj in range(n_samples):
            R_sum += R_data[jj]
        XtA_data[ii] -= X_mean_data[ii] * R_sum
        XtA_data[ii] -= beta * w_data[ii]

    if positive:
        dual_norm_XtA[0] = max(n_features, XtA_data)
    else:
        dual_norm_XtA[0] = abs_max(n_features, XtA_data)
        
    # R_norm2 = np.dot(R, R)
    R_norm2 = ddot(n_samples, R_data, 1, R_data, 1)
    
    # w_norm2 = np.dot(w, w)
    if beta>0:
        w_norm2 = ddot(n_features, w_data, 1, w_data, 1)
 
    # ytA = np.dot(y.T, R)
    ytA = ddot(n_samples, R_data, 1, y_data, n_tasks)

    # Determining the dual feasible point which is
    # proportional to R and closest to y / alpha
    if dual_norm_XtA[0] <= 0:
        if R_norm2 == 0:
            const = 1.
        else:
            const = ytA / R_norm2
    elif positive:
        const = alpha * fmin(ytA / (alpha * R_norm2),
                             1. / dual_norm_XtA[0])
    else:
        const = alpha * fmax(-1. / dual_norm_XtA[0],
                         fmin(ytA / (alpha * R_norm2),
                              1. / dual_norm_XtA[0]) )
       
    l1_norm = dasum(n_features, w_data, 1)

    gap = 0.5 * (1 + const ** 2) * (R_norm2 + beta * w_norm2) \
          + alpha * l1_norm - const * ytA

    
    return gap
        

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_accelerated_coordinate_descent(double[:] w,
                            double alpha, double beta,
                            double[:] X_data, int[:] X_indices,
                            int[:] X_indptr, double[:] y,
                            double[:] X_mean, int max_iter,
                            double tol, object rng, bint random=1,
                            bint positive=0):
    """Cython version of the accelerated coordinate descent algorithm for Elastic-Net

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
    cdef double[:] norm_cols_X = np.zeros(n_features, np.float64)

    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # additional variables
    cdef double[:] we = np.empty(n_features) # exploration variable
    cdef double[:] wc = np.zeros(n_features) # correction variable
    # initial value of the residuals
    cdef double[:] R = y.copy()             # normal variable
    cdef double[:] Re = np.empty(n_samples) # exploration variable
    cdef double[:] Rc = np.zeros(n_samples) # correction variable

    cdef double[:] XtA = np.zeros(n_features)

    cdef double tmp
    cdef double tmp2
    cdef double w_ii
    cdef double we_ii
    cdef double wc_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double X_mean_ii
    cdef double Re_sum
    cdef double Rc_sum
    cdef double normalize_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA
    cdef double theta0 = 1. / fmax(1., n_features)
    cdef double theta = theta0
    cdef double test_restart = 0.
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef unsigned int nnz
    cdef unsigned int nnze
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed
    cdef bint center = False

    if not random:
        warnings.warn("Accelerated coordinate descent with cyclic sampling"
            " may be slower than with random sampling and is discouraged.")

    with nogil:
        # center = (X_mean != 0).any()
        for ii in range(n_features):
            if X_mean[ii]:
                center = True
                break

        for ii in range(n_features):
            X_mean_ii = X_mean[ii]
            endptr = X_indptr[ii + 1]
            normalize_sum = 0.0
            w_ii = w[ii]

            for jj in range(startptr, endptr):
                normalize_sum += (X_data[jj] - X_mean_ii) ** 2
                R[X_indices[jj]] -= X_data[jj] * w_ii
            norm_cols_X[ii] = normalize_sum + \
                (n_samples - endptr + startptr) * X_mean_ii ** 2

            if center:
                for jj in range(n_samples):
                    R[jj] += X_mean_ii * w_ii
            startptr = endptr            

        # Re = R,  Rc = - np.dot(X, wc) = 0 (already initialized)
        for jj in range(n_samples):
            Re[jj] = R[jj]
        for ii in range(n_features):
            we[ii] = w[ii]

        # tol *= np.dot(y, y)
        tol *= fmax(1e-15, ddot(n_samples, <DOUBLE*>&y[0], 1, <DOUBLE*>&y[0], 1))

        # main loop
        for n_iter in range(max_iter):

            ###################
            # Accelerated Coordinate descent
            # cf. Fercoq, O., & Richtárik, P. (2013). Accelerated, parallel
            #  and proximal coordinate descent. arXiv preprint arXiv:1312.5799.

            w_max = 0.0
            d_w_max = 0.0

            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                startptr = X_indptr[ii]
                endptr = X_indptr[ii + 1]
                we_ii = we[ii]  # Store previous value
                wc_ii = wc[ii] 
                X_mean_ii = X_mean[ii]

                # compute partial derivative at wd = (theta ** 2 * wc + we) + recall at we
                # tmp = (X[:,ii]*Re).sum() 
                # tmp2 = (X[:,ii]*Rc).sum()
                # tmp = (theta**2) * tmp2 + tmp + norm_cols_X[ii]* theta/theta0 * we_ii
                tmp = 0.0
                tmp2 = 0.0
                for jj in range(startptr, endptr):
                    tmp += Re[X_indices[jj]] * X_data[jj]
                    tmp2 += Rc[X_indices[jj]] * X_data[jj]
                    
                if center:
                    Re_sum = 0.0
                    Rc_sum = 0.0
                    for jj in range(n_samples):
                        Re_sum += Re[jj]
                        Rc_sum += Rc[jj]
                    tmp -= Re_sum * X_mean_ii
                    tmp2 -= Rc_sum * X_mean_ii
                    
                tmp = (theta**2) * tmp2 + tmp \
                      + norm_cols_X[ii] * theta / theta0 * we_ii

                # soft thresholding
                if positive and tmp < 0.0:
                    we[ii] = 0.0
                else:
                    we[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                            / (norm_cols_X[ii] * theta / theta0 + beta)

                if we[ii] != we_ii:
                    # Re -=  (we[ii]-we_ii) * X[:,ii] # Update exploration residual
                    d_w_ii = we[ii] - we_ii
                    for jj in range(startptr, endptr):
                        Re[X_indices[jj]] -= X_data[jj] * d_w_ii

                    if center:
                        for jj in range(n_samples):
                            Re[jj] += X_mean_ii * d_w_ii

                    # Rc -=  (wc[ii]-wc_ii) * X[:,ii] # Update correction residual
                    d_w_ii = -(we[ii] - we_ii) * (1 - theta / theta0) / (theta**2)
                    wc[ii] += d_w_ii
                    for jj in range(startptr, endptr):
                        Rc[X_indices[jj]] -= X_data[jj] * d_w_ii

                    if center:
                        for jj in range(n_samples):
                            Rc[jj] += X_mean_ii * d_w_ii

                # update the test for restarting when minimizing a strongly convex functions
                eta = 1. / n_features
                # The formula for eta is inspired from
                #  Glasmachers, T., & Dogan, Ü. (2014). Coordinate Descent with Online
                #  Adaptation of Coordinate Frequencies. arXiv preprint arXiv:1401.3737.
                test_restart = (1-eta) * test_restart - eta * \
                               (we[ii]-we_ii)*( we[ii] - (theta**2 * wc_ii + we_ii) )


                # update the maximum absolute coefficient update
                d_w_ii = fabs((theta ** 2) * (wc[ii] - wc_ii) + (we[ii] - we_ii))
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                # update theta
                theta = 0.5 * (sqrt((theta ** 4) + 4. * (theta ** 2)) - theta ** 2)

            #######################
            # Tricks to make accelerated coordinate descent efficient for the lasso
            
            tmp = (theta ** 2) / (1 - theta) 
            for ii in range(n_features):
                w[ii] = tmp * wc[ii] + we[ii]
                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])
            for jj in range(n_samples):
                R[jj] = tmp * Rc[jj] + Re[jj]

            nnz=0
            nnze=0
            for ii in range(n_features):
                if w[ii] != 0:
                    nnz += 1
                if we[ii] != 0:
                    nnze +=1
            
            # if test_restart > 0: restart
            # (improves behaviour in presence of strong convexity)
            # cf. O’Donoghue, B., & Candes, E. (2012). Adaptive restart for accelerated
            #  gradient schemes. Foundations of Computational Mathematics, 1-18.
            if test_restart > 0:
                for ii in range(n_features):
                    wc[ii] = 0
                    we[ii] = w[ii]
                for jj in range(n_samples):
                    Rc[jj] = 0
                    Re[jj] = R[jj]
                theta = theta0

            # if nnze < nnz and F(we) < F(w) then w <- we
            # (helps promoting sparsity to w)
            elif (nnze < nnz) \
                     and ( enet_primal_value(n_samples, n_features, <DOUBLE*> &Re[0],
                                        <DOUBLE*> &we[0], alpha, beta) <
                           enet_primal_value(n_samples, n_features, <DOUBLE*> &R[0],
                                        <DOUBLE*> &w[0], alpha, beta) ):
                for i in range(n_features):
                    wc[i] = 0
                    w[i] = we[i]
                for j in range(n_samples):
                    Rc[j] = 0
                    R[j] = Re[j]
                d_w_max = 10000 * w_max * d_w_tol  # we've made a big change
            
            #######################
            # Termination criterion
                    
            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                gap = sparse_enet_duality_gap(n_samples, n_features, n_tasks,
                                         <DOUBLE*> &X_data[0], <int*> &X_indptr[0],
                                         <int*> &X_indices[0], <DOUBLE*> &X_mean[0],
                                         <DOUBLE*> &y[0], <DOUBLE*> &R[0],
                                         <DOUBLE*> &w[0], <DOUBLE*> &XtA[0],
                                         &dual_norm_XtA,
                                         alpha, beta, positive)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return w, gap, tol, n_iter + 1



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_coordinate_descent(double[:] w,
                            double alpha, double beta,
                            double[:] X_data, int[:] X_indices,
                            int[:] X_indptr, double[:] y,
                            double[:] X_mean, int max_iter,
                            double tol, object rng, bint random=0,
                            bint positive=0):
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
    cdef double[:] norm_cols_X = np.zeros(n_features, np.float64)

    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)
    
    # initial value of the residuals
    cdef double[:] R = y.copy()

    cdef double[:] XtA = np.zeros(n_features)

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double X_mean_ii
    cdef double R_sum
    cdef double normalize_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed
    cdef bint center = False

    with nogil:
        # center = (X_mean != 0).any()
        for ii in range(n_features):
            if X_mean[ii]:
                center = True
                break

        for ii in range(n_features):
            X_mean_ii = X_mean[ii]
            endptr = X_indptr[ii + 1]
            normalize_sum = 0.0
            w_ii = w[ii]

            for jj in range(startptr, endptr):
                normalize_sum += (X_data[jj] - X_mean_ii) ** 2
                R[X_indices[jj]] -= X_data[jj] * w_ii
            norm_cols_X[ii] = normalize_sum + \
                (n_samples - endptr + startptr) * X_mean_ii ** 2

            if center:
                for jj in range(n_samples):
                    R[jj] += X_mean_ii * w_ii
            startptr = endptr

        # tol *= np.dot(y, y)
        tol *= fmax(1e-15,ddot(n_samples, <DOUBLE*>&y[0], 1, <DOUBLE*>&y[0], 1))

        for n_iter in range(max_iter):

            w_max = 0.0
            d_w_max = 0.0

            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                startptr = X_indptr[ii]
                endptr = X_indptr[ii + 1]
                w_ii = w[ii]  # Store previous value
                X_mean_ii = X_mean[ii]

                # tmp = (X[:,ii]*R).sum() + norm_cols_X[ii] * w_ii
                tmp = 0.0
                for jj in range(startptr, endptr):
                    tmp += R[X_indices[jj]] * X_data[jj]
                    
                if center:
                    R_sum = 0.0
                    for jj in range(n_samples):
                        R_sum += R[jj]
                    tmp -= R_sum * X_mean_ii
                    
                tmp += norm_cols_X[ii] * w_ii

                # soft thresholding
                if positive and tmp < 0.0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                            / (norm_cols_X[ii] + beta)

                if w[ii] != w_ii:
                    # R -=  (w[ii]-w_ii) * X[:,ii] # Update residual
                    d_w_ii = w[ii] - w_ii
                    for jj in range(startptr, endptr):
                        R[X_indices[jj]] -= X_data[jj] * d_w_ii

                    if center:
                        for jj in range(n_samples):
                            R[jj] += X_mean_ii * d_w_ii

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

                gap = sparse_enet_duality_gap(n_samples, n_features, n_tasks,
                                         <DOUBLE*> &X_data[0], <int*> &X_indptr[0],
                                         <int*> &X_indices[0], <DOUBLE*> &X_mean[0],
                                         <DOUBLE*> &y[0], <DOUBLE*> &R[0],
                                         <DOUBLE*> &w[0], <DOUBLE*> &XtA[0],
                                         &dual_norm_XtA,
                                         alpha, beta, positive)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return w, gap, tol, n_iter + 1


###########################################################################
# Precomputed Gram matrices

# Function to compute the duality gap
cdef double gram_enet_duality_gap( # Data
                        unsigned int n_features,
                        unsigned int n_tasks,
                        DOUBLE * q_data,
                        double y_norm2,
                        DOUBLE * H_data,
                        DOUBLE * w_data,
                        # Variables intended to be modified
                        DOUBLE * XtA_data, 
                        double * dual_norm_XtA,
                        # Parameters
                        double alpha,
                        double beta,
                        bint positive,
                        ) nogil:

    cdef double q_dot_w
    cdef double R_norm2
    cdef double w_norm2 = 0
    cdef double l1_norm
    cdef double const
    cdef double gap
    cdef double tmp
 
    # q_dot_w = np.dot(w, q)
    # Note that increment in q is not 1 because the strides
    # vary if q is sliced from a 2-D array.
    q_dot_w = ddot(n_features, w_data, 1, q_data, n_tasks)
    
    for ii in range(n_features):
        XtA_data[ii] = q_data[ii] - H_data[ii] - beta * w_data[ii]
    if positive:
        dual_norm_XtA[0] = max(n_features, XtA_data)
    else:
        dual_norm_XtA[0] = abs_max(n_features, XtA_data)
        
    # tmp = np.sum(w * H)
    tmp = 0.0
    for ii in range(n_features):
        tmp += w_data[ii] * H_data[ii]
    R_norm2 = y_norm2 + tmp - 2.0 * q_dot_w

    # w_norm2 = np.dot(w, w)
    if beta>0:
        w_norm2 = ddot(n_features, w_data, 1, w_data, 1)

    # Determining the dual feasible point which is
    # proportional to R and closest to y / alpha
    if dual_norm_XtA[0] <= 0:
        if R_norm2 == 0:
            const = 1.
        else:
            const = (y_norm2 - q_dot_w) / R_norm2
    elif positive:
        const = alpha * fmin( (y_norm2 - q_dot_w) / (alpha * R_norm2),
                             1. / dual_norm_XtA[0])
    else:
        const = alpha * fmax(-1. / dual_norm_XtA[0],
                         fmin( (y_norm2 - q_dot_w) / (alpha * R_norm2),
                              1. / dual_norm_XtA[0]) )
        
    l1_norm = dasum(n_features, w_data, 1)

    gap = 0.5 * (1 + const ** 2) * (R_norm2 + beta * w_norm2) \
          + alpha * l1_norm - const * (y_norm2 - q_dot_w)

    return gap


# Function to compute the primal value
cdef double gram_enet_primal_value( # Data
                        unsigned int n_features,
                        unsigned int n_tasks,
                        DOUBLE * q_data,
                        double y_norm2,
                        DOUBLE * H_data,
                        DOUBLE * w_data,
                        # Parameters
                        double alpha,
                        double beta,
                        ) nogil:
    
    cdef double R_norm2
    cdef double w_norm2 = 0
    cdef double l1_norm
    cdef double fval
    cdef double tmp

    q_dot_w = ddot(n_features, w_data, 1, q_data, n_tasks)

    # tmp = np.sum(w * H)
    tmp = 0.0
    for ii in range(n_features):
        tmp += w_data[ii] * H_data[ii]
    R_norm2 = y_norm2 + tmp - 2.0 * q_dot_w

    # w_norm2 = np.dot(w, w)
    if beta>0:
        w_norm2 = ddot(n_features, w_data, 1,
                       w_data, 1)
    
    l1_norm = dasum(n_features, w_data, 1)

    fval = 0.5 * R_norm2 + alpha * l1_norm + 0.5 * beta * w_norm2 

    return fval


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_accelerated_coordinate_descent_gram(double[:] w, double alpha, double beta,
                                 double[:, :] Q, double[:] q, double[:] y,
                                 int max_iter, double tol, object rng,
                                 bint random=0, bint positive=0):
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
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # additional variables
    cdef double[:] we = w.copy()             # exploration variable
    cdef double[:] wc = np.zeros(n_features) # correction variable
    # initial value "Q w" which will be kept of up to date in the iterations
    cdef double[:] H = np.dot(Q, w)
    cdef double[:] He = H.copy()
    cdef double[:] Hc = np.zeros(n_features)

    cdef double[:] XtA = np.zeros(n_features)
    cdef double tmp
    cdef double w_ii
    cdef double we_ii
    cdef double wc_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA

    cdef double theta0 = 1. / fmax(1., n_features)
    cdef double theta = theta0
    cdef double test_restart = 0
    cdef unsigned int nnz
    cdef unsigned int nnze
    
    cdef unsigned int ii
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef double y_norm2 = np.dot(y, y)
    cdef double* Q_ptr = &Q[0, 0]
    cdef double* H_ptr = &H[0]
    cdef double* He_ptr = &He[0]
    cdef double* Hc_ptr = &Hc[0]
    tol = tol * fmax(1e-15, y_norm2)

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    with nogil:
        for n_iter in range(max_iter):
            
            # Accelerated Coordinate descent
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if Q[ii, ii] == 0.0:
                    continue

                we_ii = we[ii]  # Store previous value
                wc_ii = wc[ii]

                # partial derivative at wd = (theta ** 2 * wc + we) and recall at we
                tmp = q[ii] - (theta**2 * Hc[ii] + He[ii]) + Q[ii,ii]*theta/theta0 * we_ii

                if positive and tmp < 0:
                    we[ii] = 0.0
                else:
                    we[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (Q[ii, ii]*theta/theta0 + beta)

                if we[ii] != we_ii:
                    # He +=  (we[ii] - we_ii) * Q[ii] # Update He = X.T X we
                    daxpy(n_features, we[ii] - we_ii, Q_ptr + ii * n_features, 1,
                          He_ptr, 1)
                    
                    wc[ii] -= (we[ii] - we_ii) * (1-theta/theta0)/(theta**2)
                    # Hc +=  (wc[ii] - wc_ii) * Q[ii] # Update Hc = X.T X wc
                    daxpy(n_features, wc[ii] - wc_ii, Q_ptr + ii * n_features, 1,
                          Hc_ptr, 1)

                # update the test for restarting when minimizing a strongly convex functions
                eta = 1. / n_features 
                test_restart = (1-eta) * test_restart \
                               - eta * (we[ii]-we_ii)*( we[ii] - (theta**2 * wc_ii+we_ii) )
                
                # update the maximum absolute coefficient update
                d_w_ii = fabs((theta ** 2) * (wc[ii] - wc_ii) + (we[ii] - we_ii))
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                # update theta
                theta = 0.5 * (sqrt((theta ** 4) + 4. * (theta ** 2)) - theta ** 2)


            # Tricks to make accelerated coordinate descent efficient for the lasso
            tmp = (theta ** 2)/(1-theta) 
            for ii in range(n_features):
                w[ii] = tmp * wc[ii] + we[ii]
                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])
            for ii in range(n_features):
                H[ii] = tmp * Hc[ii] + He[ii]

            nnz=0
            nnze=0
            for ii in range(n_features):
                if w[ii] != 0:
                    nnz += 1
                if we[ii] != 0:
                    nnze +=1
            
            # if test_restart > 0: restart
            # (improves behaviour in presence of strong convexity)
            if test_restart>0:
                for ii in range(n_features):
                    wc[ii] = 0
                    we[ii] = w[ii]
                for ii in range(n_features):
                    Hc[ii] = 0
                    He[ii] = H[ii]
                    
                theta = theta0
                                
            # if nnze < nnz and F(we) < F(w) then w <- we
            # (helps promoting sparsity to w)
            elif (nnze < nnz) \
                     and ( gram_enet_primal_value(n_features, n_tasks,
                                             <DOUBLE*> &q[0], y_norm2,
                                             <DOUBLE*> &He[0], <DOUBLE*> &we[0],
                                             alpha, beta) <
                           gram_enet_primal_value(n_features, n_tasks,
                                             <DOUBLE*> &q[0], y_norm2,
                                             <DOUBLE*> &H[0], <DOUBLE*> &w[0],
                                             alpha, beta) ):
                for ii in range(n_features):
                    wc[ii] = 0
                    w[ii] = we[ii]
                for ii in range(n_features):
                    Hc[ii] = 0
                    H[ii] = He[ii]
                d_w_max = 10000 * w_max * d_w_tol  # we've made a big change

            # Termination criterion
            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                gap = gram_enet_duality_gap(n_features, n_tasks, <DOUBLE*> &q[0],
                                       y_norm2, <DOUBLE*> &H[0],
                                       <DOUBLE*> &w[0], <DOUBLE*> &XtA[0],
                                       &dual_norm_XtA,
                                       alpha, beta, positive)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return np.asarray(w), gap, tol, n_iter + 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_gram(double[:] w, double alpha, double beta,
                                 double[:, :] Q, double[:] q, double[:] y,
                                 int max_iter, double tol, object rng,
                                 bint random=0, bint positive=0):
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
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # initial value "Q w" which will be kept of up to date in the iterations
    cdef double[:] H = np.dot(Q, w)

    cdef double[:] XtA = np.zeros(n_features)
    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA
    cdef unsigned int ii
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef double y_norm2 = np.dot(y, y)
    cdef double* Q_ptr = &Q[0, 0]
    tol = tol * fmax(1e-15, y_norm2)

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    with nogil:
        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if Q[ii, ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                tmp = q[ii] - H[ii] + Q[ii, ii] * w_ii

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (Q[ii, ii] + beta)

                if w[ii] != w_ii:
                    # H +=  (w[ii] - w_ii) * Q[ii] # Update H = X.T X w
                    daxpy(n_features, w[ii] - w_ii, Q_ptr + ii * n_features, 1,
                          <DOUBLE*> &H[0], 1)

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

                gap = gram_enet_duality_gap(n_features, n_tasks, <DOUBLE*> &q[0],
                                       y_norm2, <DOUBLE*> &H[0],
                                       <DOUBLE*> &w[0], <DOUBLE*> &XtA[0],
                                       &dual_norm_XtA,
                                       alpha, beta, positive)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return np.asarray(w), gap, tol, n_iter + 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_multi_task(double[::1, :] W, double l1_reg,
                                       double l2_reg, double[::1, :] X,
                                       double[:, :] Y, int max_iter,
                                       double tol, object rng,
                                       bint random=0):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net mult-task regression

        We minimize

        1 norm(y - X w, 2)^2 + l1_reg ||w||_21 + l2_reg norm(w, 2)^2
        -                                       ----
        2                                        2

    """
    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef unsigned int n_tasks = Y.shape[1]

    # to store XtA
    cdef double[:, ::1] XtA = np.zeros((n_features, n_tasks))
    cdef double XtA_axis1norm
    cdef double dual_norm_XtA

    # initial value of the residuals
    cdef double[:, ::1] R = np.zeros((n_samples, n_tasks))

    cdef double[:] norm_cols_X = np.zeros(n_features)
    cdef double[::1] tmp = np.zeros(n_tasks, dtype=np.float)
    cdef double[:] w_ii = np.zeros(n_tasks, dtype=np.float)
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double nn
    cdef double W_ii_abs_max
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double ry_sum
    cdef double l21_norm
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef double* X_ptr = &X[0, 0]
    cdef double* W_ptr = &W[0, 0]
    cdef double* Y_ptr = &Y[0, 0]
    cdef double* wii_ptr = &w_ii[0]

    if l1_reg == 0:
        warnings.warn("Coordinate descent with l1_reg=0 may lead to unexpected"
            " results and is discouraged.")

    with nogil:
        # norm_cols_X = (np.asarray(X) ** 2).sum(axis=0)
        for ii in range(n_features):
            for jj in range(n_samples):
                norm_cols_X[ii] += X[jj, ii] ** 2

        # R = Y - np.dot(X, W.T)
        for ii in range(n_samples):
            for jj in range(n_tasks):
                R[ii, jj] = Y[ii, jj] - (
                    ddot(n_features, X_ptr + ii, n_samples, W_ptr + jj, n_tasks)
                    )

        # tol = tol * linalg.norm(Y, ord='fro') ** 2
        tol = tol * dnrm2(n_samples * n_tasks, Y_ptr, 1) ** 2

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                # w_ii = W[:, ii] # Store previous value
                dcopy(n_tasks, W_ptr + ii * n_tasks, 1, wii_ptr, 1)

                # if np.sum(w_ii ** 2) != 0.0:  # can do better
                if dnrm2(n_tasks, wii_ptr, 1) != 0.0:
                    # R += np.dot(X[:, ii][:, None], w_ii[None, :]) # rank 1 update
                    dger(CblasRowMajor, n_samples, n_tasks, 1.0,
                         X_ptr + ii * n_samples, 1,
                         wii_ptr, 1, &R[0, 0], n_tasks)

                # tmp = np.dot(X[:, ii][None, :], R).ravel()
                dgemv(CblasRowMajor, CblasTrans,
                      n_samples, n_tasks, 1.0, &R[0, 0], n_tasks,
                      X_ptr + ii * n_samples, 1, 0.0, &tmp[0], 1)

                # nn = sqrt(np.sum(tmp ** 2))
                nn = dnrm2(n_tasks, &tmp[0], 1)

                # W[:, ii] = tmp * fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg)
                dcopy(n_tasks, &tmp[0], 1, W_ptr + ii * n_tasks, 1)
                dscal(n_tasks, fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg),
                          W_ptr + ii * n_tasks, 1)

                # if np.sum(W[:, ii] ** 2) != 0.0:  # can do better
                if dnrm2(n_tasks, W_ptr + ii * n_tasks, 1) != 0.0:
                    # R -= np.dot(X[:, ii][:, None], W[:, ii][None, :])
                    # Update residual : rank 1 update
                    dger(CblasRowMajor, n_samples, n_tasks, -1.0,
                         X_ptr + ii * n_samples, 1, W_ptr + ii * n_tasks, 1,
                         &R[0, 0], n_tasks)

                # update the maximum absolute coefficient update
                d_w_ii = diff_abs_max(n_tasks, W_ptr + ii * n_tasks, wii_ptr)

                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                W_ii_abs_max = abs_max(n_tasks, W_ptr + ii * n_tasks)
                if W_ii_abs_max > w_max:
                    w_max = W_ii_abs_max

            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                # XtA = np.dot(X.T, R) - l2_reg * W.T
                for ii in range(n_features):
                    for jj in range(n_tasks):
                        XtA[ii, jj] = ddot(
                            n_samples, X_ptr + ii * n_samples, 1,
                            &R[0, 0] + jj, n_tasks
                            ) - l2_reg * W[jj, ii]

                # dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))
                dual_norm_XtA = 0.0
                for ii in range(n_features):
                    # np.sqrt(np.sum(XtA ** 2, axis=1))
                    XtA_axis1norm = dnrm2(n_tasks, &XtA[0, 0] + ii * n_tasks, 1)
                    if XtA_axis1norm > dual_norm_XtA:
                        dual_norm_XtA = XtA_axis1norm

                # TODO: use squared L2 norm directly
                # R_norm = linalg.norm(R, ord='fro')
                # w_norm = linalg.norm(W, ord='fro')
                R_norm = dnrm2(n_samples * n_tasks, &R[0, 0], 1)
                w_norm = dnrm2(n_features * n_tasks, W_ptr, 1)
                if (dual_norm_XtA > l1_reg):
                    const =  l1_reg / dual_norm_XtA
                    A_norm = R_norm * const
                    gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
                else:
                    const = 1.0
                    gap = R_norm ** 2

                # ry_sum = np.sum(R * y)
                ry_sum = 0.0
                for ii in range(n_samples):
                    for jj in range(n_tasks):
                        ry_sum += R[ii, jj] * Y[ii, jj]

                # l21_norm = np.sqrt(np.sum(W ** 2, axis=0)).sum()
                l21_norm = 0.0
                for ii in range(n_features):
                    # np.sqrt(np.sum(W ** 2, axis=0))
                    l21_norm += dnrm2(n_tasks, W_ptr + n_tasks * ii, 1)

                gap += l1_reg * l21_norm - const * ry_sum + \
                     0.5 * l2_reg * (1 + const ** 2) * (w_norm ** 2)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return np.asarray(W), gap, tol, n_iter + 1
