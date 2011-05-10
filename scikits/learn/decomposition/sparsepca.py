import time
import sys

import numpy as np
from scipy import linalg
from scikits.learn.linear_model import Lasso, lars_path
from scikits.learn.externals.joblib import Parallel, delayed


##################################
# Utilities to spread load on CPUs
def _gen_even_slices(n, n_packs):
    """Generator to create n_packs slices going up to n.

    Examples
    ========

    >>> list(_gen_even_slices(10, 1))
    [slice(0, 10, None)]
    >>> list(_gen_even_slices(10, 10)) #doctest: +ELLIPSIS
    [slice(0, 1, None), slice(1, 2, None), ..., slice(9, 10, None)]
    >>> list(_gen_even_slices(10, 5)) #doctest: +ELLIPSIS
    [slice(0, 2, None), slice(2, 4, None), ..., slice(8, 10, None)]
    >>> list(_gen_even_slices(10, 3))
    [slice(0, 4, None), slice(4, 7, None), slice(7, 10, None)]

    """
    start = 0
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            yield slice(start, end, None)
            start = end


def cpu_count():
    """ Return the number of CPUs.
    """
    # XXX: should be in joblib
    try:
        import multiprocessing
    except ImportError:
        return 1
    return multiprocessing.cpu_count()


###########
# sparsePCA
def _update_V(U, Y, V, alpha, Gram=None, method='lars', tol=1e-8):
    """ Update V in sparse_pca loop.
    """
    coef = np.empty_like(V)
    if method == 'lars':
        err_mgt = np.seterr()
        np.seterr(all='ignore')
        n_samples = U.shape[0]
        #alpha = alpha * n_samples
        XY = np.dot(U.T, Y)
        for k in range(V.shape[1]):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            _, _, coef_path_ = lars_path(U, Y[:, k], Xy=XY[:, k], Gram=Gram,
                                         alpha_min=alpha, method='lasso')
            coef[:, k] = coef_path_[:, -1]
        np.seterr(**err_mgt)
    else:
        clf = Lasso(alpha=alpha, fit_intercept=False)
        for k in range(V.shape[1]):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            clf.coef_ = V[:, k]  # Init with previous value of Vk
            clf.fit(U, Y[:, k], max_iter=1000, tol=tol)
            coef[:, k] = clf.coef_
    return coef


def _update_U(U, Y, V, verbose=False, return_r2=False):
    """ Update U in sparse_pca loop. This function modifies in-place U.
    """
    n_atoms = len(V)
    n_samples = Y.shape[0]
    R = -np.dot(U, V)  # Residuals, computed 'in-place' for efficiency
    R += Y
    R = np.asfortranarray(R)
    ger, = linalg.get_blas_funcs(('ger',), (U, V))
    for k in xrange(n_atoms):
        # R += u_k v_k^T
        R = ger(1.0, U[:, k], V[k, :], a=R, overwrite_a=True)
        U[:, k] = np.dot(R, V[k, :].T)
        # Scale Uk
        if (U[:, k] ** 2).sum() < 1e-20:
            if verbose == 1:
                sys.stdout.write("+")
                sys.stdout.flush()
            elif verbose:
                print "Adding new random atom"
            U[:, k] = np.random.randn(n_samples)
            # Setting corresponding coefs to 0
            V[k, :] = 0.0
            U[:, k] /= linalg.norm(U[:, k])
        else:
            U[:, k] /= linalg.norm(U[:, k])
            R = ger(-1.0, U[:, k], V[k, :], a=R, overwrite_a=True)
    if return_r2:
        R **= 2
        R = np.sum(R)
        return U, R
    return U


def sparse_pca(Y, n_atoms, alpha, maxit=100, tol=1e-8, method='lars',
        n_jobs=1, U_init=None, V_init=None, callback=None, verbose=False):
    """
    Compute sparse PCA with n_atoms components

    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || V ||_1
                 (U,V)
                with || U_k ||_2 = 1 for all  0<= k < n_atoms
    """
    t0 = time.time()
    n_samples = Y.shape[0]
    # Avoid integer division problems
    alpha = float(alpha)

    if n_jobs == -1:
        n_jobs = cpu_count()

    # Init U and V with SVD of Y
    if U_init is not None and V_init is not None:
        U = np.asarray(U_init).copy()
        V = np.asarray(V_init).copy()
    else:
        U, S, V = linalg.svd(Y, full_matrices=False)
        V = S[:, np.newaxis] * V
    U = U[:, :n_atoms]
    V = V[:n_atoms, :]

    def cost_function():
        return 0.5 * np.sum((Y - np.dot(U, V)) ** 2) \
                    + alpha * np.sum(np.abs(V))

    E = []
    current_cost = np.nan

    if verbose == 1:
        print '[sparse_pca]',

    for ii in xrange(maxit):
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            print ("Iteration % 3i "
                "(elapsed time: % 3is, % 4.1fmn, current cost % 7.3f)" %
                    (ii, dt, dt / 60, current_cost))

        # Update V
        Gram = np.dot(U.T, U)
        slices = list(_gen_even_slices(V.shape[1], n_jobs))
        V_views = Parallel(n_jobs=n_jobs)(
                    delayed(_update_V)(U, Y[:, this_slice], V[:, this_slice],
                                    alpha=alpha / n_samples,
                                    Gram=Gram, method=method, tol=tol)
                    for this_slice in slices)
        for this_slice, this_view in zip(slices, V_views):
            V[:, this_slice] = this_view

        # Update U
        U = _update_U(U, Y, V, verbose=verbose)

        current_cost = cost_function()
        E.append(current_cost)

        if ii > 0:
            dE = E[-2] - E[-1]
            assert(dE >= -tol * E[-1])
            if dE < tol * E[-1]:
                if verbose == 1:
                    # A line return
                    print ""
                elif verbose:
                    print "--- Convergence reached after %d iterations" % ii
                break
        if ii % 5 == 0 and callback is not None:
            callback(locals())

    return U, V, E


if __name__ == '__main__':

    # Generate toy data
    n_atoms = 3
    n_samples = 100
    img_sz = (10, 10)
    n_features = img_sz[0] * img_sz[1]

    np.random.seed(0)
    U = np.random.randn(n_samples, n_atoms)
    V = np.random.randn(n_atoms, n_features)

    centers = [(3, 3), (6, 7), (8, 1)]
    sz = [1, 2, 1]
    for k in range(n_atoms):
        img = np.zeros(img_sz)
        xmin, xmax = centers[k][0] - sz[k], centers[k][0] + sz[k]
        ymin, ymax = centers[k][1] - sz[k], centers[k][1] + sz[k]
        img[xmin:xmax][:, ymin:ymax] = 1.0
        V[k, :] = img.ravel()

    # Y is defined by : Y = UV + noise
    Y = np.dot(U, V)
    Y += 0.1 * np.random.randn(Y.shape[0], Y.shape[1])  # Add noise

    # Estimate U,V
    alpha = 0.5
    U_estimated, V_estimated, E = sparse_pca(Y, n_atoms, alpha, maxit=100,
                                             method='lasso', n_jobs=1,
                                             verbose=2, tol=1e-18)

    # View results
    import pylab as pl
    pl.close('all')

    for k in range(n_atoms):
        pl.matshow(np.reshape(V_estimated[k, :], img_sz))
        pl.title('Atom %d' % k)
        pl.colorbar()

    pl.figure()
    pl.plot(E)
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.show()
