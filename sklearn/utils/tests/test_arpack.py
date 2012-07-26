"""
Test ARPACK.  This is copied from the scipy test suite
"""

__usage__ = """
To run tests locally:
  python tests/test_arpack.py [-l<int>] [-v<int>]

"""

import warnings

import numpy as np

from numpy.testing import assert_allclose, \
        assert_array_almost_equal_nulp, TestCase, run_module_suite, dec, \
        assert_raises, verbose, assert_equal

from numpy import array, finfo, argsort, dot, round, conj, random
from scipy.linalg import eig, eigh
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix, isspmatrix
from scipy.sparse.linalg import LinearOperator, aslinearoperator

# explicitly import the scikit-learn ports for testing
from sklearn.utils.arpack import _eigs, _eigsh, _svds, \
     ArpackNoConvergence
(eigs, eigsh, svds) = (_eigs, _eigsh, _svds)

from scipy.linalg import svd


# eigs() and eigsh() are called many times, so apply a filter for the warnings
# they generate here.
_eigs_warn_msg = "Single-precision types in `eigs` and `eighs`"

def setup_module():
    warnings.filterwarnings("ignore", message=_eigs_warn_msg)

def teardown_module():
    warnings.filterwarnings("default", message=_eigs_warn_msg)


# precision for tests
_ndigits = {'f': 3, 'd': 11, 'F': 3, 'D': 11}

def _get_test_tolerance(type_char, mattype=None):
    """
    Return tolerance values suitable for a given test:

    Parameters
    ----------
    type_char : {'f', 'd', 'F', 'D'}
        Data type in ARPACK eigenvalue problem
    mattype : {csr_matrix, aslinearoperator, asarray}, optional
        Linear operator type

    Returns
    -------
    tol
        Tolerance to pass to the ARPACK routine
    rtol
        Relative tolerance for outputs
    atol
        Absolute tolerance for outputs

    """

    rtol = {'f': 3000 * np.finfo(np.float32).eps,
            'F': 3000 * np.finfo(np.float32).eps,
            'd': 2000 * np.finfo(np.float64).eps,
            'D': 2000 * np.finfo(np.float64).eps}[type_char]
    atol = rtol
    tol = 0

    if mattype is aslinearoperator and type_char in ('f', 'F'):
        # iterative methods in single precision: worse errors
        # also: bump ARPACK tolerance so that the iterative method converges
        tol = 30 * np.finfo(np.float32).eps
        rtol *= 5

    if mattype is csr_matrix and type_char in ('f', 'F'):
        # sparse in single precision: worse errors
        rtol *= 5

    return tol, rtol, atol

def generate_matrix(N, complex=False, hermitian=False,
                    pos_definite=False, sparse=False):
    M = np.random.random((N,N))
    if complex:
        M = M + 1j * np.random.random((N,N))

    if hermitian:
        if pos_definite:
            if sparse:
                i = np.arange(N)
                j = np.random.randint(N, size=N-2)
                i, j = np.meshgrid(i, j)
                M[i,j] = 0
            M = np.dot(M.conj(), M.T)
        else:
            M = np.dot(M.conj(), M.T)
            if sparse:
                i = np.random.randint(N, size=N * N / 4)
                j = np.random.randint(N, size=N * N / 4)
                ind = np.where(i == j)
                j[ind] = (j[ind] + 1) % N
                M[i,j] = 0
                M[j,i] = 0
    else:
        if sparse:
            i = np.random.randint(N, size=N * N / 2)
            j = np.random.randint(N, size=N * N / 2)
            M[i,j] = 0
    return M


def _aslinearoperator_with_dtype(m):
    m = aslinearoperator(m)
    if not hasattr(m, 'dtype'):
        x = np.zeros(m.shape[1])
        m.dtype = (m * x).dtype
    return m


def assert_allclose_cc(actual, desired, **kw):
    """Almost equal or complex conjugates almost equal"""
    try:
        assert_allclose(actual, desired, **kw)
    except:
        assert_allclose(actual, conj(desired), **kw)


def argsort_which(eval, typ, k, which,
                  sigma=None, OPpart=None, mode=None):
    """Return sorted indices of eigenvalues using the "which" keyword
    from eigs and eigsh"""
    if sigma is None:
        reval = np.round(eval, decimals=_ndigits[typ])
    else:
        if mode is None or mode=='normal':
            if OPpart is None:
                reval = 1. / (eval - sigma)
            elif OPpart == 'r':
                reval = 0.5 * (1. / (eval - sigma)
                               + 1. / (eval - np.conj(sigma)))
            elif OPpart == 'i':
                reval = -0.5j * (1. / (eval - sigma)
                                 - 1. / (eval - np.conj(sigma)))
        elif mode=='cayley':
            reval = (eval + sigma) / (eval - sigma)
        elif mode=='buckling':
            reval = eval / (eval - sigma)
        else:
            raise ValueError("mode='%s' not recognized" % mode)

        reval = np.round(reval, decimals=_ndigits[typ])

    if which in ['LM', 'SM']:
        ind = np.argsort(abs(reval))
    elif which in ['LR', 'SR', 'LA', 'SA', 'BE']:
        ind = np.argsort(np.real(reval))
    elif which in ['LI', 'SI']:
        # for LI,SI ARPACK returns largest,smallest abs(imaginary) why?
        if typ.islower():
            ind = np.argsort(abs(np.imag(reval)))
        else:
            ind = np.argsort(np.imag(reval))
    else:
        raise ValueError("which='%s' is unrecognized" % which)

    if which in ['LM', 'LA', 'LR', 'LI']:
        return ind[-k:]
    elif which in ['SM', 'SA', 'SR', 'SI']:
        return ind[:k]
    elif which == 'BE':
        return np.concatenate((ind[:k/2], ind[k/2-k:]))


def eval_evec(symmetric, d, typ, k, which, v0=None, sigma=None,
              mattype=np.asarray, OPpart=None, mode='normal'):
    general = ('bmat' in d)

    if symmetric:
        eigs_func = eigsh
    else:
        eigs_func = eigs

    if general:
        err = ("error for %s:general, typ=%s, which=%s, sigma=%s, "
               "mattype=%s, OPpart=%s, mode=%s" % (eigs_func.__name__,
                                                   typ, which, sigma,
                                                   mattype.__name__,
                                                   OPpart, mode))
    else:
        err = ("error for %s:standard, typ=%s, which=%s, sigma=%s, "
               "mattype=%s, OPpart=%s, mode=%s" % (eigs_func.__name__,
                                                   typ, which, sigma,
                                                   mattype.__name__,
                                                   OPpart, mode))

    a = d['mat'].astype(typ)
    ac = mattype(a)

    if general:
        b = d['bmat'].astype(typ.lower())
        bc = mattype(b)

    # get exact eigenvalues
    exact_eval = d['eval'].astype(typ.upper())
    ind = argsort_which(exact_eval, typ, k, which,
                        sigma, OPpart, mode)
    exact_eval_a = exact_eval
    exact_eval = exact_eval[ind]

    # compute arpack eigenvalues
    kwargs = dict(which=which, v0=v0, sigma=sigma)
    if eigs_func is eigsh:
        kwargs['mode'] = mode
    else:
        kwargs['OPpart'] = OPpart

    # compute suitable tolerances
    kwargs['tol'], rtol, atol = _get_test_tolerance(typ, mattype)

    # solve
    if general:
        try:
            eval, evec = eigs_func(ac, k, bc, **kwargs)
        except ArpackNoConvergence:
            kwargs['maxiter'] = 20*a.shape[0]
            eval, evec = eigs_func(ac, k, bc, **kwargs)
    else:
        try:
            eval, evec = eigs_func(ac, k, **kwargs)
        except ArpackNoConvergence:
            kwargs['maxiter'] = 20*a.shape[0]
            eval, evec = eigs_func(ac, k, **kwargs)

    ind = argsort_which(eval, typ, k, which,
                        sigma, OPpart, mode)
    eval_a = eval
    eval = eval[ind]
    evec = evec[:,ind]

    # check eigenvalues
    assert_allclose_cc(eval, exact_eval, rtol=rtol, atol=atol, err_msg=err)

    # check eigenvectors
    LHS = np.dot(a, evec)
    if general:
        RHS = eval * np.dot(b, evec)
    else:
        RHS = eval * evec

    assert_allclose(LHS, RHS, rtol=rtol, atol=atol, err_msg=err)

class DictWithRepr(dict):
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return "<%s>" % self.name

class SymmetricParams:
    def __init__(self):
        self.eigs = eigsh
        self.which = ['LM', 'SM', 'LA', 'SA', 'BE']
        self.mattypes = [csr_matrix, aslinearoperator, np.asarray]
        self.sigmas_modes = {None : ['normal'],
                             0.5 : ['normal', 'buckling', 'cayley']}

        #generate matrices
        # these should all be float32 so that the eigenvalues
        # are the same in float32 and float64
        N = 6
        np.random.seed(2300)
        Ar = generate_matrix(N, hermitian=True,
                             pos_definite=True).astype('f').astype('d')
        M = generate_matrix(N, hermitian=True,
                            pos_definite=True).astype('f').astype('d')
        Ac = generate_matrix(N, hermitian=True, pos_definite=True,
                             complex=True).astype('F').astype('D')
        v0 = np.random.random(N)

        # standard symmetric problem
        SS = DictWithRepr("std-symmetric")
        SS['mat'] = Ar
        SS['v0'] = v0
        SS['eval'] = eigh(SS['mat'], eigvals_only=True)

        # general symmetric problem
        GS = DictWithRepr("gen-symmetric")
        GS['mat'] = Ar
        GS['bmat'] = M
        GS['v0'] = v0
        GS['eval'] = eigh(GS['mat'], GS['bmat'], eigvals_only=True)

        # standard hermitian problem
        SH = DictWithRepr("std-hermitian")
        SH['mat'] = Ac
        SH['v0'] = v0
        SH['eval'] = eigh(SH['mat'], eigvals_only=True)

        # general hermitian problem
        GH = DictWithRepr("gen-hermitian")
        GH['mat'] = Ac
        GH['bmat'] = M
        GH['v0'] = v0
        GH['eval'] = eigh(GH['mat'], GH['bmat'], eigvals_only=True)

        self.real_test_cases = [SS, GS]
        self.complex_test_cases = [SH, GH]

class NonSymmetricParams:
    def __init__(self):
        self.eigs = eigs
        self.which = ['LM', 'LR', 'LI']#, 'SM', 'LR', 'SR', 'LI', 'SI']
        self.mattypes = [csr_matrix, aslinearoperator, np.asarray]
        self.sigmas_OPparts = {None : [None],
                               0.1 : ['r'],
                               0.1 + 0.1j : ['r', 'i']}

        #generate matrices
        # these should all be float32 so that the eigenvalues
        # are the same in float32 and float64
        N = 6
        np.random.seed(2300)
        Ar = generate_matrix(N).astype('f').astype('d')
        M = generate_matrix(N, hermitian=True,
                            pos_definite=True).astype('f').astype('d')
        Ac = generate_matrix(N, complex=True).astype('F').astype('D')
        v0 = np.random.random(N)

        # standard real nonsymmetric problem
        SNR = DictWithRepr("std-real-nonsym")
        SNR['mat'] = Ar
        SNR['v0'] = v0
        SNR['eval'] = eig(SNR['mat'], left=False, right=False)

        # general real nonsymmetric problem
        GNR = DictWithRepr("gen-real-nonsym")
        GNR['mat'] = Ar
        GNR['bmat'] = M
        GNR['v0'] = v0
        GNR['eval'] = eig(GNR['mat'], GNR['bmat'], left=False, right=False)

        # standard complex nonsymmetric problem
        SNC = DictWithRepr("std-cmplx-nonsym")
        SNC['mat'] = Ac
        SNC['v0'] = v0
        SNC['eval'] = eig(SNC['mat'], left=False, right=False)

        # general complex nonsymmetric problem
        GNC = DictWithRepr("gen-cmplx-nonsym")
        GNC['mat'] = Ac
        GNC['bmat'] = M
        GNC['v0'] = v0
        GNC['eval'] = eig(GNC['mat'], GNC['bmat'], left=False, right=False)

        self.real_test_cases = [SNR, GNR]
        self.complex_test_cases = [SNC, GNC]


def test_symmetric_modes():
    params = SymmetricParams()
    k = 2
    symmetric = True
    for D in params.real_test_cases:
        for typ in 'fd':
            for which in params.which:
                for mattype in params.mattypes:
                    for (sigma, modes) in params.sigmas_modes.iteritems():
                        for mode in modes:
                            yield  (eval_evec, symmetric, D, typ, k, which,
                                    None, sigma, mattype, None, mode)


def test_hermitian_modes():
    params = SymmetricParams()
    k = 2
    symmetric = True
    for D in params.complex_test_cases:
        for typ in 'FD':
            for which in params.which:
                if which == 'BE': continue  # BE invalid for complex
                for mattype in params.mattypes:
                    for sigma in params.sigmas_modes:
                        yield  (eval_evec, symmetric, D, typ, k, which,
                                None, sigma, mattype)


def test_symmetric_starting_vector():
    params = SymmetricParams()
    symmetric = True
    for k in [1, 2, 3, 4, 5]:
        for D in params.real_test_cases:
            for typ in 'fd':
                v0 = random.rand(len(D['v0'])).astype(typ)
                yield (eval_evec, symmetric, D, typ, k, 'LM', v0)


def test_symmetric_no_convergence():
    np.random.seed(1234)
    m = generate_matrix(30, hermitian=True, pos_definite=True)
    tol, rtol, atol = _get_test_tolerance('d')
    try:
        w, v = eigsh(m, 4, which='LM', v0=m[:, 0], maxiter=5, tol=tol)
        raise AssertionError("Spurious no-error exit")
    except ArpackNoConvergence, err:
        k = len(err.eigenvalues)
        if k <= 0:
            raise AssertionError("Spurious no-eigenvalues-found case")
        w, v = err.eigenvalues, err.eigenvectors
        assert_allclose(dot(m, v), w * v, rtol=rtol, atol=atol)


def test_real_nonsymmetric_modes():
    params = NonSymmetricParams()
    k = 2
    symmetric = False
    for D in params.real_test_cases:
        for typ in 'fd':
            for which in params.which:
                for mattype in params.mattypes:
                    for sigma, OPparts in params.sigmas_OPparts.iteritems():
                        for OPpart in OPparts:
                            yield (eval_evec, symmetric, D, typ, k, which,
                                   None, sigma, mattype, OPpart)


def test_complex_nonsymmetric_modes():
    params = NonSymmetricParams()
    k = 2
    symmetric = False
    for D in params.complex_test_cases:
        for typ in 'DF':
            for which in params.which:
                for mattype in params.mattypes:
                    for sigma in params.sigmas_OPparts:
                        yield (eval_evec, symmetric, D, typ, k, which,
                               None, sigma, mattype)


def test_standard_nonsymmetric_starting_vector():
    params = NonSymmetricParams()
    sigma = None
    symmetric = False
    for k in [1, 2, 3, 4]:
        for d in params.complex_test_cases:
            for typ in 'FD':
                A = d['mat']
                n = A.shape[0]
                v0 = random.rand(n).astype(typ)
                yield (eval_evec, symmetric, d, typ, k, "LM", v0, sigma)


def test_general_nonsymmetric_starting_vector():
    params = NonSymmetricParams()
    sigma = None
    symmetric = False
    for k in [1, 2, 3, 4]:
        for d in params.complex_test_cases:
            for typ in 'FD':
                A = d['mat']
                n = A.shape[0]
                v0 = random.rand(n).astype(typ)
                yield (eval_evec, symmetric, d, typ, k, "LM", v0, sigma)


def test_standard_nonsymmetric_no_convergence():
    np.random.seed(1234)
    m = generate_matrix(30, complex=True)
    tol, rtol, atol = _get_test_tolerance('d')
    try:
        w, v = eigs(m, 4, which='LM', v0=m[:, 0], maxiter=5, tol=tol)
        raise AssertionError("Spurious no-error exit")
    except ArpackNoConvergence, err:
        k = len(err.eigenvalues)
        if k <= 0:
            raise AssertionError("Spurious no-eigenvalues-found case")
        w, v = err.eigenvalues, err.eigenvectors
        for ww, vv in zip(w, v.T):
            assert_allclose(dot(m, vv), ww * vv, rtol=rtol, atol=atol)


def test_eigen_bad_shapes():
    # A is not square.
    A = csc_matrix(np.zeros((2, 3)))
    assert_raises(ValueError, eigs, A)


def test_eigen_bad_kwargs():
    # Test eigen on wrong keyword argument
    A = csc_matrix(np.zeros((2, 2)))
    assert_raises(ValueError, eigs, A, which='XX')

# Add underscores so this test will not be run.
# It will fail fatally unless we're using scipy 0.10.1 or greater.
def __test_ticket_1459_arpack_crash():
    for dtype in [np.float32, np.float64]:
        # XXX: this test does not seem to catch the issue for float32,
        #      but we made the same fix there, just to be sure

        N = 6
        k = 2

        np.random.seed(2301)
        A = np.random.random((N, N)).astype(dtype)
        v0 = np.array([-0.71063568258907849895, -0.83185111795729227424,
                       -0.34365925382227402451, 0.46122533684552280420,
                       -0.58001341115969040629, -0.78844877570084292984e-01],
                      dtype=dtype)

        # Should not crash:
        evals, evecs = eigs(A, k, v0=v0)


#----------------------------------------------------------------------
# sparse SVD tests

def sorted_svd(m, k):
    #Compute svd of a dense matrix m, and return singular vectors/values
    #sorted.
    if isspmatrix(m):
        m = m.todense()
    u, s, vh = svd(m)
    ii = np.argsort(s)[-k:]

    return u[:, ii], s[ii], vh[ii]


def svd_estimate(u, s, vh):
    return np.dot(u, np.dot(np.diag(s), vh))


def test_svd_simple_real():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1, 0, 2],
                  [0, 0, 1]], np.float)
    y = np.array([[1, 2, 3, 8],
                  [3, 4, 3, 5],
                  [1, 0, 2, 3],
                  [0, 0, 1, 0]], np.float)
    z = csc_matrix(x)

    for m in [x.T, x, y, z, z.T]:
        for k in range(1, min(m.shape)):
            u, s, vh = sorted_svd(m, k)
            su, ss, svh = svds(m, k)

            m_hat = svd_estimate(u, s, vh)
            sm_hat = svd_estimate(su, ss, svh)

            assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)


def test_svd_simple_complex():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1 + 1j, 0, 2],
                  [0, 0, 1]], np.complex)
    y = np.array([[1, 2, 3, 8 + 5j],
                  [3 - 2j, 4, 3, 5],
                  [1, 0, 2, 3],
                  [0, 0, 1, 0]], np.complex)
    z = csc_matrix(x)

    for m in [x, x.T.conjugate(), x.T, y, y.conjugate(), z, z.T]:
        for k in range(1, min(m.shape) - 1):
            u, s, vh = sorted_svd(m, k)
            su, ss, svh = svds(m, k)

            m_hat = svd_estimate(u, s, vh)
            sm_hat = svd_estimate(su, ss, svh)

            assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)


if __name__ == "__main__":
    run_module_suite()
