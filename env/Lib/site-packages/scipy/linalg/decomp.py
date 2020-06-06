#
# Author: Pearu Peterson, March 2002
#
# additions by Travis Oliphant, March 2002
# additions by Eric Jones,      June 2002
# additions by Johannes Loehnert, June 2006
# additions by Bart Vandereycken, June 2006
# additions by Andrew D Straw, May 2007
# additions by Tiziano Zito, November 2008
#
# April 2010: Functions for LU, QR, SVD, Schur and Cholesky decompositions were
# moved to their own files.  Still in this file are functions for eigenstuff
# and for the Hessenberg form.

from __future__ import division, print_function, absolute_import

__all__ = ['eig', 'eigvals', 'eigh', 'eigvalsh',
           'eig_banded', 'eigvals_banded',
           'eigh_tridiagonal', 'eigvalsh_tridiagonal', 'hessenberg', 'cdf2rdf']

import numpy
from numpy import (array, isfinite, inexact, nonzero, iscomplexobj, cast,
                   flatnonzero, conj, asarray, argsort, empty, newaxis,
                   argwhere, iscomplex, eye, zeros, einsum)
# Local imports
from scipy._lib.six import xrange
from scipy._lib._util import _asarray_validated
from scipy._lib.six import string_types
from .misc import LinAlgError, _datacopied, norm
from .lapack import get_lapack_funcs, _compute_lwork


_I = cast['F'](1j)


def _make_complex_eigvecs(w, vin, dtype):
    """
    Produce complex-valued eigenvectors from LAPACK DGGEV real-valued output
    """
    # - see LAPACK man page DGGEV at ALPHAI
    v = numpy.array(vin, dtype=dtype)
    m = (w.imag > 0)
    m[:-1] |= (w.imag[1:] < 0)  # workaround for LAPACK bug, cf. ticket #709
    for i in flatnonzero(m):
        v.imag[:, i] = vin[:, i+1]
        conj(v[:, i], v[:, i+1])
    return v


def _make_eigvals(alpha, beta, homogeneous_eigvals):
    if homogeneous_eigvals:
        if beta is None:
            return numpy.vstack((alpha, numpy.ones_like(alpha)))
        else:
            return numpy.vstack((alpha, beta))
    else:
        if beta is None:
            return alpha
        else:
            w = numpy.empty_like(alpha)
            alpha_zero = (alpha == 0)
            beta_zero = (beta == 0)
            beta_nonzero = ~beta_zero
            w[beta_nonzero] = alpha[beta_nonzero]/beta[beta_nonzero]
            # Use numpy.inf for complex values too since
            # 1/numpy.inf = 0, i.e. it correctly behaves as projective
            # infinity.
            w[~alpha_zero & beta_zero] = numpy.inf
            if numpy.all(alpha.imag == 0):
                w[alpha_zero & beta_zero] = numpy.nan
            else:
                w[alpha_zero & beta_zero] = complex(numpy.nan, numpy.nan)
            return w


def _geneig(a1, b1, left, right, overwrite_a, overwrite_b,
            homogeneous_eigvals):
    ggev, = get_lapack_funcs(('ggev',), (a1, b1))
    cvl, cvr = left, right
    res = ggev(a1, b1, lwork=-1)
    lwork = res[-2][0].real.astype(numpy.int)
    if ggev.typecode in 'cz':
        alpha, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr, lwork,
                                               overwrite_a, overwrite_b)
        w = _make_eigvals(alpha, beta, homogeneous_eigvals)
    else:
        alphar, alphai, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr,
                                                        lwork, overwrite_a,
                                                        overwrite_b)
        alpha = alphar + _I * alphai
        w = _make_eigvals(alpha, beta, homogeneous_eigvals)
    _check_info(info, 'generalized eig algorithm (ggev)')

    only_real = numpy.all(w.imag == 0.0)
    if not (ggev.typecode in 'cz' or only_real):
        t = w.dtype.char
        if left:
            vl = _make_complex_eigvecs(w, vl, t)
        if right:
            vr = _make_complex_eigvecs(w, vr, t)

    # the eigenvectors returned by the lapack function are NOT normalized
    for i in xrange(vr.shape[0]):
        if right:
            vr[:, i] /= norm(vr[:, i])
        if left:
            vl[:, i] /= norm(vl[:, i])

    if not (left or right):
        return w
    if left:
        if right:
            return w, vl, vr
        return w, vl
    return w, vr


def eig(a, b=None, left=False, right=True, overwrite_a=False,
        overwrite_b=False, check_finite=True, homogeneous_eigvals=False):
    """
    Solve an ordinary or generalized eigenvalue problem of a square matrix.

    Find eigenvalues w and right or left eigenvectors of a general matrix::

        a   vr[:,i] = w[i]        b   vr[:,i]
        a.H vl[:,i] = w[i].conj() b.H vl[:,i]

    where ``.H`` is the Hermitian conjugation.

    Parameters
    ----------
    a : (M, M) array_like
        A complex or real matrix whose eigenvalues and eigenvectors
        will be computed.
    b : (M, M) array_like, optional
        Right-hand side matrix in a generalized eigenvalue problem.
        Default is None, identity matrix is assumed.
    left : bool, optional
        Whether to calculate and return left eigenvectors.  Default is False.
    right : bool, optional
        Whether to calculate and return right eigenvectors.  Default is True.
    overwrite_a : bool, optional
        Whether to overwrite `a`; may improve performance.  Default is False.
    overwrite_b : bool, optional
        Whether to overwrite `b`; may improve performance.  Default is False.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    homogeneous_eigvals : bool, optional
        If True, return the eigenvalues in homogeneous coordinates.
        In this case ``w`` is a (2, M) array so that::

            w[1,i] a vr[:,i] = w[0,i] b vr[:,i]

        Default is False.

    Returns
    -------
    w : (M,) or (2, M) double or complex ndarray
        The eigenvalues, each repeated according to its
        multiplicity. The shape is (M,) unless
        ``homogeneous_eigvals=True``.
    vl : (M, M) double or complex ndarray
        The normalized left eigenvector corresponding to the eigenvalue
        ``w[i]`` is the column vl[:,i]. Only returned if ``left=True``.
    vr : (M, M) double or complex ndarray
        The normalized right eigenvector corresponding to the eigenvalue
        ``w[i]`` is the column ``vr[:,i]``.  Only returned if ``right=True``.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge.

    See Also
    --------
    eigvals : eigenvalues of general arrays
    eigh : Eigenvalues and right eigenvectors for symmetric/Hermitian arrays.
    eig_banded : eigenvalues and right eigenvectors for symmetric/Hermitian
        band matrices
    eigh_tridiagonal : eigenvalues and right eiegenvectors for
        symmetric/Hermitian tridiagonal matrices

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.array([[0., -1.], [1., 0.]])
    >>> linalg.eigvals(a)
    array([0.+1.j, 0.-1.j])

    >>> b = np.array([[0., 1.], [1., 1.]])
    >>> linalg.eigvals(a, b)
    array([ 1.+0.j, -1.+0.j])

    >>> a = np.array([[3., 0., 0.], [0., 8., 0.], [0., 0., 7.]])
    >>> linalg.eigvals(a, homogeneous_eigvals=True)
    array([[3.+0.j, 8.+0.j, 7.+0.j],
           [1.+0.j, 1.+0.j, 1.+0.j]])

    >>> a = np.array([[0., -1.], [1., 0.]])
    >>> linalg.eigvals(a) == linalg.eig(a)[0]
    array([ True,  True])
    >>> linalg.eig(a, left=True, right=False)[1] # normalized left eigenvector
    array([[-0.70710678+0.j        , -0.70710678-0.j        ],
           [-0.        +0.70710678j, -0.        -0.70710678j]])
    >>> linalg.eig(a, left=False, right=True)[1] # normalized right eigenvector
    array([[0.70710678+0.j        , 0.70710678-0.j        ],
           [0.        -0.70710678j, 0.        +0.70710678j]])



    """
    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    if b is not None:
        b1 = _asarray_validated(b, check_finite=check_finite)
        overwrite_b = overwrite_b or _datacopied(b1, b)
        if len(b1.shape) != 2 or b1.shape[0] != b1.shape[1]:
            raise ValueError('expected square matrix')
        if b1.shape != a1.shape:
            raise ValueError('a and b must have the same shape')
        return _geneig(a1, b1, left, right, overwrite_a, overwrite_b,
                       homogeneous_eigvals)

    geev, geev_lwork = get_lapack_funcs(('geev', 'geev_lwork'), (a1,))
    compute_vl, compute_vr = left, right

    lwork = _compute_lwork(geev_lwork, a1.shape[0],
                           compute_vl=compute_vl,
                           compute_vr=compute_vr)

    if geev.typecode in 'cz':
        w, vl, vr, info = geev(a1, lwork=lwork,
                               compute_vl=compute_vl,
                               compute_vr=compute_vr,
                               overwrite_a=overwrite_a)
        w = _make_eigvals(w, None, homogeneous_eigvals)
    else:
        wr, wi, vl, vr, info = geev(a1, lwork=lwork,
                                    compute_vl=compute_vl,
                                    compute_vr=compute_vr,
                                    overwrite_a=overwrite_a)
        t = {'f': 'F', 'd': 'D'}[wr.dtype.char]
        w = wr + _I * wi
        w = _make_eigvals(w, None, homogeneous_eigvals)

    _check_info(info, 'eig algorithm (geev)',
                positive='did not converge (only eigenvalues '
                         'with order >= %d have converged)')

    only_real = numpy.all(w.imag == 0.0)
    if not (geev.typecode in 'cz' or only_real):
        t = w.dtype.char
        if left:
            vl = _make_complex_eigvecs(w, vl, t)
        if right:
            vr = _make_complex_eigvecs(w, vr, t)
    if not (left or right):
        return w
    if left:
        if right:
            return w, vl, vr
        return w, vl
    return w, vr


def eigh(a, b=None, lower=True, eigvals_only=False, overwrite_a=False,
         overwrite_b=False, turbo=True, eigvals=None, type=1,
         check_finite=True):
    """
    Solve an ordinary or generalized eigenvalue problem for a complex
    Hermitian or real symmetric matrix.

    Find eigenvalues w and optionally eigenvectors v of matrix `a`, where
    `b` is positive definite::

                      a v[:,i] = w[i] b v[:,i]
        v[i,:].conj() a v[:,i] = w[i]
        v[i,:].conj() b v[:,i] = 1

    Parameters
    ----------
    a : (M, M) array_like
        A complex Hermitian or real symmetric matrix whose eigenvalues and
        eigenvectors will be computed.
    b : (M, M) array_like, optional
        A complex Hermitian or real symmetric definite positive matrix in.
        If omitted, identity matrix is assumed.
    lower : bool, optional
        Whether the pertinent array data is taken from the lower or upper
        triangle of `a`. (Default: lower)
    eigvals_only : bool, optional
        Whether to calculate only eigenvalues and no eigenvectors.
        (Default: both are calculated)
    turbo : bool, optional
        Use divide and conquer algorithm (faster but expensive in memory,
        only for generalized eigenvalue problem and if eigvals=None)
    eigvals : tuple (lo, hi), optional
        Indexes of the smallest and largest (in ascending order) eigenvalues
        and corresponding eigenvectors to be returned: 0 <= lo <= hi <= M-1.
        If omitted, all eigenvalues and eigenvectors are returned.
    type : int, optional
        Specifies the problem type to be solved:

           type = 1: a   v[:,i] = w[i] b v[:,i]

           type = 2: a b v[:,i] = w[i]   v[:,i]

           type = 3: b a v[:,i] = w[i]   v[:,i]
    overwrite_a : bool, optional
        Whether to overwrite data in `a` (may improve performance)
    overwrite_b : bool, optional
        Whether to overwrite data in `b` (may improve performance)
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    w : (N,) float ndarray
        The N (1<=N<=M) selected eigenvalues, in ascending order, each
        repeated according to its multiplicity.
    v : (M, N) complex ndarray
        (if eigvals_only == False)

        The normalized selected eigenvector corresponding to the
        eigenvalue w[i] is the column v[:,i].

        Normalization:

            type 1 and 3: v.conj() a      v  = w

            type 2: inv(v).conj() a  inv(v) = w

            type = 1 or 2: v.conj() b      v  = I

            type = 3: v.conj() inv(b) v  = I

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge,
        an error occurred, or b matrix is not definite positive. Note that
        if input matrices are not symmetric or hermitian, no error is reported
        but results will be wrong.

    See Also
    --------
    eigvalsh : eigenvalues of symmetric or Hermitian arrays
    eig : eigenvalues and right eigenvectors for non-symmetric arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eigh_tridiagonal : eigenvalues and right eiegenvectors for
        symmetric/Hermitian tridiagonal matrices

    Notes
    -----
    This function does not check the input array for being hermitian/symmetric
    in order to allow for representing arrays with only their upper/lower
    triangular parts.

    Examples
    --------
    >>> from scipy.linalg import eigh
    >>> A = np.array([[6, 3, 1, 5], [3, 0, 5, 1], [1, 5, 6, 2], [5, 1, 2, 2]])
    >>> w, v = eigh(A)
    >>> np.allclose(A @ v - v @ np.diag(w), np.zeros((4, 4)))
    True

    """
    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    if iscomplexobj(a1):
        cplx = True
    else:
        cplx = False
    if b is not None:
        b1 = _asarray_validated(b, check_finite=check_finite)
        overwrite_b = overwrite_b or _datacopied(b1, b)
        if len(b1.shape) != 2 or b1.shape[0] != b1.shape[1]:
            raise ValueError('expected square matrix')

        if b1.shape != a1.shape:
            raise ValueError("wrong b dimensions %s, should "
                             "be %s" % (str(b1.shape), str(a1.shape)))
        if iscomplexobj(b1):
            cplx = True
        else:
            cplx = cplx or False
    else:
        b1 = None

    # Set job for fortran routines
    _job = (eigvals_only and 'N') or 'V'

    # port eigenvalue range from python to fortran convention
    if eigvals is not None:
        lo, hi = eigvals
        if lo < 0 or hi >= a1.shape[0]:
            raise ValueError('The eigenvalue range specified is not valid.\n'
                             'Valid range is [%s,%s]' % (0, a1.shape[0]-1))
        lo += 1
        hi += 1
        eigvals = (lo, hi)

    # set lower
    if lower:
        uplo = 'L'
    else:
        uplo = 'U'

    # fix prefix for lapack routines
    if cplx:
        pfx = 'he'
    else:
        pfx = 'sy'

    #  Standard Eigenvalue Problem
    #  Use '*evr' routines
    # FIXME: implement calculation of optimal lwork
    #        for all lapack routines
    if b1 is None:
        driver = pfx+'evr'
        (evr,) = get_lapack_funcs((driver,), (a1,))
        if eigvals is None:
            w, v, info = evr(a1, uplo=uplo, jobz=_job, range="A", il=1,
                             iu=a1.shape[0], overwrite_a=overwrite_a)
        else:
            (lo, hi) = eigvals
            w_tot, v, info = evr(a1, uplo=uplo, jobz=_job, range="I",
                                 il=lo, iu=hi, overwrite_a=overwrite_a)
            w = w_tot[0:hi-lo+1]

    # Generalized Eigenvalue Problem
    else:
        # Use '*gvx' routines if range is specified
        if eigvals is not None:
            driver = pfx+'gvx'
            (gvx,) = get_lapack_funcs((driver,), (a1, b1))
            (lo, hi) = eigvals
            w_tot, v, ifail, info = gvx(a1, b1, uplo=uplo, iu=hi,
                                        itype=type, jobz=_job, il=lo,
                                        overwrite_a=overwrite_a,
                                        overwrite_b=overwrite_b)
            w = w_tot[0:hi-lo+1]
        # Use '*gvd' routine if turbo is on and no eigvals are specified
        elif turbo:
            driver = pfx+'gvd'
            (gvd,) = get_lapack_funcs((driver,), (a1, b1))
            v, w, info = gvd(a1, b1, uplo=uplo, itype=type, jobz=_job,
                             overwrite_a=overwrite_a,
                             overwrite_b=overwrite_b)
        # Use '*gv' routine if turbo is off and no eigvals are specified
        else:
            driver = pfx+'gv'
            (gv,) = get_lapack_funcs((driver,), (a1, b1))
            v, w, info = gv(a1, b1, uplo=uplo, itype=type, jobz=_job,
                            overwrite_a=overwrite_a,
                            overwrite_b=overwrite_b)

    # Check if we had a  successful exit
    if info == 0:
        if eigvals_only:
            return w
        else:
            return w, v
    _check_info(info, driver, positive=False)  # triage more specifically
    if info > 0 and b1 is None:
        raise LinAlgError("unrecoverable internal error.")

    # The algorithm failed to converge.
    elif 0 < info <= b1.shape[0]:
        if eigvals is not None:
            raise LinAlgError("the eigenvectors %s failed to"
                              " converge." % nonzero(ifail)-1)
        else:
            raise LinAlgError("internal fortran routine failed to converge: "
                              "%i off-diagonal elements of an "
                              "intermediate tridiagonal form did not converge"
                              " to zero." % info)

    # This occurs when b is not positive definite
    else:
        raise LinAlgError("the leading minor of order %i"
                          " of 'b' is not positive definite. The"
                          " factorization of 'b' could not be completed"
                          " and no eigenvalues or eigenvectors were"
                          " computed." % (info-b1.shape[0]))


_conv_dict = {0: 0, 1: 1, 2: 2,
              'all': 0, 'value': 1, 'index': 2,
              'a': 0, 'v': 1, 'i': 2}


def _check_select(select, select_range, max_ev, max_len):
    """Check that select is valid, convert to Fortran style."""
    if isinstance(select, string_types):
        select = select.lower()
    try:
        select = _conv_dict[select]
    except KeyError:
        raise ValueError('invalid argument for select')
    vl, vu = 0., 1.
    il = iu = 1
    if select != 0:  # (non-all)
        sr = asarray(select_range)
        if sr.ndim != 1 or sr.size != 2 or sr[1] < sr[0]:
            raise ValueError('select_range must be a 2-element array-like '
                             'in nondecreasing order')
        if select == 1:  # (value)
            vl, vu = sr
            if max_ev == 0:
                max_ev = max_len
        else:  # 2 (index)
            if sr.dtype.char.lower() not in 'hilqp':
                raise ValueError('when using select="i", select_range must '
                                 'contain integers, got dtype %s (%s)'
                                 % (sr.dtype, sr.dtype.char))
            # translate Python (0 ... N-1) into Fortran (1 ... N) with + 1
            il, iu = sr + 1
            if min(il, iu) < 1 or max(il, iu) > max_len:
                raise ValueError('select_range out of bounds')
            max_ev = iu - il + 1
    return select, vl, vu, il, iu, max_ev


def eig_banded(a_band, lower=False, eigvals_only=False, overwrite_a_band=False,
               select='a', select_range=None, max_ev=0, check_finite=True):
    """
    Solve real symmetric or complex hermitian band matrix eigenvalue problem.

    Find eigenvalues w and optionally right eigenvectors v of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    The matrix a is stored in a_band either in lower diagonal or upper
    diagonal ordered form:

        a_band[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        a_band[    i - j, j] == a[i,j]        (if lower form; i >= j)

    where u is the number of bands above the diagonal.

    Example of a_band (shape of a is (6,6), u=2)::

        upper form:
        *   *   a02 a13 a24 a35
        *   a01 a12 a23 a34 a45
        a00 a11 a22 a33 a44 a55

        lower form:
        a00 a11 a22 a33 a44 a55
        a10 a21 a32 a43 a54 *
        a20 a31 a42 a53 *   *

    Cells marked with * are not used.

    Parameters
    ----------
    a_band : (u+1, M) array_like
        The bands of the M by M matrix a.
    lower : bool, optional
        Is the matrix in the lower form. (Default is upper form)
    eigvals_only : bool, optional
        Compute only the eigenvalues and no eigenvectors.
        (Default: calculate also eigenvectors)
    overwrite_a_band : bool, optional
        Discard data in a_band (may enhance performance)
    select : {'a', 'v', 'i'}, optional
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max), optional
        Range of selected eigenvalues
    max_ev : int, optional
        For select=='v', maximum number of eigenvalues expected.
        For other values of select, has no meaning.

        In doubt, leave this parameter untouched.

    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    w : (M,) ndarray
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.
    v : (M, M) float or complex ndarray
        The normalized eigenvector corresponding to the eigenvalue w[i] is
        the column v[:,i].

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge.

    See Also
    --------
    eigvals_banded : eigenvalues for symmetric/Hermitian band matrices
    eig : eigenvalues and right eigenvectors of general arrays.
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eigh_tridiagonal : eigenvalues and right eiegenvectors for
        symmetric/Hermitian tridiagonal matrices

    Examples
    --------
    >>> from scipy.linalg import eig_banded
    >>> A = np.array([[1, 5, 2, 0], [5, 2, 5, 2], [2, 5, 3, 5], [0, 2, 5, 4]])
    >>> Ab = np.array([[1, 2, 3, 4], [5, 5, 5, 0], [2, 2, 0, 0]])
    >>> w, v = eig_banded(Ab, lower=True)
    >>> np.allclose(A @ v - v @ np.diag(w), np.zeros((4, 4)))
    True
    >>> w = eig_banded(Ab, lower=True, eigvals_only=True)
    >>> w
    array([-4.26200532, -2.22987175,  3.95222349, 12.53965359])

    Request only the eigenvalues between ``[-3, 4]``

    >>> w, v = eig_banded(Ab, lower=True, select='v', select_range=[-3, 4])
    >>> w
    array([-2.22987175,  3.95222349])

    """
    if eigvals_only or overwrite_a_band:
        a1 = _asarray_validated(a_band, check_finite=check_finite)
        overwrite_a_band = overwrite_a_band or (_datacopied(a1, a_band))
    else:
        a1 = array(a_band)
        if issubclass(a1.dtype.type, inexact) and not isfinite(a1).all():
            raise ValueError("array must not contain infs or NaNs")
        overwrite_a_band = 1

    if len(a1.shape) != 2:
        raise ValueError('expected two-dimensional array')
    select, vl, vu, il, iu, max_ev = _check_select(
        select, select_range, max_ev, a1.shape[1])
    del select_range
    if select == 0:
        if a1.dtype.char in 'GFD':
            # FIXME: implement this somewhen, for now go with builtin values
            # FIXME: calc optimal lwork by calling ?hbevd(lwork=-1)
            #        or by using calc_lwork.f ???
            # lwork = calc_lwork.hbevd(bevd.typecode, a1.shape[0], lower)
            internal_name = 'hbevd'
        else:  # a1.dtype.char in 'fd':
            # FIXME: implement this somewhen, for now go with builtin values
            #         see above
            # lwork = calc_lwork.sbevd(bevd.typecode, a1.shape[0], lower)
            internal_name = 'sbevd'
        bevd, = get_lapack_funcs((internal_name,), (a1,))
        w, v, info = bevd(a1, compute_v=not eigvals_only,
                          lower=lower, overwrite_ab=overwrite_a_band)
    else:  # select in [1, 2]
        if eigvals_only:
            max_ev = 1
        # calculate optimal abstol for dsbevx (see manpage)
        if a1.dtype.char in 'fF':  # single precision
            lamch, = get_lapack_funcs(('lamch',), (array(0, dtype='f'),))
        else:
            lamch, = get_lapack_funcs(('lamch',), (array(0, dtype='d'),))
        abstol = 2 * lamch('s')
        if a1.dtype.char in 'GFD':
            internal_name = 'hbevx'
        else:  # a1.dtype.char in 'gfd'
            internal_name = 'sbevx'
        bevx, = get_lapack_funcs((internal_name,), (a1,))
        w, v, m, ifail, info = bevx(
            a1, vl, vu, il, iu, compute_v=not eigvals_only, mmax=max_ev,
            range=select, lower=lower, overwrite_ab=overwrite_a_band,
            abstol=abstol)
        # crop off w and v
        w = w[:m]
        if not eigvals_only:
            v = v[:, :m]
    _check_info(info, internal_name)

    if eigvals_only:
        return w
    return w, v


def eigvals(a, b=None, overwrite_a=False, check_finite=True,
            homogeneous_eigvals=False):
    """
    Compute eigenvalues from an ordinary or generalized eigenvalue problem.

    Find eigenvalues of a general matrix::

        a   vr[:,i] = w[i]        b   vr[:,i]

    Parameters
    ----------
    a : (M, M) array_like
        A complex or real matrix whose eigenvalues and eigenvectors
        will be computed.
    b : (M, M) array_like, optional
        Right-hand side matrix in a generalized eigenvalue problem.
        If omitted, identity matrix is assumed.
    overwrite_a : bool, optional
        Whether to overwrite data in a (may improve performance)
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities
        or NaNs.
    homogeneous_eigvals : bool, optional
        If True, return the eigenvalues in homogeneous coordinates.
        In this case ``w`` is a (2, M) array so that::

            w[1,i] a vr[:,i] = w[0,i] b vr[:,i]

        Default is False.

    Returns
    -------
    w : (M,) or (2, M) double or complex ndarray
        The eigenvalues, each repeated according to its multiplicity
        but not in any specific order. The shape is (M,) unless
        ``homogeneous_eigvals=True``.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge

    See Also
    --------
    eig : eigenvalues and right eigenvectors of general arrays.
    eigvalsh : eigenvalues of symmetric or Hermitian arrays
    eigvals_banded : eigenvalues for symmetric/Hermitian band matrices
    eigvalsh_tridiagonal : eigenvalues of symmetric/Hermitian tridiagonal
        matrices

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.array([[0., -1.], [1., 0.]])
    >>> linalg.eigvals(a)
    array([0.+1.j, 0.-1.j])

    >>> b = np.array([[0., 1.], [1., 1.]])
    >>> linalg.eigvals(a, b)
    array([ 1.+0.j, -1.+0.j])

    >>> a = np.array([[3., 0., 0.], [0., 8., 0.], [0., 0., 7.]])
    >>> linalg.eigvals(a, homogeneous_eigvals=True)
    array([[3.+0.j, 8.+0.j, 7.+0.j],
           [1.+0.j, 1.+0.j, 1.+0.j]])

    """
    return eig(a, b=b, left=0, right=0, overwrite_a=overwrite_a,
               check_finite=check_finite,
               homogeneous_eigvals=homogeneous_eigvals)


def eigvalsh(a, b=None, lower=True, overwrite_a=False,
             overwrite_b=False, turbo=True, eigvals=None, type=1,
             check_finite=True):
    """
    Solve an ordinary or generalized eigenvalue problem for a complex
    Hermitian or real symmetric matrix.

    Find eigenvalues w of matrix a, where b is positive definite::

                      a v[:,i] = w[i] b v[:,i]
        v[i,:].conj() a v[:,i] = w[i]
        v[i,:].conj() b v[:,i] = 1


    Parameters
    ----------
    a : (M, M) array_like
        A complex Hermitian or real symmetric matrix whose eigenvalues and
        eigenvectors will be computed.
    b : (M, M) array_like, optional
        A complex Hermitian or real symmetric definite positive matrix in.
        If omitted, identity matrix is assumed.
    lower : bool, optional
        Whether the pertinent array data is taken from the lower or upper
        triangle of `a`. (Default: lower)
    turbo : bool, optional
        Use divide and conquer algorithm (faster but expensive in memory,
        only for generalized eigenvalue problem and if eigvals=None)
    eigvals : tuple (lo, hi), optional
        Indexes of the smallest and largest (in ascending order) eigenvalues
        and corresponding eigenvectors to be returned: 0 <= lo < hi <= M-1.
        If omitted, all eigenvalues and eigenvectors are returned.
    type : int, optional
        Specifies the problem type to be solved:

           type = 1: a   v[:,i] = w[i] b v[:,i]

           type = 2: a b v[:,i] = w[i]   v[:,i]

           type = 3: b a v[:,i] = w[i]   v[:,i]
    overwrite_a : bool, optional
        Whether to overwrite data in `a` (may improve performance)
    overwrite_b : bool, optional
        Whether to overwrite data in `b` (may improve performance)
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    w : (N,) float ndarray
        The N (1<=N<=M) selected eigenvalues, in ascending order, each
        repeated according to its multiplicity.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge,
        an error occurred, or b matrix is not definite positive. Note that
        if input matrices are not symmetric or hermitian, no error is reported
        but results will be wrong.

    See Also
    --------
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eigvals : eigenvalues of general arrays
    eigvals_banded : eigenvalues for symmetric/Hermitian band matrices
    eigvalsh_tridiagonal : eigenvalues of symmetric/Hermitian tridiagonal
        matrices

    Notes
    -----
    This function does not check the input array for being hermitian/symmetric
    in order to allow for representing arrays with only their upper/lower
    triangular parts.

    Examples
    --------
    >>> from scipy.linalg import eigvalsh
    >>> A = np.array([[6, 3, 1, 5], [3, 0, 5, 1], [1, 5, 6, 2], [5, 1, 2, 2]])
    >>> w = eigvalsh(A)
    >>> w
    array([-3.74637491, -0.76263923,  6.08502336, 12.42399079])

    """
    return eigh(a, b=b, lower=lower, eigvals_only=True,
                overwrite_a=overwrite_a, overwrite_b=overwrite_b,
                turbo=turbo, eigvals=eigvals, type=type,
                check_finite=check_finite)


def eigvals_banded(a_band, lower=False, overwrite_a_band=False,
                   select='a', select_range=None, check_finite=True):
    """
    Solve real symmetric or complex hermitian band matrix eigenvalue problem.

    Find eigenvalues w of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    The matrix a is stored in a_band either in lower diagonal or upper
    diagonal ordered form:

        a_band[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        a_band[    i - j, j] == a[i,j]        (if lower form; i >= j)

    where u is the number of bands above the diagonal.

    Example of a_band (shape of a is (6,6), u=2)::

        upper form:
        *   *   a02 a13 a24 a35
        *   a01 a12 a23 a34 a45
        a00 a11 a22 a33 a44 a55

        lower form:
        a00 a11 a22 a33 a44 a55
        a10 a21 a32 a43 a54 *
        a20 a31 a42 a53 *   *

    Cells marked with * are not used.

    Parameters
    ----------
    a_band : (u+1, M) array_like
        The bands of the M by M matrix a.
    lower : bool, optional
        Is the matrix in the lower form. (Default is upper form)
    overwrite_a_band : bool, optional
        Discard data in a_band (may enhance performance)
    select : {'a', 'v', 'i'}, optional
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max), optional
        Range of selected eigenvalues
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    w : (M,) ndarray
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge.

    See Also
    --------
    eig_banded : eigenvalues and right eigenvectors for symmetric/Hermitian
        band matrices
    eigvalsh_tridiagonal : eigenvalues of symmetric/Hermitian tridiagonal
        matrices
    eigvals : eigenvalues of general arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    Examples
    --------
    >>> from scipy.linalg import eigvals_banded
    >>> A = np.array([[1, 5, 2, 0], [5, 2, 5, 2], [2, 5, 3, 5], [0, 2, 5, 4]])
    >>> Ab = np.array([[1, 2, 3, 4], [5, 5, 5, 0], [2, 2, 0, 0]])
    >>> w = eigvals_banded(Ab, lower=True)
    >>> w
    array([-4.26200532, -2.22987175,  3.95222349, 12.53965359])
    """
    return eig_banded(a_band, lower=lower, eigvals_only=1,
                      overwrite_a_band=overwrite_a_band, select=select,
                      select_range=select_range, check_finite=check_finite)


def eigvalsh_tridiagonal(d, e, select='a', select_range=None,
                         check_finite=True, tol=0., lapack_driver='auto'):
    """
    Solve eigenvalue problem for a real symmetric tridiagonal matrix.

    Find eigenvalues `w` of ``a``::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    For a real symmetric matrix ``a`` with diagonal elements `d` and
    off-diagonal elements `e`.

    Parameters
    ----------
    d : ndarray, shape (ndim,)
        The diagonal elements of the array.
    e : ndarray, shape (ndim-1,)
        The off-diagonal elements of the array.
    select : {'a', 'v', 'i'}, optional
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max), optional
        Range of selected eigenvalues
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    tol : float
        The absolute tolerance to which each eigenvalue is required
        (only used when ``lapack_driver='stebz'``).
        An eigenvalue (or cluster) is considered to have converged if it
        lies in an interval of this width. If <= 0. (default),
        the value ``eps*|a|`` is used where eps is the machine precision,
        and ``|a|`` is the 1-norm of the matrix ``a``.
    lapack_driver : str
        LAPACK function to use, can be 'auto', 'stemr', 'stebz',  'sterf',
        or 'stev'. When 'auto' (default), it will use 'stemr' if ``select='a'``
        and 'stebz' otherwise. 'sterf' and 'stev' can only be used when
        ``select='a'``.

    Returns
    -------
    w : (M,) ndarray
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge.

    See Also
    --------
    eigh_tridiagonal : eigenvalues and right eiegenvectors for
        symmetric/Hermitian tridiagonal matrices

    Examples
    --------
    >>> from scipy.linalg import eigvalsh_tridiagonal, eigvalsh
    >>> d = 3*np.ones(4)
    >>> e = -1*np.ones(3)
    >>> w = eigvalsh_tridiagonal(d, e)
    >>> A = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    >>> w2 = eigvalsh(A)  # Verify with other eigenvalue routines
    >>> np.allclose(w - w2, np.zeros(4))
    True
    """
    return eigh_tridiagonal(
        d, e, eigvals_only=True, select=select, select_range=select_range,
        check_finite=check_finite, tol=tol, lapack_driver=lapack_driver)


def eigh_tridiagonal(d, e, eigvals_only=False, select='a', select_range=None,
                     check_finite=True, tol=0., lapack_driver='auto'):
    """
    Solve eigenvalue problem for a real symmetric tridiagonal matrix.

    Find eigenvalues `w` and optionally right eigenvectors `v` of ``a``::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    For a real symmetric matrix ``a`` with diagonal elements `d` and
    off-diagonal elements `e`.

    Parameters
    ----------
    d : ndarray, shape (ndim,)
        The diagonal elements of the array.
    e : ndarray, shape (ndim-1,)
        The off-diagonal elements of the array.
    select : {'a', 'v', 'i'}, optional
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max), optional
        Range of selected eigenvalues
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    tol : float
        The absolute tolerance to which each eigenvalue is required
        (only used when 'stebz' is the `lapack_driver`).
        An eigenvalue (or cluster) is considered to have converged if it
        lies in an interval of this width. If <= 0. (default),
        the value ``eps*|a|`` is used where eps is the machine precision,
        and ``|a|`` is the 1-norm of the matrix ``a``.
    lapack_driver : str
        LAPACK function to use, can be 'auto', 'stemr', 'stebz', 'sterf',
        or 'stev'. When 'auto' (default), it will use 'stemr' if ``select='a'``
        and 'stebz' otherwise. When 'stebz' is used to find the eigenvalues and
        ``eigvals_only=False``, then a second LAPACK call (to ``?STEIN``) is
        used to find the corresponding eigenvectors. 'sterf' can only be
        used when ``eigvals_only=True`` and ``select='a'``. 'stev' can only
        be used when ``select='a'``.

    Returns
    -------
    w : (M,) ndarray
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.
    v : (M, M) ndarray
        The normalized eigenvector corresponding to the eigenvalue ``w[i]`` is
        the column ``v[:,i]``.

    Raises
    ------
    LinAlgError
        If eigenvalue computation does not converge.

    See Also
    --------
    eigvalsh_tridiagonal : eigenvalues of symmetric/Hermitian tridiagonal
        matrices
    eig : eigenvalues and right eigenvectors for non-symmetric arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eig_banded : eigenvalues and right eigenvectors for symmetric/Hermitian
        band matrices

    Notes
    -----
    This function makes use of LAPACK ``S/DSTEMR`` routines.

    Examples
    --------
    >>> from scipy.linalg import eigh_tridiagonal
    >>> d = 3*np.ones(4)
    >>> e = -1*np.ones(3)
    >>> w, v = eigh_tridiagonal(d, e)
    >>> A = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
    >>> np.allclose(A @ v - v @ np.diag(w), np.zeros((4, 4)))
    True
    """
    d = _asarray_validated(d, check_finite=check_finite)
    e = _asarray_validated(e, check_finite=check_finite)
    for check in (d, e):
        if check.ndim != 1:
            raise ValueError('expected one-dimensional array')
        if check.dtype.char in 'GFD':  # complex
            raise TypeError('Only real arrays currently supported')
    if d.size != e.size + 1:
        raise ValueError('d (%s) must have one more element than e (%s)'
                         % (d.size, e.size))
    select, vl, vu, il, iu, _ = _check_select(
        select, select_range, 0, d.size)
    if not isinstance(lapack_driver, string_types):
        raise TypeError('lapack_driver must be str')
    drivers = ('auto', 'stemr', 'sterf', 'stebz', 'stev')
    if lapack_driver not in drivers:
        raise ValueError('lapack_driver must be one of %s, got %s'
                         % (drivers, lapack_driver))
    if lapack_driver == 'auto':
        lapack_driver = 'stemr' if select == 0 else 'stebz'
    func, = get_lapack_funcs((lapack_driver,), (d, e))
    compute_v = not eigvals_only
    if lapack_driver == 'sterf':
        if select != 0:
            raise ValueError('sterf can only be used when select == "a"')
        if not eigvals_only:
            raise ValueError('sterf can only be used when eigvals_only is '
                             'True')
        w, info = func(d, e)
        m = len(w)
    elif lapack_driver == 'stev':
        if select != 0:
            raise ValueError('stev can only be used when select == "a"')
        w, v, info = func(d, e, compute_v=compute_v)
        m = len(w)
    elif lapack_driver == 'stebz':
        tol = float(tol)
        internal_name = 'stebz'
        stebz, = get_lapack_funcs((internal_name,), (d, e))
        # If getting eigenvectors, needs to be block-ordered (B) instead of
        # matirx-ordered (E), and we will reorder later
        order = 'E' if eigvals_only else 'B'
        m, w, iblock, isplit, info = stebz(d, e, select, vl, vu, il, iu, tol,
                                           order)
    else:   # 'stemr'
        # ?STEMR annoyingly requires size N instead of N-1
        e_ = empty(e.size+1, e.dtype)
        e_[:-1] = e
        stemr_lwork, = get_lapack_funcs(('stemr_lwork',), (d, e))
        lwork, liwork, info = stemr_lwork(d, e_, select, vl, vu, il, iu,
                                          compute_v=compute_v)
        _check_info(info, 'stemr_lwork')
        m, w, v, info = func(d, e_, select, vl, vu, il, iu,
                             compute_v=compute_v, lwork=lwork, liwork=liwork)
    _check_info(info, lapack_driver + ' (eigh_tridiagonal)')
    w = w[:m]
    if eigvals_only:
        return w
    else:
        # Do we still need to compute the eigenvalues?
        if lapack_driver == 'stebz':
            func, = get_lapack_funcs(('stein',), (d, e))
            v, info = func(d, e, w, iblock, isplit)
            _check_info(info, 'stein (eigh_tridiagonal)',
                        positive='%d eigenvectors failed to converge')
            # Convert block-order to matrix-order
            order = argsort(w)
            w, v = w[order], v[:, order]
        else:
            v = v[:, :m]
        return w, v


def _check_info(info, driver, positive='did not converge (LAPACK info=%d)'):
    """Check info return value."""
    if info < 0:
        raise ValueError('illegal value in argument %d of internal %s'
                         % (-info, driver))
    if info > 0 and positive:
        raise LinAlgError(("%s " + positive) % (driver, info,))


def hessenberg(a, calc_q=False, overwrite_a=False, check_finite=True):
    """
    Compute Hessenberg form of a matrix.

    The Hessenberg decomposition is::

        A = Q H Q^H

    where `Q` is unitary/orthogonal and `H` has only zero elements below
    the first sub-diagonal.

    Parameters
    ----------
    a : (M, M) array_like
        Matrix to bring into Hessenberg form.
    calc_q : bool, optional
        Whether to compute the transformation matrix.  Default is False.
    overwrite_a : bool, optional
        Whether to overwrite `a`; may improve performance.
        Default is False.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    H : (M, M) ndarray
        Hessenberg form of `a`.
    Q : (M, M) ndarray
        Unitary/orthogonal similarity transformation matrix ``A = Q H Q^H``.
        Only returned if ``calc_q=True``.

    Examples
    --------
    >>> from scipy.linalg import hessenberg
    >>> A = np.array([[2, 5, 8, 7], [5, 2, 2, 8], [7, 5, 6, 6], [5, 4, 4, 8]])
    >>> H, Q = hessenberg(A, calc_q=True)
    >>> H
    array([[  2.        , -11.65843866,   1.42005301,   0.25349066],
           [ -9.94987437,  14.53535354,  -5.31022304,   2.43081618],
           [  0.        ,  -1.83299243,   0.38969961,  -0.51527034],
           [  0.        ,   0.        ,  -3.83189513,   1.07494686]])
    >>> np.allclose(Q @ H @ Q.conj().T - A, np.zeros((4, 4)))
    True
    """
    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2 or (a1.shape[0] != a1.shape[1]):
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    # if 2x2 or smaller: already in Hessenberg
    if a1.shape[0] <= 2:
        if calc_q:
            return a1, numpy.eye(a1.shape[0])
        return a1

    gehrd, gebal, gehrd_lwork = get_lapack_funcs(('gehrd', 'gebal',
                                                  'gehrd_lwork'), (a1,))
    ba, lo, hi, pivscale, info = gebal(a1, permute=0, overwrite_a=overwrite_a)
    _check_info(info, 'gebal (hessenberg)', positive=False)
    n = len(a1)

    lwork = _compute_lwork(gehrd_lwork, ba.shape[0], lo=lo, hi=hi)

    hq, tau, info = gehrd(ba, lo=lo, hi=hi, lwork=lwork, overwrite_a=1)
    _check_info(info, 'gehrd (hessenberg)', positive=False)
    h = numpy.triu(hq, -1)
    if not calc_q:
        return h

    # use orghr/unghr to compute q
    orghr, orghr_lwork = get_lapack_funcs(('orghr', 'orghr_lwork'), (a1,))
    lwork = _compute_lwork(orghr_lwork, n, lo=lo, hi=hi)

    q, info = orghr(a=hq, tau=tau, lo=lo, hi=hi, lwork=lwork, overwrite_a=1)
    _check_info(info, 'orghr (hessenberg)', positive=False)
    return h, q


def cdf2rdf(w, v):
    """
    Converts complex eigenvalues ``w`` and eigenvectors ``v`` to real
    eigenvalues in a block diagonal form ``wr`` and the associated real
    eigenvectors ``vr``, such that::

        vr @ wr = X @ vr

    continues to hold, where ``X`` is the original array for which ``w`` and
    ``v`` are the eigenvalues and eigenvectors.

    .. versionadded:: 1.1.0

    Parameters
    ----------
    w : (..., M) array_like
        Complex or real eigenvalues, an array or stack of arrays

        Conjugate pairs must not be interleaved, else the wrong result
        will be produced. So ``[1+1j, 1, 1-1j]`` will give a correct result, but
        ``[1+1j, 2+1j, 1-1j, 2-1j]`` will not.

    v : (..., M, M) array_like
        Complex or real eigenvectors, a square array or stack of square arrays.

    Returns
    -------
    wr : (..., M, M) ndarray
        Real diagonal block form of eigenvalues
    vr : (..., M, M) ndarray
        Real eigenvectors associated with ``wr``

    See Also
    --------
    eig : Eigenvalues and right eigenvectors for non-symmetric arrays
    rsf2csf : Convert real Schur form to complex Schur form

    Notes
    -----
    ``w``, ``v`` must be the eigenstructure for some *real* matrix ``X``.
    For example, obtained by ``w, v = scipy.linalg.eig(X)`` or
    ``w, v = numpy.linalg.eig(X)`` in which case ``X`` can also represent
    stacked arrays.

    .. versionadded:: 1.1.0

    Examples
    --------
    >>> X = np.array([[1, 2, 3], [0, 4, 5], [0, -5, 4]])
    >>> X
    array([[ 1,  2,  3],
           [ 0,  4,  5],
           [ 0, -5,  4]])

    >>> from scipy import linalg
    >>> w, v = linalg.eig(X)
    >>> w
    array([ 1.+0.j,  4.+5.j,  4.-5.j])
    >>> v
    array([[ 1.00000+0.j     , -0.01906-0.40016j, -0.01906+0.40016j],
           [ 0.00000+0.j     ,  0.00000-0.64788j,  0.00000+0.64788j],
           [ 0.00000+0.j     ,  0.64788+0.j     ,  0.64788-0.j     ]])

    >>> wr, vr = linalg.cdf2rdf(w, v)
    >>> wr
    array([[ 1.,  0.,  0.],
           [ 0.,  4.,  5.],
           [ 0., -5.,  4.]])
    >>> vr
    array([[ 1.     ,  0.40016, -0.01906],
           [ 0.     ,  0.64788,  0.     ],
           [ 0.     ,  0.     ,  0.64788]])

    >>> vr @ wr
    array([[ 1.     ,  1.69593,  1.9246 ],
           [ 0.     ,  2.59153,  3.23942],
           [ 0.     , -3.23942,  2.59153]])
    >>> X @ vr
    array([[ 1.     ,  1.69593,  1.9246 ],
           [ 0.     ,  2.59153,  3.23942],
           [ 0.     , -3.23942,  2.59153]])
    """
    w, v = _asarray_validated(w), _asarray_validated(v)

    # check dimensions
    if w.ndim < 1:
        raise ValueError('expected w to be at least one-dimensional')
    if v.ndim < 2:
        raise ValueError('expected v to be at least two-dimensional')
    if v.ndim != w.ndim + 1:
        raise ValueError('expected eigenvectors array to have exactly one '
                         'dimension more than eigenvalues array')

    # check shapes
    n = w.shape[-1]
    M = w.shape[:-1]
    if v.shape[-2] != v.shape[-1]:
        raise ValueError('expected v to be a square matrix or stacked square '
                         'matrices: v.shape[-2] = v.shape[-1]')
    if v.shape[-1] != n:
        raise ValueError('expected the same number of eigenvalues as '
                         'eigenvectors')

    # get indices for each first pair of complex eigenvalues
    complex_mask = iscomplex(w)
    n_complex = complex_mask.sum(axis=-1)

    # check if all complex eigenvalues have conjugate pairs
    if not (n_complex % 2 == 0).all():
        raise ValueError('expected complex-conjugate pairs of eigenvalues')

    # find complex indices
    idx = nonzero(complex_mask)
    idx_stack = idx[:-1]
    idx_elem = idx[-1]

    # filter them to conjugate indices, assuming pairs are not interleaved
    j = idx_elem[0::2]
    k = idx_elem[1::2]
    stack_ind = ()
    for i in idx_stack:
        # should never happen, assuming nonzero orders by the last axis
        assert (i[0::2] == i[1::2]).all(), "Conjugate pair spanned different arrays!"
        stack_ind += (i[0::2],)

    # all eigenvalues to diagonal form
    wr = zeros(M + (n, n), dtype=w.real.dtype)
    di = range(n)
    wr[..., di, di] = w.real

    # complex eigenvalues to real block diagonal form
    wr[stack_ind + (j, k)] = w[stack_ind + (j,)].imag
    wr[stack_ind + (k, j)] = w[stack_ind + (k,)].imag

    # compute real eigenvectors associated with real block diagonal eigenvalues
    u = zeros(M + (n, n), dtype=numpy.cdouble)
    u[..., di, di] = 1.0
    u[stack_ind + (j, j)] = 0.5j
    u[stack_ind + (j, k)] = 0.5
    u[stack_ind + (k, j)] = -0.5j
    u[stack_ind + (k, k)] = 0.5

    # multipy matrices v and u (equivalent to v @ u)
    vr = einsum('...ij,...jk->...ik', v, u).real

    return wr, vr
