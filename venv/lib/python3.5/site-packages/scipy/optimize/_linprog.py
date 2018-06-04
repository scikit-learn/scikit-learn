"""
A top-level linear programming interface. Currently this interface only
solves linear programming problems via the Simplex Method.

.. versionadded:: 0.15.0

Functions
---------
.. autosummary::
   :toctree: generated/

    linprog
    linprog_verbose_callback
    linprog_terse_callback

"""

from __future__ import division, print_function, absolute_import

import numpy as np
from .optimize import OptimizeResult, _check_unknown_options
from ._linprog_ip import _linprog_ip

__all__ = ['linprog', 'linprog_verbose_callback', 'linprog_terse_callback']

__docformat__ = "restructuredtext en"


def linprog_verbose_callback(xk, **kwargs):
    """
    A sample callback function demonstrating the linprog callback interface.
    This callback produces detailed output to sys.stdout before each iteration
    and after the final iteration of the simplex algorithm.

    Parameters
    ----------
    xk : array_like
        The current solution vector.
    **kwargs : dict
        A dictionary containing the following parameters:

        tableau : array_like
            The current tableau of the simplex algorithm.
            Its structure is defined in _solve_simplex.
        phase : int
            The current Phase of the simplex algorithm (1 or 2)
        nit : int
            The current iteration number.
        pivot : tuple(int, int)
            The index of the tableau selected as the next pivot,
            or nan if no pivot exists
        basis : array(int)
            A list of the current basic variables.
            Each element contains the name of a basic variable and its value.
        complete : bool
            True if the simplex algorithm has completed
            (and this is the final call to callback), otherwise False.
    """
    tableau = kwargs["tableau"]
    nit = kwargs["nit"]
    pivrow, pivcol = kwargs["pivot"]
    phase = kwargs["phase"]
    basis = kwargs["basis"]
    complete = kwargs["complete"]

    saved_printoptions = np.get_printoptions()
    np.set_printoptions(linewidth=500,
                        formatter={'float': lambda x: "{0: 12.4f}".format(x)})
    if complete:
        print("--------- Iteration Complete - Phase {0:d} -------\n".format(phase))
        print("Tableau:")
    elif nit == 0:
        print("--------- Initial Tableau - Phase {0:d} ----------\n".format(phase))

    else:
        print("--------- Iteration {0:d}  - Phase {1:d} --------\n".format(nit, phase))
        print("Tableau:")

    if nit >= 0:
        print("" + str(tableau) + "\n")
        if not complete:
            print("Pivot Element: T[{0:.0f}, {1:.0f}]\n".format(pivrow, pivcol))
        print("Basic Variables:", basis)
        print()
        print("Current Solution:")
        print("x = ", xk)
        print()
        print("Current Objective Value:")
        print("f = ", -tableau[-1, -1])
        print()
    np.set_printoptions(**saved_printoptions)


def linprog_terse_callback(xk, **kwargs):
    """
    A sample callback function demonstrating the linprog callback interface.
    This callback produces brief output to sys.stdout before each iteration
    and after the final iteration of the simplex algorithm.

    Parameters
    ----------
    xk : array_like
        The current solution vector.
    **kwargs : dict
        A dictionary containing the following parameters:

        tableau : array_like
            The current tableau of the simplex algorithm.
            Its structure is defined in _solve_simplex.
        vars : tuple(str, ...)
            Column headers for each column in tableau.
            "x[i]" for actual variables, "s[i]" for slack surplus variables,
            "a[i]" for artificial variables, and "RHS" for the constraint
            RHS vector.
        phase : int
            The current Phase of the simplex algorithm (1 or 2)
        nit : int
            The current iteration number.
        pivot : tuple(int, int)
            The index of the tableau selected as the next pivot,
            or nan if no pivot exists
        basics : list[tuple(int, float)]
            A list of the current basic variables.
            Each element contains the index of a basic variable and
            its value.
        complete : bool
            True if the simplex algorithm has completed
            (and this is the final call to callback), otherwise False.
    """
    nit = kwargs["nit"]

    if nit == 0:
        print("Iter:   X:")
    print("{0: <5d}   ".format(nit), end="")
    print(xk)


def _pivot_col(T, tol=1.0E-12, bland=False):
    """
    Given a linear programming simplex tableau, determine the column
    of the variable to enter the basis.

    Parameters
    ----------
    T : 2D ndarray
        The simplex tableau.
    tol : float
        Elements in the objective row larger than -tol will not be considered
        for pivoting.  Nominally this value is zero, but numerical issues
        cause a tolerance about zero to be necessary.
    bland : bool
        If True, use Bland's rule for selection of the column (select the
        first column with a negative coefficient in the objective row,
        regardless of magnitude).

    Returns
    -------
    status: bool
        True if a suitable pivot column was found, otherwise False.
        A return of False indicates that the linear programming simplex
        algorithm is complete.
    col: int
        The index of the column of the pivot element.
        If status is False, col will be returned as nan.
    """
    ma = np.ma.masked_where(T[-1, :-1] >= -tol, T[-1, :-1], copy=False)
    if ma.count() == 0:
        return False, np.nan
    if bland:
        return True, np.where(ma.mask == False)[0][0]
    return True, np.ma.where(ma == ma.min())[0][0]


def _pivot_row(T, basis, pivcol, phase, tol=1.0E-12, bland=False):
    """
    Given a linear programming simplex tableau, determine the row for the
    pivot operation.

    Parameters
    ----------
    T : 2D ndarray
        The simplex tableau.
    basis : array
        A list of the current basic variables.
    pivcol : int
        The index of the pivot column.
    phase : int
        The phase of the simplex algorithm (1 or 2).
    tol : float
        Elements in the pivot column smaller than tol will not be considered
        for pivoting.  Nominally this value is zero, but numerical issues
        cause a tolerance about zero to be necessary.
    bland : bool
        If True, use Bland's rule for selection of the row (if more than one
        row can be used, choose the one with the lowest variable index).

    Returns
    -------
    status: bool
        True if a suitable pivot row was found, otherwise False.  A return
        of False indicates that the linear programming problem is unbounded.
    row: int
        The index of the row of the pivot element.  If status is False, row
        will be returned as nan.
    """
    if phase == 1:
        k = 2
    else:
        k = 1
    ma = np.ma.masked_where(T[:-k, pivcol] <= tol, T[:-k, pivcol], copy=False)
    if ma.count() == 0:
        return False, np.nan
    mb = np.ma.masked_where(T[:-k, pivcol] <= tol, T[:-k, -1], copy=False)
    q = mb / ma
    min_rows = np.ma.where(q == q.min())[0]
    if bland:
        return True, min_rows[np.argmin(np.take(basis, min_rows))]
    return True, min_rows[0]


def _solve_simplex(T, n, basis, maxiter=1000, phase=2, callback=None,
                   tol=1.0E-12, nit0=0, bland=False):
    """
    Solve a linear programming problem in "standard maximization form" using
    the Simplex Method.

    Minimize :math:`f = c^T x`

    subject to

    .. math::

        Ax = b
        x_i >= 0
        b_j >= 0

    Parameters
    ----------
    T : array_like
        A 2-D array representing the simplex T corresponding to the
        maximization problem.  It should have the form:

        [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
         [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
         .
         .
         .
         [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
         [c[0],   c[1], ...,   c[n_total],    0]]

        for a Phase 2 problem, or the form:

        [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
         [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
         .
         .
         .
         [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
         [c[0],   c[1], ...,   c[n_total],   0],
         [c'[0],  c'[1], ...,  c'[n_total],  0]]

         for a Phase 1 problem (a Problem in which a basic feasible solution is
         sought prior to maximizing the actual objective.  T is modified in
         place by _solve_simplex.
    n : int
        The number of true variables in the problem.
    basis : array
        An array of the indices of the basic variables, such that basis[i]
        contains the column corresponding to the basic variable for row i.
        Basis is modified in place by _solve_simplex
    maxiter : int
        The maximum number of iterations to perform before aborting the
        optimization.
    phase : int
        The phase of the optimization being executed.  In phase 1 a basic
        feasible solution is sought and the T has an additional row
        representing an alternate objective function.
    callback : callable, optional
        If a callback function is provided, it will be called within each
        iteration of the simplex algorithm. The callback must have the
        signature `callback(xk, **kwargs)` where xk is the current solution
        vector and kwargs is a dictionary containing the following::
        "T" : The current Simplex algorithm T
        "nit" : The current iteration.
        "pivot" : The pivot (row, column) used for the next iteration.
        "phase" : Whether the algorithm is in Phase 1 or Phase 2.
        "basis" : The indices of the columns of the basic variables.
    tol : float
        The tolerance which determines when a solution is "close enough" to
        zero in Phase 1 to be considered a basic feasible solution or close
        enough to positive to serve as an optimal solution.
    nit0 : int
        The initial iteration number used to keep an accurate iteration total
        in a two-phase problem.
    bland : bool
        If True, choose pivots using Bland's rule [3].  In problems which
        fail to converge due to cycling, using Bland's rule can provide
        convergence at the expense of a less optimal path about the simplex.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. Possible
        values for the ``status`` attribute are:
         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded

        See `OptimizeResult` for a description of other attributes.
    """
    nit = nit0
    complete = False

    if phase == 1:
        m = T.shape[0]-2
    elif phase == 2:
        m = T.shape[0]-1
    else:
        raise ValueError("Argument 'phase' to _solve_simplex must be 1 or 2")

    if phase == 2:
        # Check if any artificial variables are still in the basis.
        # If yes, check if any coefficients from this row and a column
        # corresponding to one of the non-artificial variable is non-zero.
        # If found, pivot at this term. If not, start phase 2.
        # Do this for all artificial variables in the basis.
        # Ref: "An Introduction to Linear Programming and Game Theory"
        # by Paul R. Thie, Gerard E. Keough, 3rd Ed,
        # Chapter 3.7 Redundant Systems (pag 102)
        for pivrow in [row for row in range(basis.size)
                       if basis[row] > T.shape[1] - 2]:
            non_zero_row = [col for col in range(T.shape[1] - 1)
                            if T[pivrow, col] != 0]
            if len(non_zero_row) > 0:
                pivcol = non_zero_row[0]
                # variable represented by pivcol enters
                # variable in basis[pivrow] leaves
                basis[pivrow] = pivcol
                pivval = T[pivrow][pivcol]
                T[pivrow, :] = T[pivrow, :] / pivval
                for irow in range(T.shape[0]):
                    if irow != pivrow:
                        T[irow, :] = T[irow, :] - T[pivrow, :]*T[irow, pivcol]
                nit += 1

    if len(basis[:m]) == 0:
        solution = np.zeros(T.shape[1] - 1, dtype=np.float64)
    else:
        solution = np.zeros(max(T.shape[1] - 1, max(basis[:m]) + 1),
                            dtype=np.float64)

    while not complete:
        # Find the pivot column
        pivcol_found, pivcol = _pivot_col(T, tol, bland)
        if not pivcol_found:
            pivcol = np.nan
            pivrow = np.nan
            status = 0
            complete = True
        else:
            # Find the pivot row
            pivrow_found, pivrow = _pivot_row(T, basis, pivcol, phase, tol, bland)
            if not pivrow_found:
                status = 3
                complete = True

        if callback is not None:
            solution[:] = 0
            solution[basis[:m]] = T[:m, -1]
            callback(solution[:n], **{"tableau": T,
                                      "phase": phase,
                                      "nit": nit,
                                      "pivot": (pivrow, pivcol),
                                      "basis": basis,
                                      "complete": complete and phase == 2})

        if not complete:
            if nit >= maxiter:
                # Iteration limit exceeded
                status = 1
                complete = True
            else:
                # variable represented by pivcol enters
                # variable in basis[pivrow] leaves
                basis[pivrow] = pivcol
                pivval = T[pivrow][pivcol]
                T[pivrow, :] = T[pivrow, :] / pivval
                for irow in range(T.shape[0]):
                    if irow != pivrow:
                        T[irow, :] = T[irow, :] - T[pivrow, :]*T[irow, pivcol]
                nit += 1

    return nit, status


def _linprog_simplex(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None,
                     bounds=None, maxiter=1000, disp=False, callback=None,
                     tol=1.0E-12, bland=False, **unknown_options):
    """
    Solve the following linear programming problem via a two-phase
    simplex algorithm.::

        minimize:     c^T * x

        subject to:   A_ub * x <= b_ub
                      A_eq * x == b_eq

    Parameters
    ----------
    c : array_like
        Coefficients of the linear objective function to be minimized.
    A_ub : array_like
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``.
    b_ub : array_like
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : array_like
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``.
    b_eq : array_like
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : array_like
        The bounds for each independent variable in the solution, which can
        take one of three forms::

        None : The default bounds, all variables are non-negative.
        (lb, ub) : If a 2-element sequence is provided, the same
                  lower bound (lb) and upper bound (ub) will be applied
                  to all variables.
        [(lb_0, ub_0), (lb_1, ub_1), ...] : If an n x 2 sequence is provided,
                  each variable x_i will be bounded by lb[i] and ub[i].
        Infinite bounds are specified using -np.inf (negative)
        or np.inf (positive).

    callback : callable
        If a callback function is provide, it will be called within each
        iteration of the simplex algorithm. The callback must have the
        signature ``callback(xk, **kwargs)`` where ``xk`` is the current s
        olution vector and kwargs is a dictionary containing the following::

        "tableau" : The current Simplex algorithm tableau
        "nit" : The current iteration.
        "pivot" : The pivot (row, column) used for the next iteration.
        "phase" : Whether the algorithm is in Phase 1 or Phase 2.
        "bv" : A structured array containing a string representation of each
               basic variable and its current value.

    Options
    -------
    maxiter : int
       The maximum number of iterations to perform.
    disp : bool
        If True, print exit status message to sys.stdout
    tol : float
        The tolerance which determines when a solution is "close enough" to
        zero in Phase 1 to be considered a basic feasible solution or close
        enough to positive to serve as an optimal solution.
    bland : bool
        If True, use Bland's anti-cycling rule [3] to choose pivots to
        prevent cycling.  If False, choose pivots which should lead to a
        converged solution more quickly.  The latter method is subject to
        cycling (non-convergence) in rare instances.

    Returns
    -------
    A `scipy.optimize.OptimizeResult` consisting of the following fields:

        x : ndarray
            The independent variable vector which optimizes the linear
            programming problem.
        fun : float
            Value of the objective function.
        slack : ndarray
            The values of the slack variables.  Each slack variable corresponds
            to an inequality constraint.  If the slack is zero, then the
            corresponding constraint is active.
        success : bool
            Returns True if the algorithm succeeded in finding an optimal
            solution.
        status : int
            An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded

        nit : int
            The number of iterations performed.
        message : str
            A string descriptor of the exit status of the optimization.

    Examples
    --------
    Consider the following problem:

    Minimize: f = -1*x[0] + 4*x[1]

    Subject to: -3*x[0] + 1*x[1] <= 6
                 1*x[0] + 2*x[1] <= 4
                            x[1] >= -3

    where:  -inf <= x[0] <= inf

    This problem deviates from the standard linear programming problem.  In
    standard form, linear programming problems assume the variables x are
    non-negative.  Since the variables don't have standard bounds where
    0 <= x <= inf, the bounds of the variables must be explicitly set.

    There are two upper-bound constraints, which can be expressed as

    dot(A_ub, x) <= b_ub

    The input for this problem is as follows:

    >>> from scipy.optimize import linprog
    >>> c = [-1, 4]
    >>> A = [[-3, 1], [1, 2]]
    >>> b = [6, 4]
    >>> x0_bnds = (None, None)
    >>> x1_bnds = (-3, None)
    >>> res = linprog(c, A, b, bounds=(x0_bnds, x1_bnds))
    >>> print(res)
         fun: -22.0
     message: 'Optimization terminated successfully.'
         nit: 1
       slack: array([ 39.,   0.])
      status: 0
     success: True
           x: array([ 10.,  -3.])

    References
    ----------
    .. [1] Dantzig, George B., Linear programming and extensions. Rand
           Corporation Research Study Princeton Univ. Press, Princeton, NJ,
           1963
    .. [2] Hillier, S.H. and Lieberman, G.J. (1995), "Introduction to
           Mathematical Programming", McGraw-Hill, Chapter 4.
    .. [3] Bland, Robert G. New finite pivoting rules for the simplex method.
           Mathematics of Operations Research (2), 1977: pp. 103-107.
    """
    _check_unknown_options(unknown_options)

    status = 0
    messages = {0: "Optimization terminated successfully.",
                1: "Iteration limit reached.",
                2: "Optimization failed. Unable to find a feasible"
                   " starting point.",
                3: "Optimization failed. The problem appears to be unbounded.",
                4: "Optimization failed. Singular matrix encountered."}
    have_floor_variable = False

    cc = np.asarray(c)

    # The initial value of the objective function element in the tableau
    f0 = 0

    # The number of variables as given by c
    n = len(c)

    # Convert the input arguments to arrays (sized to zero if not provided)
    Aeq = np.asarray(A_eq) if A_eq is not None else np.empty([0, len(cc)])
    Aub = np.asarray(A_ub) if A_ub is not None else np.empty([0, len(cc)])
    beq = np.ravel(np.asarray(b_eq)) if b_eq is not None else np.empty([0])
    bub = np.ravel(np.asarray(b_ub)) if b_ub is not None else np.empty([0])

    # Analyze the bounds and determine what modifications to be made to
    # the constraints in order to accommodate them.
    L = np.zeros(n, dtype=np.float64)
    U = np.ones(n, dtype=np.float64)*np.inf
    if bounds is None or len(bounds) == 0:
        pass
    elif len(bounds) == 2 and not hasattr(bounds[0], '__len__'):
        # All bounds are the same
        a = bounds[0] if bounds[0] is not None else -np.inf
        b = bounds[1] if bounds[1] is not None else np.inf
        L = np.asarray(n*[a], dtype=np.float64)
        U = np.asarray(n*[b], dtype=np.float64)
    else:
        if len(bounds) != n:
            status = -1
            message = ("Invalid input for linprog with method = 'simplex'.  "
                       "Length of bounds is inconsistent with the length of c")
        else:
            try:
                for i in range(n):
                    if len(bounds[i]) != 2:
                        raise IndexError()
                    L[i] = bounds[i][0] if bounds[i][0] is not None else -np.inf
                    U[i] = bounds[i][1] if bounds[i][1] is not None else np.inf
            except IndexError:
                status = -1
                message = ("Invalid input for linprog with "
                           "method = 'simplex'.  bounds must be a n x 2 "
                           "sequence/array where n = len(c).")

    if np.any(L == -np.inf):
        # If any lower-bound constraint is a free variable
        # add the first column variable as the "floor" variable which
        # accommodates the most negative variable in the problem.
        n = n + 1
        L = np.concatenate([np.array([0]), L])
        U = np.concatenate([np.array([np.inf]), U])
        cc = np.concatenate([np.array([0]), cc])
        Aeq = np.hstack([np.zeros([Aeq.shape[0], 1]), Aeq])
        Aub = np.hstack([np.zeros([Aub.shape[0], 1]), Aub])
        have_floor_variable = True

    # Now before we deal with any variables with lower bounds < 0,
    # deal with finite bounds which can be simply added as new constraints.
    # Also validate bounds inputs here.
    for i in range(n):
        if(L[i] > U[i]):
            status = -1
            message = ("Invalid input for linprog with method = 'simplex'.  "
                       "Lower bound %d is greater than upper bound%d" % (i, i))

        if np.isinf(L[i]) and L[i] > 0:
            status = -1
            message = ("Invalid input for linprog with method = 'simplex'.  "
                       "Lower bound may not be +infinity")

        if np.isinf(U[i]) and U[i] < 0:
            status = -1
            message = ("Invalid input for linprog with method = 'simplex'.  "
                       "Upper bound may not be -infinity")

        if np.isfinite(L[i]) and L[i] > 0:
            # Add a new lower-bound (negative upper-bound) constraint
            Aub = np.vstack([Aub, np.zeros(n)])
            Aub[-1, i] = -1
            bub = np.concatenate([bub, np.array([-L[i]])])
            L[i] = 0

        if np.isfinite(U[i]):
            # Add a new upper-bound constraint
            Aub = np.vstack([Aub, np.zeros(n)])
            Aub[-1, i] = 1
            bub = np.concatenate([bub, np.array([U[i]])])
            U[i] = np.inf

    # Now find negative lower bounds (finite or infinite) which require a
    # change of variables or free variables and handle them appropriately
    for i in range(0, n):
        if L[i] < 0:
            if np.isfinite(L[i]) and L[i] < 0:
                # Add a change of variables for x[i]
                # For each row in the constraint matrices, we take the
                # coefficient from column i in A,
                # and subtract the product of that and L[i] to the RHS b
                beq = beq - Aeq[:, i] * L[i]
                bub = bub - Aub[:, i] * L[i]
                # We now have a nonzero initial value for the objective
                # function as well.
                f0 = f0 - cc[i] * L[i]
            else:
                # This is an unrestricted variable, let x[i] = u[i] - v[0]
                # where v is the first column in all matrices.
                Aeq[:, 0] = Aeq[:, 0] - Aeq[:, i]
                Aub[:, 0] = Aub[:, 0] - Aub[:, i]
                cc[0] = cc[0] - cc[i]

        if np.isinf(U[i]):
            if U[i] < 0:
                status = -1
                message = ("Invalid input for linprog with "
                           "method = 'simplex'.  Upper bound may not be -inf.")

    # The number of upper bound constraints (rows in A_ub and elements in b_ub)
    mub = len(bub)

    # The number of equality constraints (rows in A_eq and elements in b_eq)
    meq = len(beq)

    # The total number of constraints
    m = mub+meq

    # The number of slack variables (one for each upper-bound constraints)
    n_slack = mub

    # The number of artificial variables (one for each lower-bound and equality
    # constraint)
    n_artificial = meq + np.count_nonzero(bub < 0)

    try:
        Aub_rows, Aub_cols = Aub.shape
    except ValueError:
        raise ValueError("Invalid input.  A_ub must be two-dimensional")

    try:
        Aeq_rows, Aeq_cols = Aeq.shape
    except ValueError:
        raise ValueError("Invalid input.  A_eq must be two-dimensional")

    if Aeq_rows != meq:
        status = -1
        message = ("Invalid input for linprog with method = 'simplex'.  "
                   "The number of rows in A_eq must be equal "
                   "to the number of values in b_eq")

    if Aub_rows != mub:
        status = -1
        message = ("Invalid input for linprog with method = 'simplex'.  "
                   "The number of rows in A_ub must be equal "
                   "to the number of values in b_ub")

    if Aeq_cols > 0 and Aeq_cols != n:
        status = -1
        message = ("Invalid input for linprog with method = 'simplex'.  "
                   "Number of columns in A_eq must be equal "
                   "to the size of c")

    if Aub_cols > 0 and Aub_cols != n:
        status = -1
        message = ("Invalid input for linprog with method = 'simplex'.  "
                   "Number of columns in A_ub must be equal to the size of c")

    if status != 0:
        # Invalid inputs provided
        raise ValueError(message)

    # Create the tableau
    T = np.zeros([m+2, n+n_slack+n_artificial+1])

    # Insert objective into tableau
    T[-2, :n] = cc
    T[-2, -1] = f0

    b = T[:-2, -1]

    if meq > 0:
        # Add Aeq to the tableau
        T[:meq, :n] = Aeq
        # Add beq to the tableau
        b[:meq] = beq
    if mub > 0:
        # Add Aub to the tableau
        T[meq:meq+mub, :n] = Aub
        # At bub to the tableau
        b[meq:meq+mub] = bub
        # Add the slack variables to the tableau
        np.fill_diagonal(T[meq:m, n:n+n_slack], 1)

    # Further set up the tableau.
    # If a row corresponds to an equality constraint or a negative b (a lower
    # bound constraint), then an artificial variable is added for that row.
    # Also, if b is negative, first flip the signs in that constraint.
    slcount = 0
    avcount = 0
    basis = np.zeros(m, dtype=int)
    r_artificial = np.zeros(n_artificial, dtype=int)
    for i in range(m):
        if i < meq or b[i] < 0:
            # basic variable i is in column n+n_slack+avcount
            basis[i] = n+n_slack+avcount
            r_artificial[avcount] = i
            avcount += 1
            if b[i] < 0:
                b[i] *= -1
                T[i, :-1] *= -1
            T[i, basis[i]] = 1
            T[-1, basis[i]] = 1
        else:
            # basic variable i is in column n+slcount
            basis[i] = n+slcount
            slcount += 1

    # Make the artificial variables basic feasible variables by subtracting
    # each row with an artificial variable from the Phase 1 objective
    for r in r_artificial:
        T[-1, :] = T[-1, :] - T[r, :]

    nit1, status = _solve_simplex(T, n, basis, phase=1, callback=callback,
                                  maxiter=maxiter, tol=tol, bland=bland)

    # if pseudo objective is zero, remove the last row from the tableau and
    # proceed to phase 2
    if abs(T[-1, -1]) < tol:
        # Remove the pseudo-objective row from the tableau
        T = T[:-1, :]
        # Remove the artificial variable columns from the tableau
        T = np.delete(T, np.s_[n+n_slack:n+n_slack+n_artificial], 1)
    else:
        # Failure to find a feasible starting point
        status = 2

    if status != 0:
        message = messages[status]
        if disp:
            print(message)
        return OptimizeResult(x=np.nan, fun=-T[-1, -1], nit=nit1,
                              status=status, message=message, success=False)

    # Phase 2
    nit2, status = _solve_simplex(T, n, basis, maxiter=maxiter-nit1, phase=2,
                                  callback=callback, tol=tol, nit0=nit1,
                                  bland=bland)

    solution = np.zeros(n+n_slack+n_artificial)
    solution[basis[:m]] = T[:m, -1]
    x = solution[:n]
    slack = solution[n:n+n_slack]

    # For those variables with finite negative lower bounds,
    # reverse the change of variables
    masked_L = np.ma.array(L, mask=np.isinf(L), fill_value=0.0).filled()
    x = x + masked_L

    # For those variables with infinite negative lower bounds,
    # take x[i] as the difference between x[i] and the floor variable.
    if have_floor_variable:
        for i in range(1, n):
            if np.isinf(L[i]):
                x[i] -= x[0]
        x = x[1:]

    # Optimization complete at this point
    obj = -T[-1, -1]

    if status in (0, 1):
        if disp:
            print(messages[status])
            print("         Current function value: {0: <12.6f}".format(obj))
            print("         Iterations: {0:d}".format(nit2))
    else:
        if disp:
            print(messages[status])
            print("         Iterations: {0:d}".format(nit2))

    return OptimizeResult(x=x, fun=obj, nit=int(nit2), status=status,
                          slack=slack, message=messages[status],
                          success=(status == 0))


def linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None,
            bounds=None, method='simplex', callback=None,
            options=None):
    """
    Minimize a linear objective function subject to linear
    equality and inequality constraints.

    Linear Programming is intended to solve the following problem form::

        Minimize:     c^T * x

        Subject to:   A_ub * x <= b_ub
                      A_eq * x == b_eq

    Parameters
    ----------
    c : array_like
        Coefficients of the linear objective function to be minimized.
    A_ub : array_like, optional
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``.
    b_ub : array_like, optional
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : array_like, optional
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``.
    b_eq : array_like, optional
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for one of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative)
        If a sequence containing a single tuple is provided, then ``min`` and
        ``max`` will be applied to all variables in the problem.
    method : str, optional
        Type of solver.  :ref:`'simplex' <optimize.linprog-simplex>`
        and :ref:`'interior-point' <optimize.linprog-interior-point>`
        are supported.
    callback : callable, optional (simplex only)
        If a callback function is provide, it will be called within each
        iteration of the simplex algorithm. The callback must have the
        signature ``callback(xk, **kwargs)`` where ``xk`` is the current
        solution vector and ``kwargs`` is a dictionary containing the
        following::

            "tableau" : The current Simplex algorithm tableau
            "nit" : The current iteration.
            "pivot" : The pivot (row, column) used for the next iteration.
            "phase" : Whether the algorithm is in Phase 1 or Phase 2.
            "basis" : The indices of the columns of the basic variables.

    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.

    Returns
    -------
    A `scipy.optimize.OptimizeResult` consisting of the following fields:

        x : ndarray
            The independent variable vector which optimizes the linear
            programming problem.
        fun : float
            Value of the objective function.
        slack : ndarray
            The values of the slack variables.  Each slack variable corresponds
            to an inequality constraint.  If the slack is zero, then the
            corresponding constraint is active.
        success : bool
            Returns True if the algorithm succeeded in finding an optimal
            solution.
        status : int
            An integer representing the exit status of the optimization::

                 0 : Optimization terminated successfully
                 1 : Iteration limit reached
                 2 : Problem appears to be infeasible
                 3 : Problem appears to be unbounded

        nit : int
            The number of iterations performed.
        message : str
            A string descriptor of the exit status of the optimization.

    See Also
    --------
    show_options : Additional options accepted by the solvers

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method
    is :ref:`Simplex <optimize.linprog-simplex>`.
    :ref:`Interior point <optimize.linprog-interior-point>` is also available.

    Method *simplex* uses the simplex algorithm (as it relates to linear
    programming, NOT the Nelder-Mead simplex) [1]_, [2]_. This algorithm
    should be reasonably reliable and fast for small problems.

    .. versionadded:: 0.15.0

    Method *interior-point* uses the primal-dual path following algorithm
    as outlined in [4]_. This algorithm is intended to provide a faster
    and more reliable alternative to *simplex*, especially for large,
    sparse problems. Note, however, that the solution returned may be slightly
    less accurate than that of the simplex method and may not correspond with a
    vertex of the polytope defined by the constraints.

    References
    ----------
    .. [1] Dantzig, George B., Linear programming and extensions. Rand
           Corporation Research Study Princeton Univ. Press, Princeton, NJ,
           1963
    .. [2] Hillier, S.H. and Lieberman, G.J. (1995), "Introduction to
           Mathematical Programming", McGraw-Hill, Chapter 4.
    .. [3] Bland, Robert G. New finite pivoting rules for the simplex method.
           Mathematics of Operations Research (2), 1977: pp. 103-107.
    .. [4] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.
    .. [5] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [6] Freund, Robert M. "Primal-Dual Interior-Point Methods for Linear
           Programming based on Newton's Method." Unpublished Course Notes,
           March 2004. Available 2/25/2017 at
           https://ocw.mit.edu/courses/sloan-school-of-management/15-084j-nonlinear-programming-spring-2004/lecture-notes/lec14_int_pt_mthd.pdf
    .. [7] Fourer, Robert. "Solving Linear Programs by Interior-Point Methods."
           Unpublished Course Notes, August 26, 2005. Available 2/25/2017 at
           http://www.4er.org/CourseNotes/Book%20B/B-III.pdf
    .. [8] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
           programming." Mathematical Programming 71.2 (1995): 221-245.
    .. [9] Bertsimas, Dimitris, and J. Tsitsiklis. "Introduction to linear
           programming." Athena Scientific 1 (1997): 997.
    .. [10] Andersen, Erling D., et al. Implementation of interior point
            methods for large scale linear programming. HEC/Universite de
            Geneve, 1996.

    Examples
    --------
    Consider the following problem:

    Minimize: f = -1*x[0] + 4*x[1]

    Subject to: -3*x[0] + 1*x[1] <= 6
                 1*x[0] + 2*x[1] <= 4
                            x[1] >= -3

    where:  -inf <= x[0] <= inf

    This problem deviates from the standard linear programming problem.
    In standard form, linear programming problems assume the variables x are
    non-negative.  Since the variables don't have standard bounds where
    0 <= x <= inf, the bounds of the variables must be explicitly set.

    There are two upper-bound constraints, which can be expressed as

    dot(A_ub, x) <= b_ub

    The input for this problem is as follows:

    >>> c = [-1, 4]
    >>> A = [[-3, 1], [1, 2]]
    >>> b = [6, 4]
    >>> x0_bounds = (None, None)
    >>> x1_bounds = (-3, None)
    >>> from scipy.optimize import linprog
    >>> res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds),
    ...               options={"disp": True})
    Optimization terminated successfully.
         Current function value: -22.000000
         Iterations: 1
    >>> print(res)
         fun: -22.0
     message: 'Optimization terminated successfully.'
         nit: 1
       slack: array([39.,  0.])
      status: 0
     success: True
           x: array([10., -3.])

    Note the actual objective value is 11.428571.  In this case we minimized
    the negative of the objective function.

    """
    meth = method.lower()
    if options is None:
        options = {}

    if meth == 'simplex':
        return _linprog_simplex(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                                bounds=bounds, callback=callback, **options)
    elif meth == 'interior-point':
        return _linprog_ip(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                           bounds=bounds, callback=callback, **options)
    else:
        raise ValueError('Unknown solver %s' % method)
