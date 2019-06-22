"""
Method agnostic utility functions for linear progamming
"""

import numpy as np
import scipy.sparse as sps
from warnings import warn
from .optimize import OptimizeWarning
from scipy.optimize._remove_redundancy import (
    _remove_redundancy, _remove_redundancy_sparse, _remove_redundancy_dense
    )


def _check_sparse_inputs(options, A_ub, A_eq):
    """
    Check the provided ``A_ub`` and ``A_eq`` matrices conform to the specified
    optional sparsity variables.

    Parameters
    ----------
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.

    Returns
    -------
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.
    """
    # This is an undocumented option for unit testing sparse presolve
    _sparse_presolve = options.pop('_sparse_presolve', False)
    if _sparse_presolve and A_eq is not None:
        A_eq = sps.coo_matrix(A_eq)
    if _sparse_presolve and A_ub is not None:
        A_ub = sps.coo_matrix(A_ub)

    sparse = options.get('sparse', False)
    if not sparse and (sps.issparse(A_eq) or sps.issparse(A_ub)):
        options['sparse'] = True
        warn("Sparse constraint matrix detected; setting 'sparse':True.",
             OptimizeWarning)
    return options, A_ub, A_eq


def _format_A_constraints(A, n_x, sparse_lhs=False):
    """Format the left hand side of the constraints to a 2D array

    Parameters
    ----------
    A : 2D array
        2D array such that ``A @ x`` gives the values of the upper-bound
        (in)equality constraints at ``x``.
    n_x : int
        The number of variables in the linear programming problem.
    sparse_lhs : bool
        Whether either of `A_ub` or `A_eq` are sparse. If true return a
        coo_matrix instead of a numpy array.

    Returns
    -------
    np.ndarray or sparse.coo_matrix
        2D array such that ``A @ x`` gives the values of the upper-bound
        (in)equality constraints at ``x``.

    """
    if sparse_lhs:
        return sps.coo_matrix(
            (0, n_x) if A is None else A, dtype=float, copy=True
        )
    elif A is None:
        return np.zeros((0, n_x), dtype=float)
    else:
        return np.array(A, dtype=float, copy=True)


def _format_b_constraints(b):
    """Format the upper bounds of the constraints to a 1D array

    Parameters
    ----------
    b : 1D array
        1D array of values representing the upper-bound of each (in)equality
        constraint (row) in ``A``.

    Returns
    -------
    1D np.array
        1D array of values representing the upper-bound of each (in)equality
        constraint (row) in ``A``.

    """
    if b is None:
        return np.array([], dtype=float)
    b = np.array(b, dtype=float, copy=True).squeeze()
    return b if b.size != 1 else b.reshape((-1))


def _clean_inputs(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
                  x0=None):
    """
    Given user inputs for a linear programming problem, return the
    objective vector, upper bound constraints, equality constraints,
    and simple bounds in a preferred format.

    Parameters
    ----------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for one of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative).
        If a sequence containing a single tuple is provided, then ``min`` and
        ``max`` will be applied to all variables in the problem.
    x0 : 1D array, optional
        Starting values of the independent variables, which will be refined by
        the optimization algorithm.

    Returns
    -------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative).
    x0 : 1D array, optional
        Starting values of the independent variables, which will be refined by
        the optimization algorithm.
    """
    if c is None:
        raise TypeError

    try:
        c = np.array(c, dtype=np.float, copy=True).squeeze()
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: c must be a 1D array of numerical "
            "coefficients")
    else:
        # If c is a single value, convert it to a 1D array.
        if c.size == 1:
            c = c.reshape((-1))

        n_x = len(c)
        if n_x == 0 or len(c.shape) != 1:
            raise ValueError(
                "Invalid input for linprog: c must be a 1D array and must "
                "not have more than one non-singleton dimension")
        if not(np.isfinite(c).all()):
            raise ValueError(
                "Invalid input for linprog: c must not contain values "
                "inf, nan, or None")

    sparse_lhs = sps.issparse(A_eq) or sps.issparse(A_ub)
    try:
        A_ub = _format_A_constraints(A_ub, n_x, sparse_lhs=sparse_lhs)
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: A_ub must be a 2D array "
            "of numerical values")
    else:
        n_ub = A_ub.shape[0]
        if len(A_ub.shape) != 2 or A_ub.shape[1] != n_x:
            raise ValueError(
                "Invalid input for linprog: A_ub must have exactly two "
                "dimensions, and the number of columns in A_ub must be "
                "equal to the size of c")
        if (sps.issparse(A_ub) and not np.isfinite(A_ub.data).all()
                or not sps.issparse(A_ub) and not np.isfinite(A_ub).all()):
            raise ValueError(
                "Invalid input for linprog: A_ub must not contain values "
                "inf, nan, or None")

    try:
        b_ub = _format_b_constraints(b_ub)
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: b_ub must be a 1D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_ub")
    else:
        if b_ub.shape != (n_ub,):
            raise ValueError(
                "Invalid input for linprog: b_ub must be a 1D array; b_ub "
                "must not have more than one non-singleton dimension and "
                "the number of rows in A_ub must equal the number of values "
                "in b_ub")
        if not(np.isfinite(b_ub).all()):
            raise ValueError(
                "Invalid input for linprog: b_ub must not contain values "
                "inf, nan, or None")

    try:
        A_eq = _format_A_constraints(A_eq, n_x, sparse_lhs=sparse_lhs)
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: A_eq must be a 2D array "
            "of numerical values")
    else:
        n_eq = A_eq.shape[0]
        if len(A_eq.shape) != 2 or A_eq.shape[1] != n_x:
            raise ValueError(
                "Invalid input for linprog: A_eq must have exactly two "
                "dimensions, and the number of columns in A_eq must be "
                "equal to the size of c")

        if (sps.issparse(A_eq) and not np.isfinite(A_eq.data).all()
                or not sps.issparse(A_eq) and not np.isfinite(A_eq).all()):
            raise ValueError(
                "Invalid input for linprog: A_eq must not contain values "
                "inf, nan, or None")

    try:
        b_eq = _format_b_constraints(b_eq)
    except ValueError:
        raise TypeError(
            "Invalid input for linprog: b_eq must be a 1D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_eq")
    else:
        if b_eq.shape != (n_eq,):
            raise ValueError(
                "Invalid input for linprog: b_eq must be a 1D array; b_eq "
                "must not have more than one non-singleton dimension and "
                "the number of rows in A_eq must equal the number of values "
                "in b_eq")
        if not(np.isfinite(b_eq).all()):
            raise ValueError(
                "Invalid input for linprog: b_eq must not contain values "
                "inf, nan, or None")

    # x0 gives a (optional) starting solution to the solver. If x0 is None,
    # skip the checks. Initial solution will be generated automatically.
    if x0 is not None:
        try:
            x0 = np.array(x0, dtype=float, copy=True).squeeze()
        except ValueError:
            raise TypeError(
                "Invalid input for linprog: x0 must be a 1D array of "
                "numerical coefficients")
        if x0.ndim == 0:
            x0 = x0.reshape((-1))
        if len(x0) == 0 or x0.ndim != 1:
            raise ValueError(
                "Invalid input for linprog: x0 should be a 1D array; it "
                "must not have more than one non-singleton dimension")
        if not x0.size == c.size:
            raise ValueError(
                "Invalid input for linprog: x0 and c should contain the "
                "same number of elements")
        if not np.isfinite(x0).all():
            raise ValueError(
            "Invalid input for linprog: x0 must not contain values "
            "inf, nan, or None")

    # "If a sequence containing a single tuple is provided, then min and max
    # will be applied to all variables in the problem."
    # linprog doesn't treat this right: it didn't accept a list with one tuple
    # in it
    try:
        if isinstance(bounds, str):
            raise TypeError
        if bounds is None or len(bounds) == 0:
            bounds = [(0, None)] * n_x
        elif len(bounds) == 1:
            b = bounds[0]
            if len(b) != 2:
                raise ValueError(
                    "Invalid input for linprog: exactly one lower bound and "
                    "one upper bound must be specified for each element of x")
            bounds = [b] * n_x
        elif len(bounds) == n_x:
            try:
                len(bounds[0])
            except BaseException:
                bounds = [(bounds[0], bounds[1])] * n_x
            for i, b in enumerate(bounds):
                if len(b) != 2:
                    raise ValueError(
                        "Invalid input for linprog, bound " +
                        str(i) +
                        " " +
                        str(b) +
                        ": exactly one lower bound and one upper bound must "
                        "be specified for each element of x")
        elif (len(bounds) == 2 and np.isreal(bounds[0])
                and np.isreal(bounds[1])):
            bounds = [(bounds[0], bounds[1])] * n_x
        else:
            raise ValueError(
                "Invalid input for linprog: exactly one lower bound and one "
                "upper bound must be specified for each element of x")

        clean_bounds = []  # also creates a copy so user's object isn't changed
        for i, b in enumerate(bounds):
            if b[0] is not None and b[1] is not None and b[0] > b[1]:
                raise ValueError(
                    "Invalid input for linprog, bound " +
                    str(i) +
                    " " +
                    str(b) +
                    ": a lower bound must be less than or equal to the "
                    "corresponding upper bound")
            if b[0] == np.inf:
                raise ValueError(
                    "Invalid input for linprog, bound " +
                    str(i) +
                    " " +
                    str(b) +
                    ": infinity is not a valid lower bound")
            if b[1] == -np.inf:
                raise ValueError(
                    "Invalid input for linprog, bound " +
                    str(i) +
                    " " +
                    str(b) +
                    ": negative infinity is not a valid upper bound")
            lb = float(b[0]) if b[0] is not None and b[0] != -np.inf else None
            ub = float(b[1]) if b[1] is not None and b[1] != np.inf else None
            clean_bounds.append((lb, ub))
        bounds = clean_bounds
    except ValueError as e:
        if "could not convert string to float" in e.args[0]:
            raise TypeError
        else:
            raise e
    except TypeError as e:
        print(e)
        raise TypeError(
            "Invalid input for linprog: bounds must be a sequence of "
            "(min,max) pairs, each defining bounds on an element of x ")

    return c, A_ub, b_ub, A_eq, b_eq, bounds, x0


def _presolve(c, A_ub, b_ub, A_eq, b_eq, bounds, x0, rr, tol=1e-9):
    """
    Given inputs for a linear programming problem in preferred format,
    presolve the problem: identify trivial infeasibilities, redundancies,
    and unboundedness, tighten bounds where possible, and eliminate fixed
    variables.

    Parameters
    ----------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction.
    x0 : 1D array, optional
        Starting values of the independent variables, which will be refined by
        the optimization algorithm.
    rr : bool
        If ``True`` attempts to eliminate any redundant rows in ``A_eq``.
        Set False if ``A_eq`` is known to be of full row rank, or if you are
        looking for a potential speedup (at the expense of reliability).
    tol : float
        The tolerance which determines when a solution is "close enough" to
        zero in Phase 1 to be considered a basic feasible solution or close
        enough to positive to serve as an optimal solution.

    Returns
    -------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    c0 : 1D array
        Constant term in objective function due to fixed (and eliminated)
        variables.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. Bounds have been
        tightened where possible.
    x : 1D array
        Solution vector (when the solution is trivial and can be determined
        in presolve)
    x0 : 1D array
        Starting values of the independent variables, which will be refined by
        the optimization algorithm (if solution is not determined in presolve)
    undo: list of tuples
        (index, value) pairs that record the original index and fixed value
        for each variable removed from the problem
    complete: bool
        Whether the solution is complete (solved or determined to be infeasible
        or unbounded in presolve)
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [5] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [8] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
           programming." Mathematical Programming 71.2 (1995): 221-245.

    """
    # ideas from Reference [5] by Andersen and Andersen
    # however, unlike the reference, this is performed before converting
    # problem to standard form
    # There are a few advantages:
    #  * artificial variables have not been added, so matrices are smaller
    #  * bounds have not been converted to constraints yet. (It is better to
    #    do that after presolve because presolve may adjust the simple bounds.)
    # There are many improvements that can be made, namely:
    #  * implement remaining checks from [5]
    #  * loop presolve until no additional changes are made
    #  * implement additional efficiency improvements in redundancy removal [2]

    undo = []               # record of variables eliminated from problem
    # constant term in cost function may be added if variables are eliminated
    c0 = 0
    complete = False        # complete is True if detected infeasible/unbounded
    x = np.zeros(c.shape)   # this is solution vector if completed in presolve

    status = 0              # all OK unless determined otherwise
    message = ""

    # Standard form for bounds (from _clean_inputs) is list of tuples
    # but numpy array is more convenient here
    # In retrospect, numpy array should have been the standard
    bounds = np.array(bounds)
    lb = bounds[:, 0]
    ub = bounds[:, 1]
    lb[np.equal(lb, None)] = -np.inf
    ub[np.equal(ub, None)] = np.inf
    bounds = bounds.astype(float)
    lb = lb.astype(float)
    ub = ub.astype(float)

    m_eq, n = A_eq.shape
    m_ub, n = A_ub.shape

    if (sps.issparse(A_eq)):
        A_eq = A_eq.tolil()
        A_ub = A_ub.tolil()

        def where(A):
            return A.nonzero()

        vstack = sps.vstack
    else:
        where = np.where
        vstack = np.vstack

    # zero row in equality constraints
    zero_row = np.array(np.sum(A_eq != 0, axis=1) == 0).flatten()
    if np.any(zero_row):
        if np.any(
            np.logical_and(
                zero_row,
                np.abs(b_eq) > tol)):  # test_zero_row_1
            # infeasible if RHS is not zero
            status = 2
            message = ("The problem is (trivially) infeasible due to a row "
                       "of zeros in the equality constraint matrix with a "
                       "nonzero corresponding constraint value.")
            complete = True
            return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                    x, x0, undo, complete, status, message)
        else:  # test_zero_row_2
            # if RHS is zero, we can eliminate this equation entirely
            A_eq = A_eq[np.logical_not(zero_row), :]
            b_eq = b_eq[np.logical_not(zero_row)]

    # zero row in inequality constraints
    zero_row = np.array(np.sum(A_ub != 0, axis=1) == 0).flatten()
    if np.any(zero_row):
        if np.any(np.logical_and(zero_row, b_ub < -tol)):  # test_zero_row_1
            # infeasible if RHS is less than zero (because LHS is zero)
            status = 2
            message = ("The problem is (trivially) infeasible due to a row "
                       "of zeros in the equality constraint matrix with a "
                       "nonzero corresponding  constraint value.")
            complete = True
            return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                    x, x0, undo, complete, status, message)
        else:  # test_zero_row_2
            # if LHS is >= 0, we can eliminate this constraint entirely
            A_ub = A_ub[np.logical_not(zero_row), :]
            b_ub = b_ub[np.logical_not(zero_row)]

    # zero column in (both) constraints
    # this indicates that a variable isn't constrained and can be removed
    A = vstack((A_eq, A_ub))
    if A.shape[0] > 0:
        zero_col = np.array(np.sum(A != 0, axis=0) == 0).flatten()
        # variable will be at upper or lower bound, depending on objective
        x[np.logical_and(zero_col, c < 0)] = ub[
            np.logical_and(zero_col, c < 0)]
        x[np.logical_and(zero_col, c > 0)] = lb[
            np.logical_and(zero_col, c > 0)]
        if np.any(np.isinf(x)):  # if an unconstrained variable has no bound
            status = 3
            message = ("If feasible, the problem is (trivially) unbounded "
                       "due  to a zero column in the constraint matrices. If "
                       "you wish to check whether the problem is infeasible, "
                       "turn presolve off.")
            complete = True
            return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                    x, x0, undo, complete, status, message)
        # variables will equal upper/lower bounds will be removed later
        lb[np.logical_and(zero_col, c < 0)] = ub[
            np.logical_and(zero_col, c < 0)]
        ub[np.logical_and(zero_col, c > 0)] = lb[
            np.logical_and(zero_col, c > 0)]

    # row singleton in equality constraints
    # this fixes a variable and removes the constraint
    singleton_row = np.array(np.sum(A_eq != 0, axis=1) == 1).flatten()
    rows = where(singleton_row)[0]
    cols = where(A_eq[rows, :])[1]
    if len(rows) > 0:
        for row, col in zip(rows, cols):
            val = b_eq[row] / A_eq[row, col]
            if not lb[col] - tol <= val <= ub[col] + tol:
                # infeasible if fixed value is not within bounds
                status = 2
                message = ("The problem is (trivially) infeasible because a "
                           "singleton row in the equality constraints is "
                           "inconsistent with the bounds.")
                complete = True
                return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                        x, x0, undo, complete, status, message)
            else:
                # sets upper and lower bounds at that fixed value - variable
                # will be removed later
                lb[col] = val
                ub[col] = val
        A_eq = A_eq[np.logical_not(singleton_row), :]
        b_eq = b_eq[np.logical_not(singleton_row)]

    # row singleton in inequality constraints
    # this indicates a simple bound and the constraint can be removed
    # simple bounds may be adjusted here
    # After all of the simple bound information is combined here, get_Abc will
    # turn the simple bounds into constraints
    singleton_row = np.array(np.sum(A_ub != 0, axis=1) == 1).flatten()
    cols = where(A_ub[singleton_row, :])[1]
    rows = where(singleton_row)[0]
    if len(rows) > 0:
        for row, col in zip(rows, cols):
            val = b_ub[row] / A_ub[row, col]
            if A_ub[row, col] > 0:  # upper bound
                if val < lb[col] - tol:  # infeasible
                    complete = True
                elif val < ub[col]:  # new upper bound
                    ub[col] = val
            else:  # lower bound
                if val > ub[col] + tol:  # infeasible
                    complete = True
                elif val > lb[col]:  # new lower bound
                    lb[col] = val
            if complete:
                status = 2
                message = ("The problem is (trivially) infeasible because a "
                           "singleton row in the upper bound constraints is "
                           "inconsistent with the bounds.")
                return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                        x, x0, undo, complete, status, message)
        A_ub = A_ub[np.logical_not(singleton_row), :]
        b_ub = b_ub[np.logical_not(singleton_row)]

    # identical bounds indicate that variable can be removed
    i_f = np.abs(lb - ub) < tol   # indices of "fixed" variables
    i_nf = np.logical_not(i_f)  # indices of "not fixed" variables

    # test_bounds_equal_but_infeasible
    if np.all(i_f):  # if bounds define solution, check for consistency
        residual = b_eq - A_eq.dot(lb)
        slack = b_ub - A_ub.dot(lb)
        if ((A_ub.size > 0 and np.any(slack < 0)) or
                (A_eq.size > 0 and not np.allclose(residual, 0))):
            status = 2
            message = ("The problem is (trivially) infeasible because the "
                       "bounds fix all variables to values inconsistent with "
                       "the constraints")
            complete = True
            return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                    x, x0, undo, complete, status, message)

    ub_mod = ub
    lb_mod = lb
    if np.any(i_f):
        c0 += c[i_f].dot(lb[i_f])
        b_eq = b_eq - A_eq[:, i_f].dot(lb[i_f])
        b_ub = b_ub - A_ub[:, i_f].dot(lb[i_f])
        c = c[i_nf]
        x = x[i_nf]
        # user guess x0 stays separate from presolve solution x
        if x0 is not None:
            x0 = x0[i_nf]
        A_eq = A_eq[:, i_nf]
        A_ub = A_ub[:, i_nf]
        # record of variables to be added back in
        undo = [np.nonzero(i_f)[0], lb[i_f]]
        # don't remove these entries from bounds; they'll be used later.
        # but we _also_ need a version of the bounds with these removed
        lb_mod = lb[i_nf]
        ub_mod = ub[i_nf]

    # no constraints indicates that problem is trivial
    if A_eq.size == 0 and A_ub.size == 0:
        b_eq = np.array([])
        b_ub = np.array([])
        # test_empty_constraint_1
        if c.size == 0:
            status = 0
            message = ("The solution was determined in presolve as there are "
                       "no non-trivial constraints.")
        elif (np.any(np.logical_and(c < 0, ub_mod == np.inf)) or
              np.any(np.logical_and(c > 0, lb_mod == -np.inf))):
            # test_no_constraints()
            # test_unbounded_no_nontrivial_constraints_1
            # test_unbounded_no_nontrivial_constraints_2
            status = 3
            message = ("The problem is (trivially) unbounded "
                       "because there are no non-trivial constraints and "
                       "a) at least one decision variable is unbounded "
                       "above and its corresponding cost is negative, or "
                       "b) at least one decision variable is unbounded below "
                       "and its corresponding cost is positive. ")
        else:  # test_empty_constraint_2
            status = 0
            message = ("The solution was determined in presolve as there are "
                       "no non-trivial constraints.")
        complete = True
        x[c < 0] = ub_mod[c < 0]
        x[c > 0] = lb_mod[c > 0]
        # where c is zero, set x to a finite bound or zero
        x_zero_c = ub_mod[c == 0]
        x_zero_c[np.isinf(x_zero_c)] = ub_mod[c == 0][np.isinf(x_zero_c)]
        x_zero_c[np.isinf(x_zero_c)] = 0
        x[c == 0] = x_zero_c
        # if this is not the last step of presolve, should convert bounds back
        # to array and return here

    # *sigh* - convert bounds back to their standard form (list of tuples)
    # again, in retrospect, numpy array would be standard form
    lb[np.equal(lb, -np.inf)] = None
    ub[np.equal(ub, np.inf)] = None
    bounds = np.hstack((lb[:, np.newaxis], ub[:, np.newaxis]))
    bounds = bounds.tolist()
    for i, row in enumerate(bounds):
        for j, col in enumerate(row):
            if str(col) == "nan":
                # comparing col to float("nan") and np.nan doesn't work.
                # should use np.isnan
                bounds[i][j] = None

    # remove redundant (linearly dependent) rows from equality constraints
    n_rows_A = A_eq.shape[0]
    redundancy_warning = ("A_eq does not appear to be of full row rank. To "
                          "improve performance, check the problem formulation "
                          "for redundant equality constraints.")
    if (sps.issparse(A_eq)):
        if rr and A_eq.size > 0:  # TODO: Fast sparse rank check?
            A_eq, b_eq, status, message = _remove_redundancy_sparse(A_eq, b_eq)
            if A_eq.shape[0] < n_rows_A:
                warn(redundancy_warning, OptimizeWarning)
            if status != 0:
                complete = True
        return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
                x, x0, undo, complete, status, message)

    # This is a wild guess for which redundancy removal algorithm will be
    # faster. More testing would be good.
    small_nullspace = 5
    if rr and A_eq.size > 0:
        try:  # TODO: instead use results of first SVD in _remove_redundancy
            rank = np.linalg.matrix_rank(A_eq)
        except Exception:  # oh well, we'll have to go with _remove_redundancy_dense
            rank = 0
    if rr and A_eq.size > 0 and rank < A_eq.shape[0]:
        warn(redundancy_warning, OptimizeWarning)
        dim_row_nullspace = A_eq.shape[0]-rank
        if dim_row_nullspace <= small_nullspace:
            A_eq, b_eq, status, message = _remove_redundancy(A_eq, b_eq)
        if dim_row_nullspace > small_nullspace or status == 4:
            A_eq, b_eq, status, message = _remove_redundancy_dense(A_eq, b_eq)
        if A_eq.shape[0] < rank:
            message = ("Due to numerical issues, redundant equality "
                       "constraints could not be removed automatically. "
                       "Try providing your constraint matrices as sparse "
                       "matrices to activate sparse presolve, try turning "
                       "off redundancy removal, or try turning off presolve "
                       "altogether.")
            status = 4
        if status != 0:
            complete = True
    return (c, c0, A_ub, b_ub, A_eq, b_eq, bounds,
            x, x0, undo, complete, status, message)


def _parse_linprog(c, A_ub, b_ub, A_eq, b_eq, bounds, options, x0):
    """
    Parse the provided linear programming problem

    ``_parse_linprog`` employs two main steps ``_check_sparse_inputs`` and
    ``_clean_inputs``. ``_check_sparse_inputs`` checks for sparsity in the
    provided constraints (``A_ub`` and ``A_eq) and if these match the provided
    sparsity optional values.

    ``_clean inputs`` checks of the provided inputs. If no violations are
    identified the objective vector, upper bound constraints, equality
    constraints, and simple bounds are returned in the expected format.

    Parameters
    ----------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for one of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative). If a sequence containing a
        single tuple is provided, then ``min`` and ``max`` will be applied to
        all variables in the problem.
    options : dict
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.
    x0 : 1D array, optional
        Starting values of the independent variables, which will be refined by
        the optimization algorithm. Currently compatible only with the
        'revised simplex' method, and only if x0 is a basic feasible solution
        of the problem.

    Returns
    -------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for one of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative).
        If a sequence containing a single tuple is provided, then ``min`` and
        ``max`` will be applied to all variables in the problem.
    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:

            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        For method-specific options, see :func:`show_options('linprog')`.
    x0 : 1D array, optional
        Starting values of the independent variables, which will be refined by
        the optimization algorithm. Currently compatible only with the
        'revised simplex' method, and only if x0 is a basic feasible solution
        of the problem.
    """
    if options is None:
        options = {}

    solver_options = {k: v for k, v in options.items()}
    solver_options, A_ub, A_eq = _check_sparse_inputs(solver_options, A_ub, A_eq)
    # Convert lists to numpy arrays, etc...
    c, A_ub, b_ub, A_eq, b_eq, bounds, x0 = _clean_inputs(
        c, A_ub, b_ub, A_eq, b_eq, bounds, x0)
    return c, A_ub, b_ub, A_eq, b_eq, bounds, solver_options, x0


def _get_Abc(c, c0=0, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
             x0=None, undo=[]):
    """
    Given a linear programming problem of the form:

    Minimize::

        c @ x

    Subject to::

        A_ub @ x <= b_ub
        A_eq @ x == b_eq
         lb <= x <= ub

    where ``lb = 0`` and ``ub = None`` unless set in ``bounds``.

    Return the problem in standard form:

    Minimize::

        c @ x

    Subject to::

        A @ x == b
            x >= 0

    by adding slack variables and making variable substitutions as necessary.

    Parameters
    ----------
    c : 1D array
        Coefficients of the linear objective function to be minimized.
        Components corresponding with fixed variables have been eliminated.
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. Bounds have been
        tightened where possible.
    x0 : 1D array
        Starting values of the independent variables, which will be refined by
        the optimization algorithm
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem

    Returns
    -------
    A : 2D array
        2D array such that ``A`` @ ``x``, gives the values of the equality
        constraints at ``x``.
    b : 1D array
        1D array of values representing the RHS of each equality constraint
        (row) in A (for standard form problem).
    c : 1D array
        Coefficients of the linear objective function to be minimized (for
        standard form problem).
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.
    x0 : 1D array
        Starting values of the independent variables, which will be refined by
        the optimization algorithm

    References
    ----------
    .. [9] Bertsimas, Dimitris, and J. Tsitsiklis. "Introduction to linear
           programming." Athena Scientific 1 (1997): 997.

    """

    if sps.issparse(A_eq):
        sparse = True
        A_eq = sps.lil_matrix(A_eq)
        A_ub = sps.lil_matrix(A_ub)

        def hstack(blocks):
            return sps.hstack(blocks, format="lil")

        def vstack(blocks):
            return sps.vstack(blocks, format="lil")

        zeros = sps.lil_matrix
        eye = sps.eye
    else:
        sparse = False
        hstack = np.hstack
        vstack = np.vstack
        zeros = np.zeros
        eye = np.eye

    fixed_x = set()
    if len(undo) > 0:
        # these are indices of variables removed from the problem
        # however, their bounds are still part of the bounds list
        fixed_x = set(undo[0])
    # they are needed elsewhere, but not here
    bounds = [bounds[i] for i in range(len(bounds)) if i not in fixed_x]
    # in retrospect, the standard form of bounds should have been an n x 2
    # array. maybe change it someday.

    # modify problem such that all variables have only non-negativity bounds

    bounds = np.array(bounds)
    lbs = bounds[:, 0]
    ubs = bounds[:, 1]
    m_ub, n_ub = A_ub.shape

    lb_none = np.equal(lbs, None)
    ub_none = np.equal(ubs, None)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)

    # if preprocessing is on, lb == ub can't happen
    # if preprocessing is off, then it would be best to convert that
    # to an equality constraint, but it's tricky to make the other
    # required modifications from inside here.

    # unbounded below: substitute xi = -xi' (unbounded above)
    l_nolb_someub = np.logical_and(lb_none, ub_some)
    i_nolb = np.nonzero(l_nolb_someub)[0]
    lbs[l_nolb_someub], ubs[l_nolb_someub] = (
        -ubs[l_nolb_someub], lbs[l_nolb_someub])
    lb_none = np.equal(lbs, None)
    ub_none = np.equal(ubs, None)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)
    c[i_nolb] *= -1
    if x0 is not None:
        x0[i_nolb] *= -1
    if len(i_nolb) > 0:
        if A_ub.shape[0] > 0:  # sometimes needed for sparse arrays... weird
            A_ub[:, i_nolb] *= -1
        if A_eq.shape[0] > 0:
            A_eq[:, i_nolb] *= -1

    # upper bound: add inequality constraint
    i_newub = np.nonzero(ub_some)[0]
    ub_newub = ubs[ub_some]
    n_bounds = np.count_nonzero(ub_some)
    A_ub = vstack((A_ub, zeros((n_bounds, A_ub.shape[1]))))
    b_ub = np.concatenate((b_ub, np.zeros(n_bounds)))
    A_ub[range(m_ub, A_ub.shape[0]), i_newub] = 1
    b_ub[m_ub:] = ub_newub

    A1 = vstack((A_ub, A_eq))
    b = np.concatenate((b_ub, b_eq))
    c = np.concatenate((c, np.zeros((A_ub.shape[0],))))
    if x0 is not None:
        x0 = np.concatenate((x0, np.zeros((A_ub.shape[0],))))
    # unbounded: substitute xi = xi+ + xi-
    l_free = np.logical_and(lb_none, ub_none)
    i_free = np.nonzero(l_free)[0]
    n_free = len(i_free)
    A1 = hstack((A1, zeros((A1.shape[0], n_free))))
    c = np.concatenate((c, np.zeros(n_free)))
    if x0 is not None:
        x0 = np.concatenate((x0, np.zeros(n_free)))
    A1[:, range(n_ub, A1.shape[1])] = -A1[:, i_free]
    c[np.arange(n_ub, A1.shape[1])] = -c[i_free]
    if x0 is not None:
        i_free_neg = x0[i_free] < 0
        x0[np.arange(n_ub, A1.shape[1])[i_free_neg]] = -x0[i_free[i_free_neg]]
        x0[i_free[i_free_neg]] = 0

    # add slack variables
    A2 = vstack([eye(A_ub.shape[0]), zeros((A_eq.shape[0], A_ub.shape[0]))])
    A = hstack([A1, A2])

    # lower bound: substitute xi = xi' + lb
    # now there is a constant term in objective
    i_shift = np.nonzero(lb_some)[0]
    lb_shift = lbs[lb_some].astype(float)
    c0 += np.sum(lb_shift * c[i_shift])
    if sparse:
        b = b.reshape(-1, 1)
        A = A.tocsc()
        b -= (A[:, i_shift] * sps.diags(lb_shift)).sum(axis=1)
        b = b.ravel()
    else:
        b -= (A[:, i_shift] * lb_shift).sum(axis=1)
    if x0 is not None:
        x0[i_shift] -= lb_shift

    return A, b, c, c0, x0


def _display_summary(message, status, fun, iteration):
    """
    Print the termination summary of the linear program

    Parameters
    ----------
    message : str
            A string descriptor of the exit status of the optimization.
    status : int
        An integer representing the exit status of the optimization::

                0 : Optimization terminated successfully
                1 : Iteration limit reached
                2 : Problem appears to be infeasible
                3 : Problem appears to be unbounded
                4 : Serious numerical difficulties encountered

    fun : float
        Value of the objective function.
    iteration : iteration
        The number of iterations performed.
    """
    print(message)
    if status in (0, 1):
        print("         Current function value: {0: <12.6f}".format(fun))
    print("         Iterations: {0:d}".format(iteration))


def _postsolve(x, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
               complete=False, undo=[], tol=1e-8, copy=False):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1D array
        Solution vector to the standard-form problem.
    c : 1D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        Bounds, as modified in presolve
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack : 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    lb : 1D array
        The lower bound constraints on the original variables
    ub: 1D array
        The upper bound constraints on the original variables
    """
    # note that all the inputs are the ORIGINAL, unmodified versions
    # no rows, columns have been removed
    # the only exception is bounds; it has been modified
    # we need these modified values to undo the variable substitutions
    # in retrospect, perhaps this could have been simplified if the "undo"
    # variable also contained information for undoing variable substitutions

    n_x = len(c)

    # we don't have to undo variable substitutions for fixed variables that
    # were removed from the problem
    no_adjust = set()

    # if there were variables removed from the problem, add them back into the
    # solution vector
    if len(undo) > 0:
        no_adjust = set(undo[0])
        x = x.tolist()
        for i, val in zip(undo[0], undo[1]):
            x.insert(i, val)
        copy = True
    if copy:
        x = np.array(x, copy=True)

    # now undo variable substitutions
    # if "complete", problem was solved in presolve; don't do anything here
    if not complete and bounds is not None:  # bounds are never none, probably
        n_unbounded = 0
        for i, b in enumerate(bounds):
            if i in no_adjust:
                continue
            lb, ub = b
            if lb is None and ub is None:
                n_unbounded += 1
                x[i] = x[i] - x[n_x + n_unbounded - 1]
            else:
                if lb is None:
                    x[i] = ub - x[i]
                else:
                    x[i] += lb

    n_x = len(c)
    x = x[:n_x]  # all the rest of the variables were artificial
    fun = x.dot(c)
    slack = b_ub - A_ub.dot(x)  # report slack for ORIGINAL UB constraints
    # report residuals of ORIGINAL EQ constraints
    con = b_eq - A_eq.dot(x)

    # Patch for bug #8664. Detecting this sort of issue earlier
    # (via abnormalities in the indicators) would be better.
    bounds = np.array(bounds)  # again, this should have been the standard form
    lb = bounds[:, 0]
    ub = bounds[:, 1]
    lb[np.equal(lb, None)] = -np.inf
    ub[np.equal(ub, None)] = np.inf

    return x, fun, slack, con, lb, ub


def _check_result(x, fun, status, slack, con, lb, ub, tol, message):
    """
    Check the validity of the provided solution.

    A valid (optimal) solution satisfies all bounds, all slack variables are
    negative and all equality constraint residuals are strictly non-zero.
    Further, the lower-bounds, upper-bounds, slack and residuals contain
    no nan values.

    Parameters
    ----------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties encountered

    slack : 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    lb : 1D array
        The lower bound constraints on the original variables
    ub: 1D array
        The upper bound constraints on the original variables
    message : str
        A string descriptor of the exit status of the optimization.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.
    """
    # Somewhat arbitrary, but status 5 is very unusual
    tol = np.sqrt(tol) * 10

    contains_nans = (
        np.isnan(x).any()
        or np.isnan(fun)
        or np.isnan(slack).any()
        or np.isnan(con).any()
    )

    if contains_nans:
        is_feasible = False
    else:
        invalid_bounds = (x < lb - tol).any() or (x > ub + tol).any()
        invalid_slack = status != 3 and (slack < -tol).any()
        invalid_con = status != 3 and (np.abs(con) > tol).any()
        is_feasible = not (invalid_bounds or invalid_slack or invalid_con)

    if status == 0 and not is_feasible:
        status = 4
        message = ("The solution does not satisfy the constraints within the "
                   "required tolerance of " + "{:.2E}".format(tol) + ", yet "
                   "no errors were raised and there is no certificate of "
                   "infeasibility or unboundedness. This is known to occur "
                   "if the `presolve` option is False and the problem is "
                   "infeasible. This can also occur due to the limited "
                   "accuracy of the `interior-point` method. Check whether "
                   "the slack and constraint residuals are acceptable; "
                   "if not, consider enabling presolve, reducing option "
                   "`tol`, and/or using method `revised simplex`. "
                   "If you encounter this message under different "
                   "circumstances, please submit a bug report.")
    elif status == 0 and contains_nans:
        status = 4
        message = ("Numerical difficulties were encountered but no errors "
                   "were raised. This is known to occur if the 'presolve' "
                   "option is False, 'sparse' is True, and A_eq includes "
                   "redundant rows. If you encounter this under different "
                   "circumstances, please submit a bug report. Otherwise, "
                   "remove linearly dependent equations from your equality "
                   "constraints or enable presolve.")
    elif status == 2 and is_feasible:
        # Occurs if the simplex method exits after phase one with a very
        # nearly basic feasible solution. Postsolving can make the solution
        # basic, however, this solution is NOT optimal
        raise ValueError(message)

    return status, message


def _postprocess(x, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
                 complete=False, undo=[], status=0, message="", tol=1e-8,
                 iteration=None, disp=False):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1D array
        Solution vector to the standard-form problem.
    c : 1D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2D array, optional
        2D array such that ``A_ub @ x`` gives the values of the upper-bound
        inequality constraints at ``x``.
    b_ub : 1D array, optional
        1D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2D array, optional
        2D array such that ``A_eq @ x`` gives the values of the equality
        constraints at ``x``.
    b_eq : 1D array, optional
        1D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        Bounds, as modified in presolve
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.

    Returns
    -------
    x : 1D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack : 1D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    status : int
        An integer representing the exit status of the optimization::

             0 : Optimization terminated successfully
             1 : Iteration limit reached
             2 : Problem appears to be infeasible
             3 : Problem appears to be unbounded
             4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.

    """

    x, fun, slack, con, lb, ub = _postsolve(
        x, c, A_ub, b_ub, A_eq, b_eq,
        bounds, complete, undo, tol
    )

    status, message = _check_result(
        x, fun, status, slack, con,
        lb, ub, tol, message
    )

    if disp:
        _display_summary(message, status, fun, iteration)

    return x, fun, slack, con, status, message
