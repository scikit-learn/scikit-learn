"""
An interior-point method for linear programming.
"""
# Author: Matt Haberland

from __future__ import print_function, division, absolute_import
import numpy as np
import scipy as sp
import scipy.sparse as sps
from warnings import warn
from scipy.linalg import LinAlgError
from .optimize import OptimizeResult, OptimizeWarning, _check_unknown_options
from scipy.optimize._remove_redundancy import _remove_redundancy
from scipy.optimize._remove_redundancy import _remove_redundancy_sparse
from scipy.optimize._remove_redundancy import _remove_redundancy_dense


def _clean_inputs(
        c,
        A_ub=None,
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=None):
    """
    Given user inputs for a linear programming problem, return the
    objective vector, upper bound constraints, equality constraints,
    and simple bounds in a preferred format.

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

    Returns
    -------
    c : 1-D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``.
    b_ub : 1-D array
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``.
    b_eq : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative)

    """

    try:
        if c is None:
            raise TypeError
        try:
            c = np.asarray(c, dtype=float).copy().squeeze()
        except BaseException:  # typically a ValueError and shouldn't be, IMO
            raise TypeError
        if c.size == 1:
            c = c.reshape((-1))
        n_x = len(c)
        if n_x == 0 or len(c.shape) != 1:
            raise ValueError(
                "Invalid input for linprog: c should be a 1D array; it must "
                "not have more than one non-singleton dimension")
        if not(np.isfinite(c).all()):
            raise ValueError(
                "Invalid input for linprog: c must not contain values "
                "inf, nan, or None")
    except TypeError:
        raise TypeError(
            "Invalid input for linprog: c must be a 1D array of numerical "
            "coefficients")

    try:
        try:
            if sps.issparse(A_eq) or sps.issparse(A_ub):
                A_ub = sps.coo_matrix(
                    (0, n_x), dtype=float) if A_ub is None else sps.coo_matrix(
                    A_ub, dtype=float).copy()
            else:
                A_ub = np.zeros(
                    (0, n_x), dtype=float) if A_ub is None else np.asarray(
                    A_ub, dtype=float).copy()
        except BaseException:
            raise TypeError
        n_ub = A_ub.shape[0]
        if len(A_ub.shape) != 2 or A_ub.shape[1] != len(c):
            raise ValueError(
                "Invalid input for linprog: A_ub must have exactly two "
                "dimensions, and the number of columns in A_ub must be "
                "equal to the size of c ")
        if (sps.issparse(A_ub) and not np.isfinite(A_ub.data).all()
                or not sps.issparse(A_ub) and not np.isfinite(A_ub).all()):
            raise ValueError(
                "Invalid input for linprog: A_ub must not contain values "
                "inf, nan, or None")
    except TypeError:
        raise TypeError(
            "Invalid input for linprog: A_ub must be a numerical 2D array "
            "with each row representing an upper bound inequality constraint")

    try:
        try:
            b_ub = np.array(
                [], dtype=float) if b_ub is None else np.asarray(
                b_ub, dtype=float).copy().squeeze()
        except BaseException:
            raise TypeError
        if b_ub.size == 1:
            b_ub = b_ub.reshape((-1))
        if len(b_ub.shape) != 1:
            raise ValueError(
                "Invalid input for linprog: b_ub should be a 1D array; it "
                "must not have more than one non-singleton dimension")
        if len(b_ub) != n_ub:
            raise ValueError(
                "Invalid input for linprog: The number of rows in A_ub must "
                "be equal to the number of values in b_ub")
        if not(np.isfinite(b_ub).all()):
            raise ValueError(
                "Invalid input for linprog: b_ub must not contain values "
                "inf, nan, or None")
    except TypeError:
        raise TypeError(
            "Invalid input for linprog: b_ub must be a 1D array of "
            "numerical values, each representing the upper bound of an "
            "inequality constraint (row) in A_ub")

    try:
        try:
            if sps.issparse(A_eq) or sps.issparse(A_ub):
                A_eq = sps.coo_matrix(
                    (0, n_x), dtype=float) if A_eq is None else sps.coo_matrix(
                    A_eq, dtype=float).copy()
            else:
                A_eq = np.zeros(
                    (0, n_x), dtype=float) if A_eq is None else np.asarray(
                    A_eq, dtype=float).copy()
        except BaseException:
            raise TypeError
        n_eq = A_eq.shape[0]
        if len(A_eq.shape) != 2 or A_eq.shape[1] != len(c):
            raise ValueError(
                "Invalid input for linprog: A_eq must have exactly two "
                "dimensions, and the number of columns in A_eq must be "
                "equal to the size of c ")

        if (sps.issparse(A_eq) and not np.isfinite(A_eq.data).all()
                or not sps.issparse(A_eq) and not np.isfinite(A_eq).all()):
            raise ValueError(
                "Invalid input for linprog: A_eq must not contain values "
                "inf, nan, or None")
    except TypeError:
        raise TypeError(
            "Invalid input for linprog: A_eq must be a 2D array with each "
            "row representing an equality constraint")

    try:
        try:
            b_eq = np.array(
                [], dtype=float) if b_eq is None else np.asarray(
                b_eq, dtype=float).copy().squeeze()
        except BaseException:
            raise TypeError
        if b_eq.size == 1:
            b_eq = b_eq.reshape((-1))
        if len(b_eq.shape) != 1:
            raise ValueError(
                "Invalid input for linprog: b_eq should be a 1D array; it "
                "must not have more than one non-singleton dimension")
        if len(b_eq) != n_eq:
            raise ValueError(
                "Invalid input for linprog: the number of rows in A_eq "
                "must be equal to the number of values in b_eq")
        if not(np.isfinite(b_eq).all()):
            raise ValueError(
                "Invalid input for linprog: b_eq must not contain values "
                "inf, nan, or None")
    except TypeError:
        raise TypeError(
            "Invalid input for linprog: b_eq must be a 1D array of "
            "numerical values, each representing the right hand side of an "
            "equality constraints (row) in A_eq")

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

    return c, A_ub, b_ub, A_eq, b_eq, bounds


def _presolve(c, A_ub, b_ub, A_eq, b_eq, bounds, rr):
    """
    Given inputs for a linear programming problem in preferred format,
    presolve the problem: identify trivial infeasibilities, redundancies,
    and unboundedness, tighten bounds where possible, and eliminate fixed
    variables.

    Parameters
    ----------
    c : 1-D array
        Coefficients of the linear objective function to be minimized.
    A_ub : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``.
    b_ub : 1-D array
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``.
    A_eq : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``.
    b_eq : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction.

    Returns
    -------
    c : 1-D array
        Coefficients of the linear objective function to be minimized.
    c0 : 1-D array
        Constant term in objective function due to fixed (and eliminated)
        variables.
    A_ub : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``. Unnecessary
        rows/columns have been removed.
    b_ub : 1-D array
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``. Unnecessary elements have been removed.
    A_eq : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``. Unnecessary rows/columns have been
        removed.
    b_eq : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``. Unnecessary elements have been removed.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. Bounds have been
        tightened where possible.
    x : 1-D array
        Solution vector (when the solution is trivial and can be determined
        in presolve)
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

    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [2] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [5] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
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

    tol = 1e-9    # tolerance for equality. should this be exposed to user?

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
                    x, undo, complete, status, message)
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
                    x, undo, complete, status, message)
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
                    x, undo, complete, status, message)
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
                        x, undo, complete, status, message)
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
                        x, undo, complete, status, message)
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
                    x, undo, complete, status, message)

    ub_mod = ub
    lb_mod = lb
    if np.any(i_f):
        c0 += c[i_f].dot(lb[i_f])
        b_eq = b_eq - A_eq[:, i_f].dot(lb[i_f])
        b_ub = b_ub - A_ub[:, i_f].dot(lb[i_f])
        c = c[i_nf]
        x = x[i_nf]
        A_eq = A_eq[:, i_nf]
        A_ub = A_ub[:, i_nf]
        # record of variables to be added back in
        undo = [np.where(i_f)[0], lb[i_f]]
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
        elif (np.any(np.logical_and(c < 0, ub == np.inf)) or
                np.any(np.logical_and(c > 0, lb == -np.inf))):
                # test_no_constraints()
            status = 3
            message = ("If feasible, the problem is (trivially) unbounded "
                       "because there are no constraints and at least one "
                       "element of c is negative. If you wish to check "
                       "whether the problem is infeasible, turn presolve "
                       "off.")
        else:  # test_empty_constraint_2
            status = 0
            message = ("The solution was determined in presolve as there are "
                       "no non-trivial constraints.")
        complete = True
        x[c < 0] = ub_mod[c < 0]
        x[c > 0] = lb_mod[c > 0]
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
            if str(
                    col) == "nan":  # comparing col to float("nan") and
                                    # np.nan doesn't work. should use np.isnan
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
                x, undo, complete, status, message)

    # This is a wild guess for which redundancy removal algorithm will be
    # faster. More testing would be good.
    small_nullspace = 5
    if rr and A_eq.size > 0:
        try:  # TODO: instead use results of first SVD in _remove_redundancy
            rank = np.linalg.matrix_rank(A_eq)
        except:  # oh well, we'll have to go with _remove_redundancy_dense
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
            x, undo, complete, status, message)


def _get_Abc(
        c,
        c0=0,
        A_ub=None,
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=None,
        undo=[]):
    """
    Given a linear programming problem of the form:

    minimize:     c^T * x

    subject to:   A_ub * x <= b_ub
                  A_eq * x == b_eq
                  bounds[i][0] < x_i < bounds[i][1]

    return the problem in standard form:
    minimize:     c'^T * x'

    subject to:   A * x' == b
                  0 < x' < oo

    by adding slack variables and making variable substitutions as necessary.

    Parameters
    ----------
    c : 1-D array
        Coefficients of the linear objective function to be minimized.
        Components corresponding with fixed variables have been eliminated.
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.
    A_ub : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the upper-bound inequality constraints at ``x``. Unnecessary
        rows/columns have been removed.
    b_ub : 1-D array
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in ``A_ub``. Unnecessary elements have been removed.
    A_eq : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x``. Unnecessary rows/columns have been
        removed.
    b_eq : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A_eq``. Unnecessary elements have been removed.
    bounds : sequence of tuples
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for each of ``min`` or
        ``max`` when there is no bound in that direction. Bounds have been
        tightened where possible.
    undo: list of tuples
        (`index`, `value`) pairs that record the original index and fixed value
        for each variable removed from the problem

    Returns
    -------
    A : 2-D array
        2-D array which, when matrix-multiplied by x, gives the values of the
        equality constraints at x (for standard form problem).
    b : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in A (for standard form problem).
    c : 1-D array
        Coefficients of the linear objective function to be minimized (for
        standard form problem).
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables.

    References
    ----------
    .. [6] Bertsimas, Dimitris, and J. Tsitsiklis. "Introduction to linear
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
    i_nolb = np.where(l_nolb_someub)[0]
    lbs[l_nolb_someub], ubs[l_nolb_someub] = (
        -ubs[l_nolb_someub], lbs[l_nolb_someub])
    lb_none = np.equal(lbs, None)
    ub_none = np.equal(ubs, None)
    lb_some = np.logical_not(lb_none)
    ub_some = np.logical_not(ub_none)
    c[i_nolb] *= -1
    if len(i_nolb) > 0:
        if A_ub.shape[0] > 0:  # sometimes needed for sparse arrays... weird
            A_ub[:, i_nolb] *= -1
        if A_eq.shape[0] > 0:
            A_eq[:, i_nolb] *= -1

    # upper bound: add inequality constraint
    i_newub = np.where(ub_some)[0]
    ub_newub = ubs[ub_some]
    n_bounds = np.count_nonzero(ub_some)
    A_ub = vstack((A_ub, zeros((n_bounds, A_ub.shape[1]))))
    b_ub = np.concatenate((b_ub, np.zeros(n_bounds)))
    A_ub[range(m_ub, A_ub.shape[0]), i_newub] = 1
    b_ub[m_ub:] = ub_newub

    A1 = vstack((A_ub, A_eq))
    b = np.concatenate((b_ub, b_eq))
    c = np.concatenate((c, np.zeros((A_ub.shape[0],))))

    # unbounded: substitute xi = xi+ + xi-
    l_free = np.logical_and(lb_none, ub_none)
    i_free = np.where(l_free)[0]
    n_free = len(i_free)
    A1 = hstack((A1, zeros((A1.shape[0], n_free))))
    c = np.concatenate((c, np.zeros(n_free)))
    A1[:, range(n_ub, A1.shape[1])] = -A1[:, i_free]
    c[np.arange(n_ub, A1.shape[1])] = -c[i_free]

    # add slack variables
    A2 = vstack([eye(A_ub.shape[0]), zeros((A_eq.shape[0], A_ub.shape[0]))])
    A = hstack([A1, A2])

    # lower bound: substitute xi = xi' + lb
    # now there is a constant term in objective
    i_shift = np.where(lb_some)[0]
    lb_shift = lbs[lb_some].astype(float)
    c0 += np.sum(lb_shift * c[i_shift])
    if sparse:
        b = b.reshape(-1, 1)
        A = A.tocsc()
        b -= (A[:, i_shift] * sps.diags(lb_shift)).sum(axis=1)
        b = b.ravel()
    else:
        b -= (A[:, i_shift] * lb_shift).sum(axis=1)

    return A, b, c, c0


def _postprocess(
        x,
        c,
        A_ub=None,
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=None,
        complete=False,
        undo=[],
        status=0,
        message="",
        tol=1e-8):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1-D array
        Solution vector to the standard-form problem.
    c : 1-D array
        Original coefficients of the linear objective function to be minimized.
    A_ub : 2-D array
        Original upper bound constraint matrix.
    b_ub : 1-D array
        Original upper bound constraint vector.
    A_eq : 2-D array
        Original equality constraint matrix.
    b_eq : 1-D array
        Original equality constraint vector.
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
    x : 1-D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack: 1-D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub * x``
    con : 1-D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq * x``
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
        x = np.array(x)

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
    tol = np.sqrt(tol)  # Somewhat arbitrary, but status 5 is very unusual
    if status == 0 and ((slack < -tol).any() or (np.abs(con) > tol).any() or
                        (x < lb - tol).any() or (x > ub + tol).any()):
        status = 4
        message = ("The solution does not satisfy the constraints, yet "
                   "no errors were raised and there is no certificate of "
                   "infeasibility or unboundedness. This is known to occur "
                   "if the `presolve` option is False and the problem is "
                   "infeasible. If you uncounter this under different "
                   "circumstances, please submit a bug report. Otherwise, "
                   "please enable presolve.")
    elif status == 0 and (np.isnan(x).any() or np.isnan(fun) or
                          np.isnan(slack).any() or np.isnan(con).any()):
        status = 4
        message = ("Numerical difficulties were encountered but no errors "
                   "were raised. This is known to occur if the 'presolve' "
                   "option is False, 'sparse' is True, and A_eq includes "
                   "redundant rows. If you encounter this under different "
                   "circumstances, please submit a bug report. Otherwise, "
                   "remove linearly dependent equations from your equality "
                   "constraints or enable presolve.")

    return x, fun, slack, con, status, message


def _get_solver(sparse=False, lstsq=False, sym_pos=True, cholesky=True):
    """
    Given solver options, return a handle to the appropriate linear system
    solver.

    Parameters
    ----------
    sparse : bool
        True if the system to be solved is sparse. This is typically set
        True when the original ``A_ub`` and ``A_eq`` arrays are sparse.
    lstsq : bool
        True if the system is ill-conditioned and/or (nearly) singular and
        thus a more robust least-squares solver is desired. This is sometimes
        needed as the solution is approached.
    sym_pos : bool
        True if the system matrix is symmetric positive definite
        Sometimes this needs to be set false as the solution is approached,
        even when the system should be symmetric positive definite, due to
        numerical difficulties.
    cholesky : bool
        True if the system is to be solved by Cholesky, rather than LU,
        decomposition. This is typically faster unless the problem is very
        small or prone to numerical difficulties.

    Returns
    -------
    solve : function
        Handle to the appropriate solver function

    """
    if sparse:
        if lstsq or not(sym_pos):
            def solve(M, r, sym_pos=False):
                return sps.linalg.lsqr(M, r)[0]
        else:
            # this is not currently used; it is replaced by splu solve
            # TODO: expose use of this as an option
            def solve(M, r):
                return sps.linalg.spsolve(M, r, permc_spec="MMD_AT_PLUS_A")

    else:
        if lstsq:  # sometimes necessary as solution is approached
            def solve(M, r):
                return sp.linalg.lstsq(M, r)[0]
        elif cholesky:
            solve = sp.linalg.cho_solve
        else:
            # this seems to cache the matrix factorization, so solving
            # with multiple right hand sides is much faster
            def solve(M, r, sym_pos=sym_pos):
                return sp.linalg.solve(M, r, sym_pos=sym_pos)

    return solve


def _get_delta(
    A,
    b,
    c,
    x,
    y,
    z,
    tau,
    kappa,
    gamma,
    eta,
    sparse=False,
    lstsq=False,
    sym_pos=True,
    cholesky=True,
    pc=True,
    ip=False,
        permc_spec='MMD_AT_PLUS_A'):
    """
    Given standard form problem defined by ``A``, ``b``, and ``c``;
    current variable estimates ``x``, ``y``, ``z``, ``tau``, and ``kappa``;
    algorithmic parameters ``gamma and ``eta;
    and options ``sparse``, ``lstsq``, ``sym_pos``, ``cholesky``, ``pc``
    (predictor-corrector), and ``ip`` (initial point improvement),
    get the search direction for increments to the variable estimates.

    Parameters
    ----------
    As defined in [1], except:
    sparse : bool
        True if the system to be solved is sparse. This is typically set
        True when the original ``A_ub`` and ``A_eq`` arrays are sparse.
    lstsq : bool
        True if the system is ill-conditioned and/or (nearly) singular and
        thus a more robust least-squares solver is desired. This is sometimes
        needed as the solution is approached.
    sym_pos : bool
        True if the system matrix is symmetric positive definite
        Sometimes this needs to be set false as the solution is approached,
        even when the system should be symmetric positive definite, due to
        numerical difficulties.
    cholesky : bool
        True if the system is to be solved by Cholesky, rather than LU,
        decomposition. This is typically faster unless the problem is very
        small or prone to numerical difficulties.
    pc : bool
        True if the predictor-corrector method of Mehrota is to be used. This
        is almost always (if not always) beneficial. Even though it requires
        the solution of an additional linear system, the factorization
        is typically (implicitly) reused so solution is efficient, and the
        number of algorithm iterations is typically reduced.
    ip : bool
        True if the improved initial point suggestion due to [1] section 4.3
        is desired. It's unclear whether this is beneficial.
    permc_spec : str (default = 'MMD_AT_PLUS_A')
        (Has effect only with ``sparse = True``, ``lstsq = False``, ``sym_pos =
        True``.) A matrix is factorized in each iteration of the algorithm.
        This option specifies how to permute the columns of the matrix for
        sparsity preservation. Acceptable values are:

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering.

        This option can impact the convergence of the
        interior point algorithm; test different values to determine which
        performs best for your problem. For more information, refer to
        ``scipy.sparse.linalg.splu``.

    Returns
    -------
    Search directions as defined in [1]

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """

    if A.shape[0] == 0:
        # If there are no constraints, some solvers fail (understandably)
        # rather than returning empty solution. This gets the job done.
        sparse, lstsq, sym_pos, cholesky = False, False, True, False
    solve = _get_solver(sparse, lstsq, sym_pos, cholesky)
    n_x = len(x)

    # [1] Equation 8.8
    r_P = b * tau - A.dot(x)
    r_D = c * tau - A.T.dot(y) - z
    r_G = c.dot(x) - b.transpose().dot(y) + kappa
    mu = (x.dot(z) + tau * kappa) / (n_x + 1)

    #  Assemble M from [1] Equation 8.31
    Dinv = x / z
    splu = False
    if sparse and not lstsq:
        # sparse requires Dinv to be diag matrix
        M = A.dot(sps.diags(Dinv, 0, format="csc").dot(A.T))
        try:
            # TODO: should use linalg.factorized instead, but I don't have
            #       umfpack and therefore cannot test its performance
            solve = sps.linalg.splu(M, permc_spec=permc_spec).solve
            splu = True
        except:
            lstsq = True
            solve = _get_solver(sparse, lstsq, sym_pos, cholesky)
    else:
        # dense does not; use broadcasting
        M = A.dot(Dinv.reshape(-1, 1) * A.T)

    # For some small problems, calling sp.linalg.solve w/ sym_pos = True
    # may be faster. I am pretty certain it caches the factorization for
    # multiple uses and checks the incoming matrix to see if it's the same as
    # the one it already factorized. (I can't explain the speed otherwise.)
    if cholesky:
        try:
            L = sp.linalg.cho_factor(M)
        except:
            cholesky = False
            solve = _get_solver(sparse, lstsq, sym_pos, cholesky)

    # pc: "predictor-corrector" [1] Section 4.1
    # In development this option could be turned off
    # but it always seems to improve performance substantially
    n_corrections = 1 if pc else 0

    i = 0
    alpha, d_x, d_z, d_tau, d_kappa = 0, 0, 0, 0, 0
    while i <= n_corrections:
        # Reference [1] Eq. 8.6
        rhatp = eta(gamma) * r_P
        rhatd = eta(gamma) * r_D
        rhatg = np.array(eta(gamma) * r_G).reshape((1,))

        # Reference [1] Eq. 8.7
        rhatxs = gamma * mu - x * z
        rhattk = np.array(gamma * mu - tau * kappa).reshape((1,))

        if i == 1:
            if ip:  # if the correction is to get "initial point"
                # Reference [1] Eq. 8.23
                rhatxs = ((1 - alpha) * gamma * mu -
                          x * z - alpha**2 * d_x * d_z)
                rhattk = np.array(
                    (1 -
                     alpha) *
                    gamma *
                    mu -
                    tau *
                    kappa -
                    alpha**2 *
                    d_tau *
                    d_kappa).reshape(
                    (1,
                     ))
            else:  # if the correction is for "predictor-corrector"
                # Reference [1] Eq. 8.13
                rhatxs -= d_x * d_z
                rhattk -= d_tau * d_kappa

        # sometimes numerical difficulties arise as the solution is approached
        # this loop tries to solve the equations using a sequence of functions
        # for solve. For dense systems, the order is:
        # 1. scipy.linalg.cho_factor/scipy.linalg.cho_solve,
        # 2. scipy.linalg.solve w/ sym_pos = True,
        # 3. scipy.linalg.solve w/ sym_pos = False, and if all else fails
        # 4. scipy.linalg.lstsq
        # For sparse systems, the order is:
        # 1. scipy.sparse.linalg.splu
        # 2. scipy.sparse.linalg.lsqr
        # TODO: if umfpack is installed, use factorized instead of splu.
        #       Can't do that now because factorized doesn't pass permc_spec
        #       to splu if umfpack isn't installed. Also, umfpack not tested.
        solved = False
        while(not solved):
            try:
                solve_this = L if cholesky else M
                # [1] Equation 8.28
                p, q = _sym_solve(Dinv, solve_this, A, c, b, solve, splu)
                # [1] Equation 8.29
                u, v = _sym_solve(Dinv, solve_this, A, rhatd -
                                  (1 / x) * rhatxs, rhatp, solve, splu)
                if np.any(np.isnan(p)) or np.any(np.isnan(q)):
                    raise LinAlgError
                solved = True
            except (LinAlgError, ValueError) as e:
                # Usually this doesn't happen. If it does, it happens when
                # there are redundant constraints or when approaching the
                # solution. If so, change solver.
                cholesky = False
                if not lstsq:
                    if sym_pos:
                        warn(
                            "Solving system with option 'sym_pos':True "
                            "failed. It is normal for this to happen "
                            "occasionally, especially as the solution is "
                            "approached. However, if you see this frequently, "
                            "consider setting option 'sym_pos' to False.",
                            OptimizeWarning)
                        sym_pos = False
                    else:
                        warn(
                            "Solving system with option 'sym_pos':False "
                            "failed. This may happen occasionally, "
                            "especially as the solution is "
                            "approached. However, if you see this frequently, "
                            "your problem may be numerically challenging. "
                            "If you cannot improve the formulation, consider "
                            "setting 'lstsq' to True. Consider also setting "
                            "`presolve` to True, if it is not already.",
                            OptimizeWarning)
                        lstsq = True
                else:
                    raise e
                solve = _get_solver(sparse, lstsq, sym_pos)
        # [1] Results after 8.29
        d_tau = ((rhatg + 1 / tau * rhattk - (-c.dot(u) + b.dot(v))) /
                 (1 / tau * kappa + (-c.dot(p) + b.dot(q))))
        d_x = u + p * d_tau
        d_y = v + q * d_tau

        # [1] Relations between  after 8.25 and 8.26
        d_z = (1 / x) * (rhatxs - z * d_x)
        d_kappa = 1 / tau * (rhattk - kappa * d_tau)

        # [1] 8.12 and "Let alpha be the maximal possible step..." before 8.23
        alpha = _get_step(x, d_x, z, d_z, tau, d_tau, kappa, d_kappa, 1)
        if ip:  # initial point - see [1] 4.4
            gamma = 10
        else:  # predictor-corrector, [1] definition after 8.12
            beta1 = 0.1  # [1] pg. 220 (Table 8.1)
            gamma = (1 - alpha)**2 * min(beta1, (1 - alpha))
        i += 1

    return d_x, d_y, d_z, d_tau, d_kappa


def _sym_solve(Dinv, M, A, r1, r2, solve, splu=False):
    """
    An implementation of [1] equation 8.31 and 8.32

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    # [1] 8.31
    r = r2 + A.dot(Dinv * r1)
    if splu:
        v = solve(r)
    else:
        v = solve(M, r)
    # [1] 8.32
    u = Dinv * (A.T.dot(v) - r1)
    return u, v


def _get_step(x, d_x, z, d_z, tau, d_tau, kappa, d_kappa, alpha0):
    """
    An implementation of [1] equation 8.21

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    # [1] 4.3 Equation 8.21, ignoring 8.20 requirement
    # same step is taken in primal and dual spaces
    # alpha0 is basically beta3 from [1] Table 8.1, but instead of beta3
    # the value 1 is used in Mehrota corrector and initial point correction
    i_x = d_x < 0
    i_z = d_z < 0
    alpha_x = alpha0 * np.min(x[i_x] / -d_x[i_x]) if np.any(i_x) else 1
    alpha_tau = alpha0 * tau / -d_tau if d_tau < 0 else 1
    alpha_z = alpha0 * np.min(z[i_z] / -d_z[i_z]) if np.any(i_z) else 1
    alpha_kappa = alpha0 * kappa / -d_kappa if d_kappa < 0 else 1
    alpha = np.min([1, alpha_x, alpha_tau, alpha_z, alpha_kappa])
    return alpha


def _get_message(status):
    """
    Given problem status code, return a more detailed message.

    Parameters
    ----------
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered.

    Returns
    -------
    message : str
        A string descriptor of the exit status of the optimization.

    """
    messages = (
        ["Optimization terminated successfully.",
         "The iteration limit was reached before the algorithm converged.",
         "The algorithm terminated successfully and determined that the "
         "problem is infeasible.",
         "The algorithm terminated successfully and determined that the "
         "problem is unbounded.",
         "Numerical difficulties were encountered before the problem "
         "converged. Please check your problem formulation for errors, "
         "independence of linear equality constraints, and reasonable "
         "scaling and matrix condition numbers. If you continue to "
         "encounter this error, please submit a bug report."
         ])
    return messages[status]


def _do_step(x, y, z, tau, kappa, d_x, d_y, d_z, d_tau, d_kappa, alpha):
    """
    An implementation of [1] Equation 8.9

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    x = x + alpha * d_x
    tau = tau + alpha * d_tau
    z = z + alpha * d_z
    kappa = kappa + alpha * d_kappa
    y = y + alpha * d_y
    return x, y, z, tau, kappa


def _get_blind_start(shape):
    """
    Return the starting point from [1] 4.4

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    m, n = shape
    x0 = np.ones(n)
    y0 = np.zeros(m)
    z0 = np.ones(n)
    tau0 = 1
    kappa0 = 1
    return x0, y0, z0, tau0, kappa0


def _indicators(A, b, c, c0, x, y, z, tau, kappa):
    """
    Implementation of several equations from [1] used as indicators of
    the status of optimization.

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """

    # residuals for termination are relative to initial values
    x0, y0, z0, tau0, kappa0 = _get_blind_start(A.shape)

    # See [1], Section 4 - The Homogeneous Algorithm, Equation 8.8
    def r_p(x, tau):
        return b * tau - A.dot(x)

    def r_d(y, z, tau):
        return c * tau - A.T.dot(y) - z

    def r_g(x, y, kappa):
        return kappa + c.dot(x) - b.dot(y)

    # np.dot unpacks if they are arrays of size one
    def mu(x, tau, z, kappa):
        return (x.dot(z) + np.dot(tau, kappa)) / (len(x) + 1)

    obj = c.dot(x / tau) + c0

    def norm(a):
        return np.linalg.norm(a)

    # See [1], Section 4.5 - The Stopping Criteria
    r_p0 = r_p(x0, tau0)
    r_d0 = r_d(y0, z0, tau0)
    r_g0 = r_g(x0, y0, kappa0)
    mu_0 = mu(x0, tau0, z0, kappa0)
    rho_A = norm(c.T.dot(x) - b.T.dot(y)) / (tau + norm(b.T.dot(y)))
    rho_p = norm(r_p(x, tau)) / max(1, norm(r_p0))
    rho_d = norm(r_d(y, z, tau)) / max(1, norm(r_d0))
    rho_g = norm(r_g(x, y, kappa)) / max(1, norm(r_g0))
    rho_mu = mu(x, tau, z, kappa) / mu_0
    return rho_p, rho_d, rho_A, rho_g, rho_mu, obj


def _display_iter(rho_p, rho_d, rho_g, alpha, rho_mu, obj, header=False):
    """
    Print indicators of optimization status to the console.

    Parameters
    ----------
    rho_p : float
        The (normalized) primal feasibility, see [1] 4.5
    rho_d : float
        The (normalized) dual feasibility, see [1] 4.5
    rho_g : float
        The (normalized) duality gap, see [1] 4.5
    alpha : float
        The step size, see [1] 4.3
    rho_mu : float
        The (normalized) path parameter, see [1] 4.5
    obj : float
        The objective function value of the current iterate
    header : bool
        True if a header is to be printed

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    if header:
        print("Primal Feasibility ",
              "Dual Feasibility   ",
              "Duality Gap        ",
              "Step            ",
              "Path Parameter     ",
              "Objective          ")

    # no clue why this works
    fmt = '{0:<20.13}{1:<20.13}{2:<20.13}{3:<17.13}{4:<20.13}{5:<20.13}'
    print(fmt.format(
        rho_p,
        rho_d,
        rho_g,
        alpha,
        rho_mu,
        obj))


def _ip_hsd(A, b, c, c0, alpha0, beta, maxiter, disp, tol,
            sparse, lstsq, sym_pos, cholesky, pc, ip, permc_spec):
    r"""
    Solve a linear programming problem in standard form:

    minimize:     c'^T * x'

    subject to:   A * x' == b
                  0 < x' < oo

    using the interior point method of [1].

    Parameters
    ----------
    A : 2-D array
        2-D array which, when matrix-multiplied by ``x``, gives the values of
        the equality constraints at ``x`` (for standard form problem).
    b : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A`` (for standard form problem).
    c : 1-D array
        Coefficients of the linear objective function to be minimized (for
        standard form problem).
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables. (Purely for display.)
    alpha0 : float
        The maximal step size for Mehrota's predictor-corrector search
        direction; see :math:`\beta_3`of [1] Table 8.1
    beta : float
        The desired reduction of the path parameter :math:`\mu` (see  [3]_)
    maxiter : int
        The maximum number of iterations of the algorithm.
    disp : bool
        Set to ``True`` if indicators of optimization status are to be printed
        to the console each iteration.
    tol : float
        Termination tolerance; see [1]_ Section 4.5.
    sparse : bool
        Set to ``True`` if the problem is to be treated as sparse. However,
        the inputs ``A_eq`` and ``A_ub`` should nonetheless be provided as
        (dense) arrays rather than sparse matrices.
    lstsq : bool
        Set to ``True`` if the problem is expected to be very poorly
        conditioned. This should always be left as ``False`` unless severe
        numerical difficulties are frequently encountered, and a better option
        would be to improve the formulation of the problem.
    sym_pos : bool
        Leave ``True`` if the problem is expected to yield a well conditioned
        symmetric positive definite normal equation matrix (almost always).
    cholesky : bool
        Set to ``True`` if the normal equations are to be solved by explicit
        Cholesky decomposition followed by explicit forward/backward
        substitution. This is typically faster for moderate, dense problems
        that are numerically well-behaved.
    pc : bool
        Leave ``True`` if the predictor-corrector method of Mehrota is to be
        used. This is almost always (if not always) beneficial.
    ip : bool
        Set to ``True`` if the improved initial point suggestion due to [1]_
        Section 4.3 is desired. It's unclear whether this is beneficial.
    permc_spec : str (default = 'MMD_AT_PLUS_A')
        (Has effect only with ``sparse = True``, ``lstsq = False``, ``sym_pos =
        True``.) A matrix is factorized in each iteration of the algorithm.
        This option specifies how to permute the columns of the matrix for
        sparsity preservation. Acceptable values are:

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering.

        This option can impact the convergence of the
        interior point algorithm; test different values to determine which
        performs best for your problem. For more information, refer to
        ``scipy.sparse.linalg.splu``.

    Returns
    -------
    x_hat : float
        Solution vector (for standard form problem).
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered.

    message : str
        A string descriptor of the exit status of the optimization.
    iteration : int
        The number of iterations taken to solve the problem

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.
    .. [3] Freund, Robert M. "Primal-Dual Interior-Point Methods for Linear
           Programming based on Newton's Method." Unpublished Course Notes,
           March 2004. Available 2/25/2017 at:
           https://ocw.mit.edu/courses/sloan-school-of-management/15-084j-nonlinear-programming-spring-2004/lecture-notes/lec14_int_pt_mthd.pdf

    """

    iteration = 0

    # default initial point
    x, y, z, tau, kappa = _get_blind_start(A.shape)

    # first iteration is special improvement of initial point
    ip = ip if pc else False

    # [1] 4.5
    rho_p, rho_d, rho_A, rho_g, rho_mu, obj = _indicators(
        A, b, c, c0, x, y, z, tau, kappa)
    go = rho_p > tol or rho_d > tol or rho_A > tol  # we might get lucky : )

    if disp:
        _display_iter(rho_p, rho_d, rho_g, "-", rho_mu, obj, header=True)

    status = 0
    message = "Optimization terminated successfully."

    if sparse:
        A = sps.csc_matrix(A)
        A.T = A.transpose()  # A.T is defined for sparse matrices but is slow
        # Redefine it to avoid calculating again
        # This is fine as long as A doesn't change

    while go:

        iteration += 1

        if ip:  # initial point
            # [1] Section 4.4
            gamma = 1

            def eta(g):
                return 1
        else:
            # gamma = 0 in predictor step according to [1] 4.1
            # if predictor/corrector is off, use mean of complementarity [3]
            # 5.1 / [4] Below Figure 10-4
            gamma = 0 if pc else beta * np.mean(z * x)
            # [1] Section 4.1

            def eta(g=gamma):
                return 1 - g

        try:
            # Solve [1] 8.6 and 8.7/8.13/8.23
            d_x, d_y, d_z, d_tau, d_kappa = _get_delta(
                A, b, c, x, y, z, tau, kappa, gamma, eta,
                sparse, lstsq, sym_pos, cholesky, pc, ip, permc_spec)

            if ip:  # initial point
                # [1] 4.4
                # Formula after 8.23 takes a full step regardless if this will
                # take it negative
                alpha = 1.0
                x, y, z, tau, kappa = _do_step(
                    x, y, z, tau, kappa, d_x, d_y,
                    d_z, d_tau, d_kappa, alpha)
                x[x < 1] = 1
                z[z < 1] = 1
                tau = max(1, tau)
                kappa = max(1, kappa)
                ip = False  # done with initial point
            else:
                # [1] Section 4.3
                alpha = _get_step(x, d_x, z, d_z, tau,
                                  d_tau, kappa, d_kappa, alpha0)
                # [1] Equation 8.9
                x, y, z, tau, kappa = _do_step(
                    x, y, z, tau, kappa, d_x, d_y, d_z, d_tau, d_kappa, alpha)

        except (LinAlgError, FloatingPointError,
                ValueError, ZeroDivisionError):
            # this can happen when sparse solver is used and presolve
            # is turned off. Also observed ValueError in AppVeyor Python 3.6
            # Win32 build (PR #8676). I've never seen it otherwise.
            status = 4
            message = _get_message(status)
            break

        # [1] 4.5
        rho_p, rho_d, rho_A, rho_g, rho_mu, obj = _indicators(
            A, b, c, c0, x, y, z, tau, kappa)
        go = rho_p > tol or rho_d > tol or rho_A > tol

        if disp:
            _display_iter(rho_p, rho_d, rho_g, alpha, float(rho_mu), obj)

        # [1] 4.5
        inf1 = (rho_p < tol and rho_d < tol and rho_g < tol and tau < tol *
                max(1, kappa))
        inf2 = rho_mu < tol and tau < tol * min(1, kappa)
        if inf1 or inf2:
            # [1] Lemma 8.4 / Theorem 8.3
            if b.transpose().dot(y) > tol:
                status = 2
            else:  # elif c.T.dot(x) < tol: ? Probably not necessary.
                status = 3
            message = _get_message(status)
            break
        elif iteration >= maxiter:
            status = 1
            message = _get_message(status)
            break

    if disp:
        print(message)

    x_hat = x / tau
    # [1] Statement after Theorem 8.2
    return x_hat, status, message, iteration


def _linprog_ip(
        c,
        A_ub=None,
        b_ub=None,
        A_eq=None,
        b_eq=None,
        bounds=None,
        callback=None,
        alpha0=.99995,
        beta=0.1,
        maxiter=1000,
        disp=False,
        tol=1e-8,
        sparse=False,
        lstsq=False,
        sym_pos=True,
        cholesky=None,
        pc=True,
        ip=False,
        presolve=True,
        permc_spec='MMD_AT_PLUS_A',
        rr=True,
        _sparse_presolve=False,
        **unknown_options):
    r"""
    Minimize a linear objective function subject to linear
    equality constraints, linear inequality constraints, and simple bounds
    using the interior point method of [1]_.

    Linear programming is intended to solve problems of the following form::

        Minimize:     c^T * x

        Subject to:   A_ub * x <= b_ub
                      A_eq * x == b_eq
                      bounds[i][0] < x_i < bounds[i][1]

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
        1-D array of values representing the right hand side of each equality
        constraint (row) in ``A_eq``.
    bounds : sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use ``None`` for one of ``min`` or
        ``max`` when there is no bound in that direction. By default
        bounds are ``(0, None)`` (non-negative).
        If a sequence containing a single tuple is provided, then ``min`` and
        ``max`` will be applied to all variables in the problem.

    Options
    -------
    maxiter : int (default = 1000)
        The maximum number of iterations of the algorithm.
    disp : bool (default = False)
        Set to ``True`` if indicators of optimization status are to be printed
        to the console each iteration.
    tol : float (default = 1e-8)
        Termination tolerance to be used for all termination criteria;
        see [1]_ Section 4.5.
    alpha0 : float (default = 0.99995)
        The maximal step size for Mehrota's predictor-corrector search
        direction; see :math:`\beta_{3}` of [1]_ Table 8.1.
    beta : float (default = 0.1)
        The desired reduction of the path parameter :math:`\mu` (see [3]_)
        when Mehrota's predictor-corrector is not in use (uncommon).
    sparse : bool (default = False)
        Set to ``True`` if the problem is to be treated as sparse after
        presolve. If either ``A_eq`` or ``A_ub`` is a sparse matrix,
        this option will automatically be set ``True``, and the problem
        will be treated as sparse even during presolve. If your constraint
        matrices contain mostly zeros and the problem is not very small (less
        than about 100 constraints or variables), consider setting ``True``
        or providing ``A_eq`` and ``A_ub`` as sparse matrices.
    lstsq : bool (default = False)
        Set to ``True`` if the problem is expected to be very poorly
        conditioned. This should always be left ``False`` unless severe
        numerical difficulties are encountered. Leave this at the default
        unless you receive a warning message suggesting otherwise.
    sym_pos : bool (default = True)
        Leave ``True`` if the problem is expected to yield a well conditioned
        symmetric positive definite normal equation matrix
        (almost always). Leave this at the default unless you receive
        a warning message suggesting otherwise.
    cholesky : bool (default = True)
        Set to ``True`` if the normal equations are to be solved by explicit
        Cholesky decomposition followed by explicit forward/backward
        substitution. This is typically faster for moderate, dense problems
        that are numerically well-behaved.
    pc : bool (default = True)
        Leave ``True`` if the predictor-corrector method of Mehrota is to be
        used. This is almost always (if not always) beneficial.
    ip : bool (default = False)
        Set to ``True`` if the improved initial point suggestion due to [1]_
        Section 4.3 is desired. Whether this is beneficial or not
        depends on the problem.
    presolve : bool (default = True)
        Leave ``True`` if presolve routine should be run. The presolve routine
        is almost always useful because it can detect trivial infeasibilities
        and unboundedness, eliminate fixed variables, and remove redundancies.
        One circumstance in which it might be turned off (set ``False``) is
        when it detects that the problem is trivially unbounded; it is possible
        that that the problem is truly infeasibile but this has not been
        detected.
    rr : bool (default = True)
        Default ``True`` attempts to eliminate any redundant rows in ``A_eq``.
        Set ``False`` if ``A_eq`` is known to be of full row rank, or if you
        are looking for a potential speedup (at the expense of reliability).
    permc_spec : str (default = 'MMD_AT_PLUS_A')
        (Has effect only with ``sparse = True``, ``lstsq = False``, ``sym_pos =
        True``.) A matrix is factorized in each iteration of the algorithm.
        This option specifies how to permute the columns of the matrix for
        sparsity preservation. Acceptable values are:

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering.

        This option can impact the convergence of the
        interior point algorithm; test different values to determine which
        performs best for your problem. For more information, refer to
        ``scipy.sparse.linalg.splu``.

    Returns
    -------
    A ``scipy.optimize.OptimizeResult`` consisting of the following fields:

        x : ndarray
            The independent variable vector which optimizes the linear
            programming problem.
        fun : float
            The optimal value of the objective function
        con : float
            The residuals of the equality constraints (nominally zero).
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
                 4 : Serious numerical difficulties encountered

        nit : int
            The number of iterations performed.
        message : str
            A string descriptor of the exit status of the optimization.

    Notes
    -----

    This method implements the algorithm outlined in [1]_ with ideas from [5]_
    and a structure inspired by the simpler methods of [3]_ and [4]_.

    First, a presolve procedure based on [5]_ attempts to identify trivial
    infeasibilities, trivial unboundedness, and potential problem
    simplifications. Specifically, it checks for:

    - rows of zeros in ``A_eq`` or ``A_ub``, representing trivial constraints;
    - columns of zeros in ``A_eq`` `and` ``A_ub``, representing unconstrained
      variables;
    - column singletons in ``A_eq``, representing fixed variables; and
    - column singletons in ``A_ub``, representing simple bounds.

    If presolve reveals that the problem is unbounded (e.g. an unconstrained
    and unbounded variable has negative cost) or infeasible (e.g. a row of
    zeros in ``A_eq`` corresponds with a nonzero in ``b_eq``), the solver
    terminates with the appropriate status code. Note that presolve terminates
    as soon as any sign of unboundedness is detected; consequently, a problem
    may be reported as unbounded when in reality the problem is infeasible
    (but infeasibility has not been detected yet). Therefore, if the output
    message states that unboundedness is detected in presolve and it is
    necessary to know whether the problem is actually infeasible, set option
    ``presolve=False``.

    If neither infeasibility nor unboundedness are detected in a single pass
    of the presolve check, bounds are tightened where possible and fixed
    variables are removed from the problem. Then, linearly dependent rows
    of the ``A_eq`` matrix are removed, (unless they represent an
    infeasibility) to avoid numerical difficulties in the primary solve
    routine. Note that rows that are nearly linearly dependent (within a
    prescibed tolerance) may also be removed, which can change the optimal
    solution in rare cases. If this is a concern, eliminate redundancy from
    your problem formulation and run with option ``rr=False`` or
    ``presolve=False``.

    Several potential improvements can be made here: additional presolve
    checks outlined in [5]_ should be implemented, the presolve routine should
    be run multiple times (until no further simplifications can be made), and
    more of the efficiency improvements from [2]_ should be implemented in the
    redundancy removal routines.

    After presolve, the problem is transformed to standard form by converting
    the (tightened) simple bounds to upper bound constraints, introducing
    non-negative slack variables for inequality constraints, and expressing
    unbounded variables as the difference between two non-negative variables.

    The primal-dual path following method begins with initial 'guesses' of
    the primal and dual variables of the standard form problem and iteratively
    attempts to solve the (nonlinear) Karush-Kuhn-Tucker conditions for the
    problem with a gradually reduced logarithmic barrier term added to the
    objective. This particular implementation uses a homogeneous self-dual
    formulation, which provides certificates of infeasibility or unboundedness
    where applicable.

    The default initial point for the primal and dual variables is that
    defined in [1]_ Section 4.4 Equation 8.22. Optionally (by setting initial
    point option ``ip=True``), an alternate (potentially improved) starting
    point can be calculated according to the additional recommendations of
    [1]_ Section 4.4.

    A search direction is calculated using the predictor-corrector method
    (single correction) proposed by Mehrota and detailed in [1]_ Section 4.1.
    (A potential improvement would be to implement the method of multiple
    corrections described in [1]_ Section 4.2.) In practice, this is
    accomplished by solving the normal equations, [1]_ Section 5.1 Equations
    8.31 and 8.32, derived from the Newton equations [1]_ Section 5 Equations
    8.25 (compare to [1]_ Section 4 Equations 8.6-8.8). The advantage of
    solving the normal equations rather than 8.25 directly is that the
    matrices involved are symmetric positive definite, so Cholesky
    decomposition can be used rather than the more expensive LU factorization.

    With the default ``cholesky=True``, this is accomplished using
    ``scipy.linalg.cho_factor`` followed by forward/backward substitutions
    via ``scipy.linalg.cho_solve``. With ``cholesky=False`` and
    ``sym_pos=True``, Cholesky decomposition is performed instead by
    ``scipy.linalg.solve``. Based on speed tests, this also appears to retain
    the Cholesky decomposition of the matrix for later use, which is beneficial
    as the same system is solved four times with different right hand sides
    in each iteration of the algorithm.

    In problems with redundancy (e.g. if presolve is turned off with option
    ``presolve=False``) or if the matrices become ill-conditioned (e.g. as the
    solution is approached and some decision variables approach zero),
    Cholesky decomposition can fail. Should this occur, successively more
    robust solvers (``scipy.linalg.solve`` with ``sym_pos=False`` then
    ``scipy.linalg.lstsq``) are tried, at the cost of computational efficiency.
    These solvers can be used from the outset by setting the options
    ``sym_pos=False`` and ``lstsq=True``, respectively.

    Note that with the option ``sparse=True``, the normal equations are solved
    using ``scipy.sparse.linalg.spsolve``. Unfortunately, this uses the more
    expensive LU decomposition from the outset, but for large, sparse problems,
    the use of sparse linear algebra techniques improves the solve speed
    despite the use of LU rather than Cholesky decomposition. A simple
    improvement would be to use the sparse Cholesky decomposition of
    ``CHOLMOD`` via ``scikit-sparse`` when available.

    Other potential improvements for combatting issues associated with dense
    columns in otherwise sparse problems are outlined in [1]_ Section 5.3 and
    [7]_ Section 4.1-4.2; the latter also discusses the alleviation of
    accuracy issues associated with the substitution approach to free
    variables.

    After calculating the search direction, the maximum possible step size
    that does not activate the non-negativity constraints is calculated, and
    the smaller of this step size and unity is applied (as in [1]_ Section
    4.1.) [1]_ Section 4.3 suggests improvements for choosing the step size.

    The new point is tested according to the termination conditions of [1]_
    Section 4.5. The same tolerance, which can be set using the ``tol`` option,
    is used for all checks. (A potential improvement would be to expose
    the different tolerances to be set independently.) If optimality,
    unboundedness, or infeasibility is detected, the solve procedure
    terminates; otherwise it repeats.

    If optimality is achieved, a postsolve procedure undoes transformations
    associated with presolve and converting to standard form. It then
    calculates the residuals (equality constraint violations, which should
    be very small) and slacks (difference between the left and right hand
    sides of the upper bound constraints) of the original problem, which are
    returned with the solution in an ``OptimizeResult`` object.

    References
    ----------
    .. [1] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.
    .. [2] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [3] Freund, Robert M. "Primal-Dual Interior-Point Methods for Linear
           Programming based on Newton's Method." Unpublished Course Notes,
           March 2004. Available 2/25/2017 at
           https://ocw.mit.edu/courses/sloan-school-of-management/15-084j-nonlinear-programming-spring-2004/lecture-notes/lec14_int_pt_mthd.pdf
    .. [4] Fourer, Robert. "Solving Linear Programs by Interior-Point Methods."
           Unpublished Course Notes, August 26, 2005. Available 2/25/2017 at
           http://www.4er.org/CourseNotes/Book%20B/B-III.pdf
    .. [5] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
           programming." Mathematical Programming 71.2 (1995): 221-245.
    .. [6] Bertsimas, Dimitris, and J. Tsitsiklis. "Introduction to linear
           programming." Athena Scientific 1 (1997): 997.
    .. [7] Andersen, Erling D., et al. Implementation of interior point methods
           for large scale linear programming. HEC/Universite de Geneve, 1996.

    """

    _check_unknown_options(unknown_options)

    if callback is not None:
        raise NotImplementedError("method 'interior-point' does not support "
                                  "callback functions.")

    # This is an undocumented option for unit testing sparse presolve
    if _sparse_presolve and A_eq is not None:
        A_eq = sp.sparse.coo_matrix(A_eq)
    if _sparse_presolve and A_ub is not None:
        A_ub = sp.sparse.coo_matrix(A_ub)

    # These should be warnings, not errors
    if not sparse and (sp.sparse.issparse(A_eq) or sp.sparse.issparse(A_ub)):
        sparse = True
        warn("Sparse constraint matrix detected; setting 'sparse':True.",
             OptimizeWarning)

    if sparse and lstsq:
        warn("Invalid option combination 'sparse':True "
             "and 'lstsq':True; Sparse least squares is not recommended.",
             OptimizeWarning)

    if sparse and not sym_pos:
        warn("Invalid option combination 'sparse':True "
             "and 'sym_pos':False; the effect is the same as sparse least "
             "squares, which is not recommended.",
             OptimizeWarning)

    if sparse and cholesky:
        # Cholesky decomposition is not available for sparse problems
        warn("Invalid option combination 'sparse':True "
             "and 'cholesky':True; sparse Colesky decomposition is not "
             "available.",
             OptimizeWarning)

    if lstsq and cholesky:
        warn("Invalid option combination 'lstsq':True "
             "and 'cholesky':True; option 'cholesky' has no effect when "
             "'lstsq' is set True.",
             OptimizeWarning)

    valid_permc_spec = ('NATURAL', 'MMD_ATA', 'MMD_AT_PLUS_A', 'COLAMD')
    if permc_spec.upper() not in valid_permc_spec:
        warn("Invalid permc_spec option: '" + str(permc_spec) + "'. "
             "Acceptable values are 'NATURAL', 'MMD_ATA', 'MMD_AT_PLUS_A', "
             "and 'COLAMD'. Reverting to default.",
             OptimizeWarning)
        permc_spec = 'MMD_AT_PLUS_A'

    # This can be an error
    if not sym_pos and cholesky:
        raise ValueError(
            "Invalid option combination 'sym_pos':False "
            "and 'cholesky':True: Cholesky decomposition is only possible "
            "for symmetric positive definite matrices.")

    cholesky = cholesky is None and sym_pos and not sparse and not lstsq

    iteration = 0
    complete = False    # will become True if solved in presolve
    undo = []

    # Convert lists to numpy arrays, etc...
    c, A_ub, b_ub, A_eq, b_eq, bounds = _clean_inputs(
        c, A_ub, b_ub, A_eq, b_eq, bounds)

    # Keep the original arrays to calculate slack/residuals for original
    # problem.
    c_o, A_ub_o, b_ub_o, A_eq_o, b_eq_o = c.copy(
    ), A_ub.copy(), b_ub.copy(), A_eq.copy(), b_eq.copy()

    # Solve trivial problem, eliminate variables, tighten bounds, etc...
    c0 = 0  # we might get a constant term in the objective
    if presolve is True:
        (c, c0, A_ub, b_ub, A_eq, b_eq, bounds, x, undo, complete, status,
            message) = _presolve(c, A_ub, b_ub, A_eq, b_eq, bounds, rr)

    # If not solved in presolve, solve it
    if not complete:
        # Convert problem to standard form
        A, b, c, c0 = _get_Abc(c, c0, A_ub, b_ub, A_eq, b_eq, bounds, undo)
        # Solve the problem
        x, status, message, iteration = _ip_hsd(A, b, c, c0, alpha0, beta,
                                                maxiter, disp, tol, sparse,
                                                lstsq, sym_pos, cholesky,
                                                pc, ip, permc_spec)

    # Eliminate artificial variables, re-introduce presolved variables, etc...
    # need modified bounds here to translate variables appropriately
    x, fun, slack, con, status, message = _postprocess(
        x, c_o, A_ub_o, b_ub_o, A_eq_o, b_eq_o,
        bounds, complete, undo, status, message, tol)

    sol = {
        'x': x,
        'fun': fun,
        'slack': slack,
        'con': con,
        'status': status,
        'message': message,
        'nit': iteration,
        "success": status == 0}

    return OptimizeResult(sol)
