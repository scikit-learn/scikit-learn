# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

cpdef enum:
    # N.B.: Negative values indicate constraints, postive values indicate penalties
    # No OPeration / do nothing
    NOP = 1987

    # Lasso: \sum_j ||Wj||_1
    L1_PENALTY = 1
    L11_PENALTY = 1

    # Group Lasso: penalty = \sum_j ||Wj||_2
    L21_PENALTY = 21

    # L2 ball constraint: ||Wj|_2 \le reg for all j in [0...n_tasks - 1]
    L2_CONSTRAINT = -2
    L2INF_CONSTRAINT = -2

    # L1 ball constraint: ||Wj|_1 \le reg for all j in [0...n_tasks - 1]
    L1INF_CONSTRAINT = -1
    L1_CONSTRAINT = -1

