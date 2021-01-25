# distutils: language=c++
# cython: language_level=3

cdef extern from "HConst.h" nogil:

    const int HIGHS_CONST_I_INF
    const double HIGHS_CONST_INF
    const double HIGHS_CONST_TINY
    const double HIGHS_CONST_ZERO
    const int HIGHS_THREAD_LIMIT

    cdef enum HighsBasisStatus:
        LOWER "HighsBasisStatus::LOWER" = 0, # (slack) variable is at its lower bound [including fixed variables]
        BASIC "HighsBasisStatus::BASIC" # (slack) variable is basic
        UPPER "HighsBasisStatus::UPPER" # (slack) variable is at its upper bound
        ZERO "HighsBasisStatus::ZERO" # free variable is non-basic and set to zero
        NONBASIC "HighsBasisStatus::NONBASIC" # nonbasic with no specific bound information - useful for users and postsolve
        SUPER "HighsBasisStatus::SUPER"

    cdef enum SolverOption:
        SOLVER_OPTION_SIMPLEX "SolverOption::SOLVER_OPTION_SIMPLEX" = -1
        SOLVER_OPTION_CHOOSE "SolverOption::SOLVER_OPTION_CHOOSE"
        SOLVER_OPTION_IPM "SolverOption::SOLVER_OPTION_IPM"

    cdef enum PrimalDualStatus:
        PrimalDualStatusSTATUS_NOT_SET "PrimalDualStatus::STATUS_NOT_SET" = -1
        PrimalDualStatusSTATUS_MIN "PrimalDualStatus::STATUS_MIN" = PrimalDualStatusSTATUS_NOT_SET
        PrimalDualStatusSTATUS_NO_SOLUTION "PrimalDualStatus::STATUS_NO_SOLUTION"
        PrimalDualStatusSTATUS_UNKNOWN "PrimalDualStatus::STATUS_UNKOWN"
        PrimalDualStatusSTATUS_INFEASIBLE_POINT "PrimalDualStatus::STATUS_INFEASIBLE_POINT"
        PrimalDualStatusSTATUS_FEASIBLE_POINT "PrimalDualStatus::STATUS_FEASIBLE_POINT"
        PrimalDualStatusSTATUS_MAX "PrimalDualStatus::STATUS_MAX" = PrimalDualStatusSTATUS_FEASIBLE_POINT

    cdef enum HighsOptionType:
        HighsOptionTypeBOOL "HighsOptionType::BOOL" = 0
        HighsOptionTypeINT "HighsOptionType::INT"
        HighsOptionTypeDOUBLE "HighsOptionType::DOUBLE"
        HighsOptionTypeSTRING "HighsOptionType::STRING"
