# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from HConst cimport HighsBasisStatus

cdef extern from "HighsLp.h" nogil:
    # From HiGHS/src/lp_data/HighsLp.h
    cdef cppclass HighsLp:
        int numCol_
        int numRow_

        vector[int] Astart_
        vector[int] Aindex_
        vector[double] Avalue_
        vector[double] colCost_
        vector[double] colLower_
        vector[double] colUpper_
        vector[double] rowLower_
        vector[double] rowUpper_

        ObjSense sense_
        double offset_

        string model_name_
        string lp_name_

        vector[string] row_names_
        vector[string] col_names_

        vector[int] integrality_

    ctypedef enum HighsModelStatus:
        HighsModelStatusNOTSET "HighsModelStatus::NOTSET" = 0
        HighsModelStatusHIGHS_MODEL_STATUS_MIN "HighsModelStatus::HIGHS_MODEL_STATUS_MIN" = HighsModelStatusNOTSET
        HighsModelStatusLOAD_ERROR "HighsModelStatus::LOAD_ERROR"
        HighsModelStatusMODEL_ERROR "HighsModelStatus::MODEL_ERROR"
        HighsModelStatusPRESOLVE_ERROR "HighsModelStatus::PRESOLVE_ERROR"
        HighsModelStatusSOLVE_ERROR "HighsModelStatus::SOLVE_ERROR"
        HighsModelStatusPOSTSOLVE_ERROR "HighsModelStatus::POSTSOLVE_ERROR"
        HighsModelStatusMODEL_EMPTY "HighsModelStatus::MODEL_EMPTY"
        HighsModelStatusPRIMAL_INFEASIBLE "HighsModelStatus::PRIMAL_INFEASIBLE"
        HighsModelStatusPRIMAL_UNBOUNDED "HighsModelStatus::PRIMAL_UNBOUNDED"
        HighsModelStatusOPTIMAL "HighsModelStatus::OPTIMAL"
        HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND "HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND"
        HighsModelStatusREACHED_TIME_LIMIT "HighsModelStatus::REACHED_TIME_LIMIT"
        HighsModelStatusREACHED_ITERATION_LIMIT "HighsModelStatus::REACHED_ITERATION_LIMIT"
        HighsModelStatusHIGHS_MODEL_STATUS_MAX "HighsModelStatus::HIGHS_MODEL_STATUS_MAX" = HighsModelStatusREACHED_ITERATION_LIMIT

    ctypedef enum ObjSense:
        ObjSenseMINIMIZE "ObjSense::MINIMIZE" = 1
        ObjSenseMAXIMIZE "ObjSense::MAXIMIZE" = -1

    cdef cppclass HighsSolution:
        vector[double] col_value
        vector[double] col_dual
        vector[double] row_value
        vector[double] row_dual

    cdef cppclass HighsBasis:
        bool valid_
        vector[HighsBasisStatus] col_status
        vector[HighsBasisStatus] row_status
