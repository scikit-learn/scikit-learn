# distutils: language=c++
# cython: language_level=3

from libcpp.string cimport string

from HighsLp cimport HighsModelStatus

cdef extern from "HighsModelUtils.h" nogil:
    # From HiGHS/src/lp_data/HighsModelUtils.h
    string utilHighsModelStatusToString(const HighsModelStatus model_status)
