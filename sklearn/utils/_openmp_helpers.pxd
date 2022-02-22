# Helpers to access OpenMP threads information
#
# Those interfaces act as indirections which allows the non-support of OpenMP
# for implementations which have been written for it.

cdef int _openmp_thread_num() nogil
