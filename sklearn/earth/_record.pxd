cimport numpy as cnp
ctypedef cnp.float64_t FLOAT_t
ctypedef cnp.intp_t INT_t
ctypedef cnp.ulong_t INDEX_t
ctypedef cnp.uint8_t BOOL_t
from _basis cimport Basis

cdef class Record:
    cdef list iterations
    cdef int num_samples
    cdef int num_variables
    cdef FLOAT_t penalty
    cdef FLOAT_t sst #Sum of squares total
    
    cpdef append(Record self, Iteration iteration)
    
    cpdef FLOAT_t mse(Record self, INDEX_t iteration)
    
    cpdef FLOAT_t rsq(Record self, INDEX_t iteration)
    
    cpdef FLOAT_t gcv(Record self, INDEX_t iteration)
    
    cpdef FLOAT_t grsq(Record self, INDEX_t iteration)

cdef class PruningPassRecord(Record):
    cdef INDEX_t selected

    cpdef set_selected(PruningPassRecord self, INDEX_t selected)

    cpdef INDEX_t get_selected(PruningPassRecord self)

    cpdef roll_back(PruningPassRecord self, Basis basis)
	
cdef class ForwardPassRecord(Record):
    cdef int stopping_condition
    
    cpdef set_stopping_condition(ForwardPassRecord self, int stopping_condition)
    
cdef class Iteration:
    cdef FLOAT_t mse
    cdef INDEX_t size
    
    cpdef FLOAT_t get_mse(Iteration self)
    
    cpdef INDEX_t get_size(Iteration self)
    
cdef class PruningPassIteration(Iteration):
    cdef INDEX_t pruned
    
    cpdef INDEX_t get_pruned(PruningPassIteration self)
    
cdef class FirstPruningPassIteration(PruningPassIteration):
    pass
    
cdef class ForwardPassIteration(Iteration):
    cdef INDEX_t parent
    cdef INDEX_t variable
    cdef FLOAT_t knot
    cdef int code
    cdef bint no_candidates
    
    cpdef set_no_candidates(ForwardPassIteration self, bint value)
        
    cpdef no_further_candidates(ForwardPassIteration self)
    
cdef class FirstForwardPassIteration(ForwardPassIteration):
    cpdef INDEX_t get_size(FirstForwardPassIteration self)