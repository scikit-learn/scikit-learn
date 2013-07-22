cimport numpy as cnp
import numpy as np
ctypedef cnp.float64_t FLOAT_t
ctypedef cnp.intp_t INT_t
ctypedef cnp.ulong_t INDEX_t
ctypedef cnp.uint8_t BOOL_t
from _basis cimport Basis
from _record cimport ForwardPassRecord

ctypedef enum StoppingCondition:
    MAXTERMS=0,
    MAXRSQ=1,
    NOIMPRV=2,
    LOWGRSQ=3,
    NOCAND=4

cdef class ForwardPasser:

	#User selected parameters
    cdef int endspan
    cdef int minspan
    cdef FLOAT_t endspan_alpha
    cdef FLOAT_t minspan_alpha
    cdef int max_terms
    cdef int max_degree
    cdef FLOAT_t thresh
    cdef FLOAT_t penalty
    cdef int check_every
    cdef int min_search_points
    cdef list xlabels
    cdef FLOAT_t zero_tol
    
    #Input data
    cdef cnp.ndarray X
    cdef cnp.ndarray y
    cdef cnp.ndarray weights
    cdef INDEX_t m
    cdef INDEX_t n
    cdef FLOAT_t sst
    cdef FLOAT_t y_squared
    
    #Working floating point data 
    cdef cnp.ndarray B #Data matrix in basis space
    cdef cnp.ndarray B_orth #Orthogonalized version of B
    cdef cnp.ndarray c
    cdef cnp.ndarray norms
    cdef cnp.ndarray u
    cdef cnp.ndarray B_orth_times_parent_cum
    cdef FLOAT_t c_squared
    
    #Working integer data
    cdef cnp.ndarray sort_tracker
    cdef cnp.ndarray sorting
    cdef cnp.ndarray mwork
#    cdef cnp.ndarray x_order
    cdef cnp.ndarray linear_variables
    
    #Object construction
    cdef ForwardPassRecord record
    cdef Basis basis
    
    cpdef Basis get_basis(ForwardPasser self)
    
#    cpdef init_x_order(ForwardPasser self)
    
    cpdef init_linear_variables(ForwardPasser self)
    
    cpdef run(ForwardPasser self)
    
    cdef stop_check(ForwardPasser self)
    
    cpdef int orthonormal_update(ForwardPasser self, INDEX_t k)
    
    cpdef orthonormal_downdate(ForwardPasser self, INDEX_t k)
    
    cdef next_pair(ForwardPasser self)
    
    cdef best_knot(ForwardPasser self, INDEX_t parent, INDEX_t variable, INDEX_t k, cnp.ndarray[INT_t,ndim=1] candidates, cnp.ndarray[INT_t,ndim=1] order, FLOAT_t * mse, FLOAT_t * knot, INDEX_t * knot_idx)

    