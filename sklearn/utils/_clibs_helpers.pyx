cimport openmp

from ._joblib import effective_n_jobs


cpdef int effective_n_jobs_openmp(n_jobs=None):
    """Determine the number of threads for OpenMP
    
    The behavior differs from joblib's effective_n_jobs for n_jobs=None in
    which case effective_n_jobs_openmp returns openmp.omp_get_max_threads(). 
    """
    if n_jobs is None:
        return openmp.omp_get_max_threads()
    
    return effective_n_jobs(n_jobs)