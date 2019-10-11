def _openmp_supported():
    """
    Return True if OpenMP is not supported and OpenMP support has not been
    explicitly disabled during the build
    """
    return SKLEARN_OPENMP_SUPPORTED
