def _openmp_warn_unsupported():
    """
    Return True if OpenMP is not supported and OpenMP support has not been
    explicitly disabled during the build
    """
    return not (SKLEARN_OPENMP_SUPPORTED or SKLEARN_OPENMP_EXPLICITLY_DISABLED)
