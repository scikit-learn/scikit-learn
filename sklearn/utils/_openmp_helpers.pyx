def _openmp_supported():
    """Determines whether scikit-learn has been built with OpenMP support"""
    # SKLEARN_OPENMP_SUPPORTED is resolved at compile time during cythonization
    return SKLEARN_OPENMP_SUPPORTED
