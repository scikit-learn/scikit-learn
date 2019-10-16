def _openmp_supported():
    """Determines whether scikit-learn has been built with OpenMP support"""
    # SKLEARN_OPENMP_SUPPORTED is resolved at compile time during cythonization
    # The cythonize parameter `compile_time_env` works like a #define in C.
    return SKLEARN_OPENMP_SUPPORTED
