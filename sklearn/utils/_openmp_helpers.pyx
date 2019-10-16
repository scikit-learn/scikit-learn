def _openmp_supported():
    """Determines whether scikit-learn has been built with OpenMP support
    
    It allows to retrieve at runtime the information gathered at compile time.
    """
    # SKLEARN_OPENMP_SUPPORTED is resolved at compile time during cythonization
    # It is defined via the `compile_time_env` kwarg of the `cythonize` call
    # and behaves like the `-D` option of the C preprocessor.
    return SKLEARN_OPENMP_SUPPORTED
