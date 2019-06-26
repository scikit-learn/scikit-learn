def openmp_status():
    if SKLEARN_OPENMP_SUPPORTED:
        return "supported"
    elif OPENMP_EXPLICIT_DISABLED:
        return "explicit disabled"
    return "unsupported"
