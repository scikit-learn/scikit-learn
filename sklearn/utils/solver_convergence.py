def _check_convergence_params (tol, max_iter, solver='liblinear'):
    default_tol = 1e-4
    default_max_iter = 1000
    if ((tol == 'auto' and max_iter != 'auto') or 
            (tol != 'auto' and max_iter == 'auto')):
        raise ValueError("'auto' value on tol and max_iter can only be "
                            "used simultaneously. Either provide a numeric "
                            "for both values or set both to 'auto'")
    
    if tol == 'auto' and max_iter == 'auto':
        tol = default_tol
        max_iter = default_max_iter

    return tol, max_iter