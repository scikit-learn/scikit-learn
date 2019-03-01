DEFAULT_CONV = {'liblinear': {0: (1e-2, 1000), 1: (1e-1, 1000), 2: (1e-2, 1000),
                              3: (1e-1, 1000), 4: (1e-1, 1000), 5: (1e-2, 1000),
                              6: (1e-2, 1000), 7: (1e-1, 1000),
                              11: (1e-3, 1000), 12: (1e-1, 1000),
                              13: (1e-1, 1000)},
                'lbfgs': {0: (1e-4, 15000)},
                'sag': {0: (1e-3, 1000)},
                'saga': {0: (1e-3, 1000)},
                'newton-cg': {0: (1e-4, 100)},
                'logistic': {0: (1e-4, 1e2)},
                'svm': {0: (1e-4, 1e3)}}


def _check_convergence_params(solver, tol, max_iter):
    if ((tol == 'auto' and max_iter != 'auto') or 
            (tol != 'auto' and max_iter == 'auto')):
        raise ValueError("'auto' value on tol and max_iter can only be "
                            "used simultaneously. Either provide a numeric "
                            "for both values or set both to 'auto'")
    
    if tol == 'auto' and max_iter == 'auto':
        tol = DEFAULT_CONV[solver][0][0]
        max_iter = DEFAULT_CONV[solver][0][1]

    return tol, max_iter