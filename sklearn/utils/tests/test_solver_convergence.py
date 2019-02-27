from ..solver_convergence import _check_convergence_params
from sklearn.utils.testing import assert_raises


def test_check_convergence_params():
    default_tol = 1e-4
    default_max_iter = 1000
    all_solvers = ['liblinear', 'newton-cg', 'lbfgs', 'sag', 'saga']

    # Check if inputed tol and max_iter are correctly returned
    tol = 1e-3
    max_iter = 100
    for solver in all_solvers :
        assert((tol, max_iter) == 
               _check_convergence_params(tol, max_iter, solver))
    
    # Check if when tol and max_iter are auto, they are consequently and 
    # automatically set
    tol = 'auto'
    max_iter = 'auto'
    for solver in all_solvers :
        assert((default_tol, default_max_iter) == 
               _check_convergence_params(tol, max_iter, solver))


def test_errors_check_convergence_params():
    all_solvers = ['liblinear', 'newton-cg', 'lbfgs', 'sag', 'saga']

    # Check if when tol is inputed and max_iter is auto an error is raised
    tol = 1e-4
    max_iter = 'auto'
    for solver in all_solvers :
        assert_raises(ValueError, _check_convergence_params, tol,
                      max_iter, solver)

    # Check if when tol is auto and max_iter is inputed an error is raised
    tol = 'auto'
    max_iter = 100
    for solver in all_solvers :
        assert_raises(ValueError, _check_convergence_params, tol,
                      max_iter, solver)