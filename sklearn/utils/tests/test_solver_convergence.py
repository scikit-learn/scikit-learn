import pytest

from ..solver_convergence import _check_convergence_params
#from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
#from sklearn.svm import LinearSVC, LinearSVR

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
ALL_SOLVERS = ['liblinear', 'lbfgs', 'sag', 'saga', 'newton-cg']

@pytest.mark.parametrize("solver", ALL_SOLVERS)
def test_ok_check_convergence_params(solver):
    # Check if when tol and max_iter are auto, they are consequently and 
    # automatically set
    tol = 'auto'
    max_iter = 'auto'
    default_tol = DEFAULT_CONV[solver][0][0]
    default_max_iter = DEFAULT_CONV[solver][0][1]
    assert((default_tol, default_max_iter) == _check_convergence_params(solver,
            tol, max_iter))

    # Check if inputed tol and max_iter are correctly returned
    tol = 1e-2
    max_iter = 10000
    assert((tol, max_iter) == _check_convergence_params(solver, tol, max_iter))


@pytest.mark.parametrize("solver", ALL_SOLVERS)
def test_ko_check_convergence_params(solver):
    # Check if when tol is inputed and max_iter is auto an error is raised
    tol = 1e-4
    max_iter = 'auto'
    with pytest.raises(ValueError):
        _check_convergence_params(solver, tol, max_iter)

    # Check if when tol is auto and max_iter is inputed an error is raised
    tol = 'auto'
    max_iter = 100
    with pytest.raises(ValueError):
        _check_convergence_params(solver, tol, max_iter)


@pytest.mark.parametrize("internal_solver",
        list(DEFAULT_CONV['liblinear'].keys()))
def test_internal_solver_check_convergence_params(internal_solver):
    solver = 'liblinear'
    default_tol, default_max_iter = DEFAULT_CONV[solver][internal_solver]
    assert((default_tol, default_max_iter),
           _check_convergence_params('liblinear', 'auto', 'auto',
                                     internal_solver))
