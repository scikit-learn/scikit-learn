import inspect

from ctypes import *
import numpy as N

__all__ = [
    'svm_node_dtype'
    ]

_libsvm = N.ctypes_load_library('libsvm_', __file__)

svm_node_dtype = \
    N.dtype({'names' : ['index', 'value'],
             'formats' : [N.intc, N.float64]},
            align=True)

# svm types
C_SVC = 0
NU_SVC = 1
ONE_CLASS = 2
EPSILON_SVR = 3
NU_SVR = 4

# kernel types
LINEAR = 0
POLY = 1
RBF = 2
SIGMOID = 3
PRECOMPUTED = 4

class svm_node(Structure):
    _fields_ = [
        ('index', c_int),
        ('value', c_double)
        ]

class svm_parameter(Structure):
    _fields_ = [
	('svm_type', c_int),
	('kernel_type', c_int),
        # degree in polynomial kernel function
        ('degree', c_int),
        # gamma in kernel function
        ('gamma', c_double),
        # coef0 in kernel function
        ('coef0', c_double),
        # kernel evaluation cache memory size in MB
        ('cache_size', c_double),
        # tolerance of termination criterion
        ('eps', c_double),
        # parameter C of C-SVC, epsilon-SVR and nu-SVR
        ('C', c_double),
        ('nr_weight', c_int),
        ('weight_label', POINTER(c_int)),
        # parameter weight of C-SVC and nu-SVC
        # actual weight of class is this weight * C
        ('weight', POINTER(c_double)),
        # parameter nu of nu-SVC, one-class SVM and nu-SVR
        ('nu', c_double),
        # epsilon in loss function of epsilon-SVR
        ('p', c_double),
        # whether to use the shrinking heuristics
	('shrinking', c_int),
        # whether to train a SVC or SVR model for probability estimates
        ('probability', c_int)
        ]

class svm_problem(Structure):
    _fields_ = [
        # length (number of elements)
        ('l', c_int),
        # label for classification or response variable for regression
        ('y', POINTER(c_double)),
        ('x', POINTER(POINTER(svm_node)))
        ]

class svm_model(Structure):
    _fields_ = [
        # parameters used to train the model
        ('param', svm_parameter),
        ('nr_class', c_int),
        # el
        ('l', c_int),
        # support vectors (length el)
        ('SV', POINTER(POINTER(svm_node))),
        # length nr_class-1, each length el
        ('sv_coef', POINTER(POINTER(c_double))),
        # length nr_class*(nr_class-1)/2
        ('rho', POINTER(c_double)),
        # length nr_class*(nr_class-1)/2 for classification
        # length 1 for regression
        ('probA', POINTER(c_double)),
        # length nr_class*(nr_class-1)/2 for classification
        ('probB', POINTER(c_double)),
        # support vectors labels for classifcation (length nr_class)
        ('labels', POINTER(c_int)),
        # length nr_class; nSV[0] + nSV[1] + ... + nSV[n-1] = el
        ('nSV', POINTER(c_int)),
        # whether svm_destroy_model should free support vectors
        ('free_sv', c_int)
        ]

libsvm_api = {
    'svm_check_parameter' :
    (c_char_p, [POINTER(svm_problem), POINTER(svm_parameter)]),
    'svm_train' :
    (POINTER(svm_model), [POINTER(svm_problem), POINTER(svm_parameter)]),
    'svm_check_probability_model' :
    (None, [POINTER(svm_model)]),
    'svm_predict' :
    (c_double, [POINTER(svm_model), POINTER(svm_node)]),
    'svm_predict_values' :
    (None, [POINTER(svm_model), POINTER(svm_node), POINTER(c_double)]),
    'svm_predict_probability' :
    (c_double, [POINTER(svm_model), POINTER(svm_node), POINTER(c_double)]),
    'svm_get_svr_probability' :
    (c_double, [POINTER(svm_model)]),
    'svm_cross_validation' :
    (None,
     [POINTER(svm_problem), POINTER(svm_parameter), c_int,
      POINTER(c_double)]),
    'svm_destroy_model' :
    (None, [POINTER(svm_model)])
    }

for f, (restype, argtypes) in libsvm_api.iteritems():
    func = getattr(_libsvm, f)
    func.restype = restype
    func.argtypes = argtypes
    inspect.currentframe().f_locals[f] = func

def create_svm_problem(data):
    problem = svm_problem()
    problem.l = len(data)
    y = (c_double*problem.l)()
    x = (POINTER(svm_node)*problem.l)()
    for i, (yi, xi) in enumerate(data):
        y[i] = yi
        x[i] = xi.ctypes.data_as(POINTER(svm_node))
    problem.x = x
    problem.y = y
    return problem
