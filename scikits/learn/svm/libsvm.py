__all__ = [
    'svm_node_dtype',
    'C_SVC',
    'NU_SVC',
    'ONE_CLASS',
    'EPSILON_SVR',
    'NU_SVR',
    'LINEAR',
    'POLY',
    'RBF',
    'SIGMOID',
    'PRECOMPUTED',
    'svm_node',
    'svm_parameter',
    'svm_problem',
    'svm_model',
    'svm_check_parameter',
    'svm_train',
    'svm_check_probability_model',
    'svm_predict',
    'svm_predict_values',
    'svm_predict_probability',
    'svm_get_svr_probability',
    'svm_cross_validation',
    'svm_destroy_model',
    'svm_destroy_param'
    ]

import numpy as N

svm_node_dtype = \
    N.dtype({'names' : ['index', 'value'],
             'formats' : [N.intc, N.float64]},
            align=1)

import utils
_libsvm = utils.load_ctypes_library('libsvm_', __file__)

from ctypes import c_int, c_double, POINTER, Structure, c_char_p

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
        # parameter weight is weight*C, of C-SVC and nu-SVC
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
        ('param', svm_parameter),
        ('nr_class', c_int),
        ('l', c_int),
        # support vectors (length el)
        ('SV', POINTER(POINTER(svm_node))),
        # length nr_class-1, each length el
        ('sv_coef', POINTER(POINTER(c_double))),
        # length nr_class*(nr_class-1)/2
        ('rho', POINTER(c_double)),
        # length nr_class*(nr_class-1)/2 for classification
        # length 1 for regression
        # length 0 for one-class
        ('probA', POINTER(c_double)),
        # ??? for classification (length nr_class*(nr_class-1)/2)
        ('probB', POINTER(c_double)),
        # support vectors labels for classifcation (length nr_class)
        ('labels', POINTER(c_int)),
        # length nr_class; nSV[0] + nSV[1] + ... + nSV[n-1] = el
        ('nSV', POINTER(c_int)),
        # whether svm_destroy_model should free support vectors
        ('free_sv', c_int)
        ]

# svm_check_parameter
_libsvm.svm_check_parameter.restype = c_char_p
_libsvm.svm_check_parameter.argtypes = \
    [POINTER(svm_problem), POINTER(svm_parameter)]
svm_check_parameter = _libsvm.svm_check_parameter

# svm_train
_libsvm.svm_train.argtypes = \
    [POINTER(svm_problem), POINTER(svm_parameter)]
_libsvm.svm_train.restype = POINTER(svm_model)
svm_train = _libsvm.svm_train

# svm_check_probability_model
_libsvm.svm_check_probability_model.argtypes = [POINTER(svm_model)]
svm_check_probability_model = _libsvm.svm_check_probability_model

# svm_predict
_libsvm.svm_predict.argtypes = [POINTER(svm_model), POINTER(svm_node)]
_libsvm.svm_predict.restype = c_double
svm_predict = _libsvm.svm_predict

# svm_predict_values
_libsvm.svm_predict_values.argtypes = \
    [POINTER(svm_model), POINTER(svm_node), POINTER(c_double)]
svm_predict_values = _libsvm.svm_predict_values

# svm_predict_probability
_libsvm.svm_predict_probability.restype = c_double
_libsvm.svm_predict_probability.argtypes = \
    [POINTER(svm_model), POINTER(svm_node), POINTER(c_double)]
svm_predict_probability = _libsvm.svm_predict_probability

# svm_get_svr_probability
_libsvm.svm_get_svr_probability.restype = c_double
_libsvm.svm_get_svr_probability.argtypes = [POINTER(svm_model)]
svm_get_svr_probability = _libsvm.svm_get_svr_probability

# svm_cross_validation
_libsvm.svm_cross_validation.argtypes = [
    POINTER(svm_problem),
    POINTER(svm_parameter),
    c_int,
    POINTER(c_double)
    ]
svm_cross_validation = _libsvm.svm_cross_validation

# svm_destroy_model
_libsvm.svm_destroy_model.argtypes = [POINTER(svm_model)]
svm_destroy_model = _libsvm.svm_destroy_model

# svm_destroy_param
_libsvm.svm_destroy_param.argtypes = [POINTER(svm_parameter)]
svm_destroy_param = _libsvm.svm_destroy_param
