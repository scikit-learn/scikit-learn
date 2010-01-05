from ctypes import POINTER, c_double, c_int

from kernel import *
import libsvm

__all__ = [
    'LibSvmModel'
    ]

class LibSvmModel:
    def __init__(self, kernel,
                 tolerance=0.001, shrinking=True, cache_size=40):
        """
        Parameters:

        - `kernel`: XXX
        - `tolerance`: tolerance of termination criterion
        - `shrinking`: whether to use the shrinking heuristics
        - `cache_size` kernel evaluation cache size (MB)
        """
        self.kernel = kernel
        self.tolerance = tolerance
        self.shrinking = shrinking
        self.cache_size = cache_size

        # kernel parameters
        param = libsvm.svm_parameter()
        param.kernel_type = kernel.kernel_type
        param.degree = getattr(kernel, 'degree', 0)
        param.gamma = getattr(kernel, 'gamma', 0.0)
        param.coef0 = getattr(kernel, 'coef0', 0.0)

        # other parameters
        param.eps = tolerance
        param.shrinking = shrinking
        param.cache_size = cache_size

        # defaults for optional parameters
        param.nr_weight = 0
        param.weight = POINTER(c_double)()
        param.weight_label = POINTER(c_int)()
        param.probability = False

        self.param = param

    def fit(self, dataset, ResultType, PredictorType):
        problem = dataset._create_svm_problem()
        dataset._update_svm_parameter(self.param)
        self._check_problem_param(problem, self.param)
        model = libsvm.svm_train(problem, self.param)
        return ResultType(model, dataset, PredictorType)

    def _check_problem_param(self, problem, param):
        error_msg = libsvm.svm_check_parameter(problem, param)
        if error_msg:
            raise ValueError, error_msg
