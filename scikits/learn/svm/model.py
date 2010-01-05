__all__ = [
    'LibSvmModel'
    ]

from ctypes import *

from kernel import *
import libsvm

class LibSvmModel:
    def __init__(self, kernel, tolerance=0.001, shrinking=True, cache_size=40):
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

        param = libsvm.svm_parameter()

        if isinstance(kernel, LinearKernel):
            param.kernel_type = libsvm.LINEAR
        elif isinstance(kernel, PolynomialKernel):
            param.kernel_type = libsvm.POLY
            param.degree = kernel.degree
            param.gamma = kernel.gamma
            param.coef0 = kernel.coef0
        elif isinstance(kernel, RBFKernel):
            param.kernel_type = libsvm.RBF
            param.gamma = kernel.gamma
        elif isinstance(kernel, SigmoidKernel):
            param.kernel_type = libsvm.SIGMOID
            param.gamma = kernel.gamma
            param.coef0 = kernel.coef0
        else:
            raise ValueError, 'unknown kernel type'

        param.eps = tolerance
        param.shrinking = shrinking
        param.cache_size = cache_size
        # set defaults for optional parameters
        param.nr_weight = 0
        param.weight = None
        param.weight_label = None
        param.probability = False

        self.param = param

    def fit(self, dataset):
        # XXX don't poke around in dataset's internals

        # no reference to the svm_problem is kept because a svm_model
        # only requires some parameters and the support vectors chosen
        # from the dataset
        problem = libsvm.svm_problem()
        problem.l = len(dataset.data)
        y = (c_double*problem.l)()
        x = (POINTER(libsvm.svm_node)*problem.l)()
        for i, (yi, xi) in enumerate(dataset.data):
            y[i] = yi
            x[i] = cast(xi.ctypes.data, POINTER(libsvm.svm_node))
        problem.x = cast(addressof(x), POINTER(POINTER(libsvm.svm_node)))
        problem.y = cast(addressof(y), POINTER(c_double))
        self._check_problem_param(problem, self.param)

        model = libsvm.svm_train(problem, self.param)

        # XXX because libsvm only does a shallow copy of the
        # svm_parameter into the model, we have to make sure that a
        # reference to weight labels and weights are kept somewhere

        return self.Results(model, dataset)

    def _check_problem_param(self, problem, param):
        error_msg = libsvm.svm_check_parameter(problem, param)
        if error_msg:
            raise ValueError, error_msg
