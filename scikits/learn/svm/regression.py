__all__ = [
    'LibSvmEpsilonRegressionModel',
    'LibSvmNuRegressionModel'
    ]

from ctypes import cast, POINTER

from model import LibSvmModel
import libsvm

# XXX document why get_svr_probability could be useful

class LibSvmRegressionResults:
    def __init__(self, model, dataset):
        self.model = model
        self.dataset = dataset
        model = model.contents
        self.rho = model.rho[0]
        self.sv_coef = model.sv_coef[0][:model.l]
        self.sigma = model.probA[0]

    def __del__(self):
        libsvm.svm_destroy_model(self.model)

    def predict(self, dataset):
        """
        This function does regression on a test vector x and returns
        the function value of x calculated using the model.
        """
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            return libsvm.svm_predict(self.model, xptr)
        return map(p, dataset.data)

    def get_svr_probability(self):
        """
        This function returns a value sigma > 0, which is a parameter
        in the probability model for the test data.

        For test data, we consider the probability model:

        target value = predicted value + z

        where z is the Laplace distribution: e^(-|x|/sigma)/(2*sigma).
        """
        return self.sigma

class LibSvmEpsilonRegressionModel(LibSvmModel):
    """
    A model for epsilon-SV regression.

    See also:

    - Smola, Schoelkopf. A Tutorial on Support Vector Regression.
    - Gunn. Support Vector Machines for Classification and Regression.
    - Mueller, Vapnik. Using Support Vector Machines for Time Series
      Prediction.
    """

    Results = LibSvmRegressionResults

    def __init__(self, kernel, epsilon=0.1, cost=1.0, **kwargs):
        LibSvmModel.__init__(self, kernel, **kwargs)
        self.epsilon = epsilon
        self.cost = cost
        self.param.svm_type = libsvm.EPSILON_SVR
        self.param.p = epsilon
        self.param.C = cost
        self.param.probability = True

class LibSvmNuRegressionModel(LibSvmModel):
    """
    A model for nu-SV regression.

    See also: Schoelkopf, et al. New Support Vector Algorithms.
    """

    Results = LibSvmRegressionResults

    def __init__(self, kernel, nu=0.5, cost=1.0, **kwargs):
        LibSvmModel.__init__(self, kernel, **kwargs)
        self.nu = nu
        self.cost = cost
        self.param.svm_type = libsvm.NU_SVR
        self.param.nu = nu
        self.param.C = cost
        self.param.probability = True
