from ctypes import cast, POINTER, c_double
import numpy as N

from model import LibSvmModel
import libsvm

__all__ = [
    'LibSvmEpsilonRegressionModel',
    'LibSvmNuRegressionModel',
    'LibSvmRegressionResults'
    ]

# XXX document why get_svr_probability could be useful

class LibSvmRegressionResults:
    def __init__(self, model, traindataset, kernel, PredictorType):
        modelc = model.contents
        if modelc.param.svm_type not in [libsvm.EPSILON_SVR, libsvm.NU_SVR]:
            raise TypeError, '%s is for regression problems' % \
                str(self.__class__)
        self.rho = modelc.rho[0]
        self.sv_coef = modelc.sv_coef[0][:modelc.l]
        if modelc.probA:
            self.sigma = modelc.probA[0]
        self.predictor = PredictorType(model, traindataset, kernel)

    def predict(self, dataset):
        """
        This function does regression on a test vector x and returns
        the function value of x calculated using the model.
        """
        return [self.predictor.predict(x) for x in dataset]

    def get_svr_probability(self):
        """
        This function returns a value sigma > 0, which is a parameter
        in the probability model for the test data.

        For test data, we consider the probability model:

        target value = predicted value + z

        where z is the Laplace distribution: e^(-|x|/sigma)/(2*sigma).
        """
        return self.sigma

class LibSvmRegressionModel(LibSvmModel):
    ResultsType = LibSvmRegressionResults

    def __init__(self, kernel, **kwargs):
        LibSvmModel.__init__(self, kernel, **kwargs)

    def cross_validate(self, dataset, nr_fold):
        """
        Perform stratified cross-validation to determine the
        suitability of chosen model parameters.

        Data are separated to nr_fold folds. Each fold is validated
        against a model trained using the data from the remaining
        (nr_fold-1) folds.

        This function returns a 2-tuple containing the mean squared
        error and the squared correlation coefficient.
        """

        problem = dataset._create_svm_problem()
        target = N.empty((len(dataset.data),), dtype=N.float64)
        tp = cast(target.ctypes.data, POINTER(c_double))
        libsvm.svm_cross_validation(problem, self.param, nr_fold, tp)

        total_error = sumv = sumy = sumvv = sumyy = sumvy = 0.
        for i in range(len(dataset.data)):
            v = target[i]
            y = dataset.data[i][0]
            sumv = sumv + v
            sumy = sumy + y
            sumvv = sumvv + v * v
            sumyy = sumyy + y * y
            sumvy = sumvy + v * y
            total_error = total_error + (v - y) * (v - y)

        # mean squared error
        mse = total_error / len(dataset.data)
        # squared correlation coefficient
        l = len(dataset.data)
        scc = ((l * sumvy - sumv * sumy) * (l * sumvy - sumv * sumy)) / \
            ((l * sumvv - sumv*sumv) * (l * sumyy - sumy * sumy))

        return mse, scc

class LibSvmEpsilonRegressionModel(LibSvmRegressionModel):
    """
    A model for epsilon-SV regression.

    See also:

    - Smola, Schoelkopf. A Tutorial on Support Vector Regression.
    - Gunn. Support Vector Machines for Classification and Regression.
    - Mueller, Vapnik. Using Support Vector Machines for Time Series
      Prediction.
    """

    def __init__(self, kernel, epsilon=0.1, cost=1.0,
                 probability=False, **kwargs):
        LibSvmRegressionModel.__init__(self, kernel, **kwargs)
        self.epsilon = epsilon
        self.cost = cost
        self.param.svm_type = libsvm.EPSILON_SVR
        self.param.p = epsilon
        self.param.C = cost
        self.param.probability = probability

class LibSvmNuRegressionModel(LibSvmRegressionModel):
    """
    A model for nu-SV regression.

    See also: Schoelkopf, et al. New Support Vector Algorithms.
    """

    def __init__(self, kernel,
                 nu=0.5, cost=1.0, probability=False, **kwargs):
        LibSvmRegressionModel.__init__(self, kernel, **kwargs)
        self.nu = nu
        self.cost = cost
        self.param.svm_type = libsvm.NU_SVR
        self.param.nu = nu
        self.param.C = cost
        self.param.probability = probability
