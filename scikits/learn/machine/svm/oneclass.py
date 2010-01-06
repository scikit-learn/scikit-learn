from model import Model
import libsvm

__all__ = [
    'OneClassModel',
    'OneClassResults'
    ]

class OneClassResults:
    def __init__(self, model, traindataset, kernel, PredictorType):
        modelc = model.contents
        if modelc.param.svm_type != libsvm.ONE_CLASS:
            raise TypeError, '%s is for one-class problems' % \
                str(self.__class__)
        self.rho = modelc.rho[0]
        self.sv_coef = modelc.sv_coef[0][:modelc.l]
        self.predictor = PredictorType(model, traindataset, kernel)

    def predict(self, dataset):
        """
        This function returns a list of boolean values which indicate
        whether the test vectors form part of the distribution.
        """
        return [self.predictor.predict(x) > 0 for x in dataset]

    def predict_values(self, dataset):
        """
        This function returns a list of floating point values which
        indicate whether the test vectors form part of the
        distribution.

        A positive value indicates that the test vector is part of the
        distribution, while a non-positive value indicates that is is
        not.
        """
        return [self.predictor.predict_values(x, 1) for x in dataset]

    def compact(self):
        self.predictor.compact()

class OneClassModel(Model):
    """
    A model for distribution estimation (one-class SVM).

    See also: Scholkopf, et al. Estimating the Support of a
    High-Dimensional Distribution.
    """

    ResultsType = OneClassResults

    def __init__(self, kernel, nu=0.5, **kwargs):
        """
        Parameters:

        - `nu`: XXX
        """
        Model.__init__(self, kernel, **kwargs)
        self.nu = nu
        self.param.svm_type = libsvm.ONE_CLASS
        self.param.nu = nu
