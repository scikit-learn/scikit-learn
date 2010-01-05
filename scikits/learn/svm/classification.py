__all__ = [
    'LibSvmCClassificationModel',
    'LibSvmNuClassificationModel',
    ]

import numpy as N
from ctypes import cast, POINTER, c_int, c_double

from model import LibSvmModel
import libsvm

class LibSvmClassificationResults:
    def __init__(self, model, dataset):
        self.model = model
        # keep a reference to the dataset here because the support
        # vectors point to some or all of the vectors in the dataset
        self.dataset = dataset
        # XXX this is probably suboptimal when training many models
        # that each only have a few support vectors. look at some
        # options when we do memory optimization.

        model = model.contents
        self.nr_class = model.nr_class
        self.labels = model.labels[:self.nr_class]
        self.rho = model.rho[:self.nr_class*(self.nr_class-1)/2]
        self.nSV = model.nSV[:self.nr_class]
        sv_coef = N.empty((self.nr_class-1, model.l), dtype=N.float64)
        for i, c in enumerate(model.sv_coef[:self.nr_class-1]):
            sv_coef[i,:] = c[:model.l]
        self.sv_coef = sv_coef

    def __del__(self):
        libsvm.svm_destroy_model(self.model)

    def predict(self, dataset):
        """
        This function does classification on a test vector x and
        returns the label of the predicted class.
        """
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            return int(libsvm.svm_predict(self.model, xptr))
        return map(p, dataset.data)

    def predict_values(self, dataset):
        """
        This function does classification on a test vector x and
        returns decision values.

        For training data with nr_class classes, this function returns
        nr_class*(nr_class-1)/2 decision values in a dictionary. The
        keys of the dictionary are 2-tuples, one for each combination
        of two class labels.
        """
        def p(x):
            xptr = cast(x.ctypes.data, POINTER(libsvm.svm_node))
            n = self.nr_class*(self.nr_class-1)/2
            v = N.empty((n,), dtype=N.float64)
            vptr = cast(v.ctypes.data, POINTER(c_double))
            libsvm.svm_predict_values(self.model, xptr, vptr)
            count = 0
            d = {}
            for i in range(len(self.labels)):
                for j in range(i+1, len(self.labels)):
                    d[self.labels[i], self.labels[j]] = v[count]
                    d[self.labels[j], self.labels[i]] = -v[count]
                    count += 1
            return d
        return map(p, dataset.data)

    def predict_probability(self, x):
        """
        This function does classification on a test vector x for a
        model with probability information.

        This function returns a 2-tuple. The first item is the label
        of the class with the highest probability. The second item is
        a list of all the class probabilities.
        """
        raise NotImplementedError

class LibSvmClassificationModel(LibSvmModel):
    """
    A model for support vector classification.

    Classification models can predict a class label, decision values
    over all classes or a posterior class probability.

    See also:

    - Platt. Probabilistic Outputs for Support Vector Machines and
      Comparisons to Regularized Likelihood Methods.
    - Lin. A Note on Platt's Probabilistic Outputs for Support Vector
      Machines.
    """

    Results = LibSvmClassificationResults

    def __init__(self, kernel, weights, **kwargs):
        LibSvmModel.__init__(self, kernel, **kwargs)
        if weights is not None:
            # XXX check whether labels need to be sorted
            self.weight_labels = N.empty((len(weights),), dtype=N.intp)
            self.weights = N.empty((len(weights),), dtype=N.float64)
            for i, (label, weight) in enumerate(weights):
                self.weight_labels[i] = label
                self.weights[i] = weight

            self.param.nr_weight = len(weights)
            self.param.weight_label = \
                cast(self.weight_labels.ctypes.data, POINTER(c_int))
            self.param.weight = \
                cast(self.weights.ctypes.data, POINTER(c_double))

class LibSvmCClassificationModel(LibSvmClassificationModel):
    """
    A model for C-SV classification.

    See also:

    - Hsu, et al. A Practical Guide to Support Vector Classification.
    - Gunn. Support Vector Machines for Classification and Regression.
    - Burges. A Tutorial on Support Vector Machines for Pattern
      Recognition.
    """

    def __init__(self, kernel, cost=1.0, weights=None, **kwargs):
        """
        Parameters:

        - `cost`: XXX
        - `weights`: XXX
        """
        LibSvmClassificationModel.__init__(self, kernel, weights, **kwargs)
        self.cost = cost
        self.param.svm_type = libsvm.C_SVC
        self.param.C = cost
        # always train probability model parameters
        # XXX review this decision when we do performance optimization
        self.param.probability = True

class LibSvmNuClassificationModel(LibSvmClassificationModel):
    """
    A model for nu-SV classification.

    See also:

    - Chen, et al. A Tutorial on nu-Support Vector Machines.
    - Scholkopf, et al. New Support Vector Algorithms.
    """

    def __init__(self, kernel, nu=0.5, weights=None, **kwargs):
        """
        Parameters:

        - `nu`: XXX
        - `weights`: XXX
        """
        LibSvmClassificationModel.__init__(self, kernel, weights, **kwargs)
        self.nu = nu
        self.param.svm_type = libsvm.NU_SVC
        self.param.nu = nu
        # always train probability model parameters
        # XXX review this decision when we do performance optimization
        self.param.probability = True
