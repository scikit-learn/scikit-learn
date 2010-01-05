from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.dataset import *
from svm.kernel import LinearKernel
from svm.predict import *
from svm.regression import *
restore_path()

class test_precomputed(NumpyTestCase):
    def check_precomputed(self):
        ModelType = LibSvmEpsilonRegressionModel
        ResultType = LibSvmRegressionResults

        kernel = LinearKernel()

        # this dataset remains constant
        y1 = N.random.randn(50)
        x1 = N.random.randn(len(y1), 10)
        data1 = LibSvmRegressionDataSet(zip(y1, x1))
        pcdata1 = data1.precompute(kernel)

        # in a typical problem, this dataset would be smaller than the
        # part that remains constant and would differ for each model
        y2 = N.random.randn(5)
        x2 = N.random.randn(len(y2), x1.shape[1])
        data2 = LibSvmRegressionDataSet(zip(y2, x2))

        pcdata12 = pcdata1.combine(data2)
        model = LibSvmEpsilonRegressionModel(kernel)
        results = model.fit(pcdata12, ResultType, LibSvmPredictor)

        # reference model, calculated without involving the
        # precomputed Gram matrix
        refy = N.concatenate([y1, y2])
        refx = N.vstack([x1, x2])
        refdata = LibSvmRegressionDataSet(zip(refy, refx))
        model = ModelType(kernel)
        #refresults = model.fit(refdata, ResultType,
        #                       LibSvmPrecomputedPredictor)

        #self.assertAlmostEqual(results.rho, refresults.rho)
        #assert_array_almost_equal(results.sv_coef, refresults.sv_coef)

        # XXX sigmas don't match yet. need to find out why.
        #self.assertAlmostEqual(results.sigma, refresults.sigma)

if __name__ == '__main__':
    NumpyTest().run()
