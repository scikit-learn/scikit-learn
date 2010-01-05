from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.classification import *
from svm.dataset import *
from svm.kernel import *
from svm.predict import *
from svm.regression import *
restore_path()

class test_precomputed(NumpyTestCase):
    def _create_regression_data(self):
        y1 = N.random.randn(50)
        x1 = N.random.randn(len(y1), 10)
        data1 = LibSvmRegressionDataSet(zip(y1, x1))
        pcdata1 = data1.precompute(kernel)
        y2 = N.random.randn(5)
        x2 = N.random.randn(len(y2), x1.shape[1])
        data2 = LibSvmRegressionDataSet(zip(y2, x2))
        pcdata12 = pcdata1.combine(data2)
        refy = N.concatenate([y1, y2])
        refx = N.vstack([x1, x2])
        testdata = LibSvmTestDataSet(refx)

    def xcheck_precomputed_classification(self):
        ModelType = LibSvmCClassificationModel
        ResultType = LibSvmClassificationResults
        kernel = LinearKernel()

        labels1 = ([0] * 10) + ([1] * 10) + ([2] * 10)
        x1 = N.random.randn(len(labels1), 10)
        data1 = LibSvmClassificationDataSet(zip(labels1, x1))
        pcdata1 = data1.precompute(kernel)

        labels2 = ([0] * 5) + ([1] * 5) + ([2] * 5)
        x2 = N.random.randn(len(labels2), x1.shape[1])
        data2 = LibSvmClassificationDataSet(zip(labels2, x2))

        pcdata12 = pcdata1.combine(data2)
        model = LibSvmCClassificationModel(kernel)
        results = model.fit(pcdata12, ResultType,LibSvmPrecomputedPredictor)

        reflabels = labels1 + labels2
        refx = N.vstack([x1, x2])
        refdata = LibSvmClassificationDataSet(zip(reflabels, refx))
        model = ModelType(kernel)
        refresults = model.fit(refdata, ResultType, LibSvmPredictor)

        assert_array_almost_equal(results.rho, refresults.rho)
        assert_array_almost_equal(results.sv_coef, refresults.sv_coef)

        testdata = LibSvmTestDataSet(refx)
        p = results.predict(testdata)
        refp = refresults.predict(testdata)
        assert_array_almost_equal(p, refp)

        pv = results.predict_values(testdata)
        refpv = refresults.predict_values(testdata)
        for v, refv in zip(pv, refpv):
            for key, value in refv.iteritems():
                self.assertAlmostEqual(v[key], value)

        pp = results.predict_probability(testdata)
        refpp = refresults.predict_probability(testdata)
        for (lbl, p), (reflbl, refp) in zip(pp, refpp):
            self.assertEqual(lbl, reflbl)
            assert_array_almost_equal(p, refp)

    def check_precomputed_regression(self):
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
        results = model.fit(pcdata12, ResultType, LibSvmPrecomputedPredictor)

        # reference model, calculated without involving the
        # precomputed Gram matrix
        refy = N.concatenate([y1, y2])
        refx = N.vstack([x1, x2])
        refdata = LibSvmRegressionDataSet(zip(refy, refx))
        model = ModelType(kernel)
        refresults = model.fit(refdata, ResultType, LibSvmPredictor)

        self.assertAlmostEqual(results.rho, refresults.rho)
        assert_array_almost_equal(results.sv_coef, refresults.sv_coef)

        # XXX sigmas don't match exactly. need to find out why.
        #self.assertAlmostEqual(results.sigma, refresults.sigma)

        testdata = LibSvmTestDataSet(refx)
        p = results.predict(testdata)
        refp = refresults.predict(testdata)
        assert_array_almost_equal(p, refp)

    def check_precomputed_regression_py(self):
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
        results = model.fit(pcdata12, ResultType, LibSvmPythonPredictor)

        # reference model, calculated without involving the
        # precomputed Gram matrix
        refy = N.concatenate([y1, y2])
        refx = N.vstack([x1, x2])
        refdata = LibSvmRegressionDataSet(zip(refy, refx))
        model = ModelType(kernel)
        refresults = model.fit(refdata, ResultType, LibSvmPredictor)

        self.assertAlmostEqual(results.rho, refresults.rho)
        assert_array_almost_equal(results.sv_coef, refresults.sv_coef)

        # XXX sigmas don't match exactly. need to find out why.
        #self.assertAlmostEqual(results.sigma, refresults.sigma)

        testdata = LibSvmTestDataSet(refx)
        p = results.predict(testdata)
        refp = refresults.predict(testdata)
        assert_array_almost_equal(p, refp)

if __name__ == '__main__':
    NumpyTest().run()
