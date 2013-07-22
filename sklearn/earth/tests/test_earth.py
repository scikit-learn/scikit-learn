'''
Created on Feb 24, 2013

@author: jasonrudy
'''
import numpy
from pyearth._basis import Basis, ConstantBasisFunction, HingeBasisFunction, LinearBasisFunction
from pyearth import Earth
import pickle
import copy
import os
from testing_utils import if_statsmodels, if_pandas, if_patsy
from nose.tools import assert_equal, assert_not_equal, assert_true, assert_false, \
    assert_almost_equal, assert_list_equal
    
class TestEarth(object):

    def __init__(self):
        numpy.random.seed(0)
        self.basis = Basis()
        constant = ConstantBasisFunction()
        self.basis.append(constant)
        bf1 = HingeBasisFunction(constant,0.1,10,1,False,'x1')
        bf2 = HingeBasisFunction(constant,0.1,10,1,True,'x1')
        bf3 = LinearBasisFunction(bf1,2,'x2')
        self.basis.append(bf1)
        self.basis.append(bf2)
        self.basis.append(bf3)
        self.X = numpy.random.normal(size=(100,10))
        self.B = numpy.empty(shape=(100,4),dtype=numpy.float64)
        self.basis.transform(self.X,self.B)
        self.beta = numpy.random.normal(size=4)
        self.y = numpy.empty(shape=100,dtype=numpy.float64)
        self.y[:] = numpy.dot(self.B,self.beta) + numpy.random.normal(size=100)
        self.earth = Earth(penalty=1)
        
    def test_get_params(self):
        assert_equal(Earth().get_params(), {'penalty': None, 'min_search_points': None, 
                                            'endspan_alpha': None, 'check_every': None, 
                                            'max_terms': None, 'xlabels': None, 'max_degree': None, 
                                            'minspan_alpha': None, 'linvars': None, 'thresh': None, 
                                            'minspan': None, 'endspan': None})
        assert_equal(Earth(max_degree=3).get_params(), {'penalty': None, 'min_search_points': None, 
                                            'endspan_alpha': None, 'check_every': None, 
                                            'max_terms': None, 'xlabels': None, 'max_degree': 3, 
                                            'minspan_alpha': None, 'linvars': None, 'thresh': None, 
                                            'minspan': None, 'endspan': None})
    
    @if_statsmodels
    def test_linear_fit(self):
        from statsmodels.regression.linear_model import GLS, OLS
        self.earth.fit(self.X, self.y)
        self.earth.linear_fit(self.X, self.y)
        soln = OLS(self.y, self.earth.transform(self.X)).fit().params
        assert_almost_equal(numpy.mean((self.earth.coef_-soln)**2), 0.0)
        
        weights = 1.0 / (numpy.random.normal(size=self.y.shape) ** 2)
        self.earth.fit(self.X, self.y)
        self.earth.linear_fit(self.X, self.y, weights)
        soln = GLS(self.y, self.earth.transform(self.X), 1.0 / weights).fit().params
        assert_almost_equal(numpy.mean((self.earth.coef_-soln)**2), 0.0)
        
    def test_weights(self):
        group = numpy.random.binomial(1,.5,size=1000)  == 1
        weights =  1 / (group * 100 + 1.0)
        x = numpy.random.uniform(-10,10,size=1000)
        y = numpy.abs(x)
        y[group] = numpy.abs(x[group] - 5)
        y += numpy.random.normal(0,1,size=1000)
        model = Earth().fit(x,y,weights = weights)
        
        #Check that the model fits better for the more heavily weighted group
        assert_true(model.score(x[group],y[group]) < model.score(x[numpy.logical_not(group)],y[numpy.logical_not(group)]))
        
        #Make sure that the score function gives the same answer as the trace
        assert_almost_equal(model.score(x,y,weights=weights), model.pruning_trace().rsq(model.pruning_trace().get_selected()))
        
        #Uncomment below to see what this test situation looks like
#        from matplotlib import pyplot
#        print model.summary()
#        print model.score(x,y,weights = weights)
#        pyplot.figure()
#        pyplot.plot(x,y,'b.')
#        pyplot.plot(x,model.predict(x),'r.')
#        pyplot.show()
        
    def test_fit(self):
        self.earth.fit(self.X, self.y)
        res = str(self.earth.trace()) + '\n' + self.earth.summary()
#        with open('earth_regress.txt','w') as fl:
#            fl.write(res)
        with open(os.path.join(os.path.dirname(__file__),'earth_regress.txt'),'r') as fl:
            prev = fl.read()
        assert_equal(res,prev)
        
    def test_score(self):
        model = self.earth.fit(self.X, self.y)
        record = model.pruning_trace()
        rsq = record.rsq(record.get_selected())
        assert_almost_equal(rsq,model.score(self.X,self.y))

    @if_pandas
    def test_pandas_compatibility(self):
        import pandas
        X = pandas.DataFrame(self.X)
        y = pandas.DataFrame(self.y)
        colnames = ['xx'+str(i) for i in range(X.shape[1])]
        X.columns = colnames
        model = self.earth.fit(X,y)
        assert_list_equal(colnames,model.xlabels)
        
    @if_patsy
    @if_pandas
    def test_patsy_compatibility(self):
        import pandas
        import patsy
        X = pandas.DataFrame(self.X)
        y = pandas.DataFrame(self.y)
        colnames = ['xx'+str(i) for i in range(X.shape[1])]
        X.columns = colnames
        X['y'] = y
        y, X = patsy.dmatrices('y ~ xx0 + xx1 + xx2 + xx3 + xx4 + xx5 + xx6 + xx7 + xx8 + xx9 - 1',data=X)
        model = self.earth.fit(X,y)
        assert_list_equal(colnames,model.xlabels)
        
    def test_pickle_compatibility(self):
        model = self.earth.fit(self.X, self.y)
        model_copy = pickle.loads(pickle.dumps(model))
        assert_true(model_copy == model)
        assert_true(numpy.all(model.predict(self.X) == model_copy.predict(self.X)))
        assert_true(model.basis_[0] is model.basis_[1]._get_root())
        assert_true(model_copy.basis_[0] is model_copy.basis_[1]._get_root())
        
    def test_copy_compatibility(self):
        model = self.earth.fit(self.X, self.y)
        model_copy = copy.copy(model)
        assert_true(model_copy == model)
        assert_true(numpy.all(model.predict(self.X) == model_copy.predict(self.X)))
        assert_true(model.basis_[0] is model.basis_[1]._get_root())
        assert_true(model_copy.basis_[0] is model_copy.basis_[1]._get_root())
        
if __name__ == '__main__':
    import nose
    nose.run(argv=[__file__, '-s', '-v'])
