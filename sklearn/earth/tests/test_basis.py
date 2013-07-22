'''
Created on Feb 17, 2013

@author: jasonrudy
'''
from nose.tools import assert_true, assert_false, assert_equal
from pyearth._basis import Basis, ConstantBasisFunction, HingeBasisFunction, LinearBasisFunction
import numpy
import pickle
import os

class BaseTestClass(object):
    def __init__(self):
        numpy.random.seed(0)
        data = numpy.genfromtxt(os.path.join(os.path.dirname(__file__),'test_data.csv'),
                                 delimiter=',', skip_header=1)
        self.y = numpy.array(data[:,5])
        self.X = numpy.array(data[:,0:5])

class TestConstantBasisFunction(BaseTestClass):
        
    def __init__(self):
        super(self.__class__,self).__init__()
        self.bf = ConstantBasisFunction()
        
    def test_apply(self):
        m,n = self.X.shape
        B = numpy.empty(shape=(m,10))
        
        assert_false(numpy.all(B[:,0] == 1))
        self.bf.apply(self.X,B[:,0])
        assert_true(numpy.all(B[:,0] == 1))
        
    def test_pickle_compatibility(self):
        bf_copy = pickle.loads(pickle.dumps(self.bf))
        assert_true(self.bf == bf_copy)

class TestHingeBasisFunction(BaseTestClass):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.parent = ConstantBasisFunction()
        self.bf = HingeBasisFunction(self.parent,1.0,10,1,False)
        
    def test_getters(self):
        assert self.bf.get_reverse() == False
        assert self.bf.get_knot() == 1.0
        assert self.bf.get_variable() == 1
        assert self.bf.get_knot_idx() == 10
        assert self.bf.get_parent() == self.parent
        
    def test_apply(self):
        m,n = self.X.shape
        B = numpy.ones(shape=(m,10))
        self.bf.apply(self.X,B[:,0])
        assert_true(numpy.all(B[:,0] == (self.X[:,1] - 1.0) * (self.X[:,1] > 1.0)))
    
    def test_degree(self):
        assert_equal(self.bf.degree(),1)
        
    def test_pickle_compatibility(self):
        bf_copy = pickle.loads(pickle.dumps(self.bf))
        assert_true(self.bf == bf_copy)
        
class TestLinearBasisFunction(BaseTestClass):
    def __init__(self):
        super(self.__class__,self).__init__()
        parent = ConstantBasisFunction()
        self.bf = LinearBasisFunction(parent,1)
        
    def test_apply(self):
        m,n = self.X.shape
        B = numpy.ones(shape=(m,10))
        self.bf.apply(self.X,B[:,0])
        assert_true(numpy.all(B[:,0] == self.X[:,1]))
    
    def test_degree(self):
        assert_equal(self.bf.degree(),1)
        
    def test_pickle_compatibility(self):
        bf_copy = pickle.loads(pickle.dumps(self.bf))
        assert_true(self.bf == bf_copy)
    
class TestBasis(BaseTestClass):
    def __init__(self):
        super(self.__class__,self).__init__()
        self.basis = Basis()
        self.parent = ConstantBasisFunction()
        self.bf = HingeBasisFunction(self.parent,1.0,10,1,False)
        self.basis.append(self.parent)
        self.basis.append(self.bf)
        
    def test_add(self):
        assert_equal(len(self.basis),2)
    
    def test_translate_and_scale(self):
        m,n = self.X.shape
        numpy.random.seed(1)
        B = numpy.empty(shape=(m,self.basis.plen()))
        self.basis.transform(self.X,B)
        B_ = numpy.empty(shape=(m,self.basis.plen()))
        mu = numpy.mean(self.X,axis=0)
        sigma = numpy.std(self.X,axis=0)
        coeff = numpy.random.normal(size=B.shape[1])
        X_ = self.X * sigma + mu
        coeff_ = coeff.copy()
        self.basis.translate(sigma,mu)
        self.basis.scale(sigma,mu,coeff_)
        self.basis.transform(X_,B_)
        assert_true(numpy.all((numpy.dot(B,coeff) - numpy.dot(B_,coeff_))**2 < 1e-12))
    
    def test_pickle_compat(self):
        basis_copy = pickle.loads(pickle.dumps(self.basis))
        assert_true(self.basis == basis_copy)

if __name__ == '__main__':
    import nose
    nose.run(argv=[__file__, '-s', '-v'])
