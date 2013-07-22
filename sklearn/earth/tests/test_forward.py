'''
Created on Feb 16, 2013

@author: jasonrudy
'''
from pyearth._forward import ForwardPasser
from pyearth._basis import Basis, ConstantBasisFunction, HingeBasisFunction, LinearBasisFunction
import numpy
import os
from nose.tools import assert_true, assert_equal
        
class TestForwardPasser(object):
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
        self.forwardPasser = ForwardPasser(self.X,self.y,numpy.ones(self.y.shape),max_terms = 1000, penalty=1)
        
    def test_orthonormal_update(self):
        numpy.set_printoptions(precision=4)
        m,n = self.X.shape
        B_orth = self.forwardPasser.get_B_orth()
        v = numpy.random.normal(size=m)
        for i in range(1,10):
            v_ = numpy.random.normal(size=m)
            B_orth[:,i] = 10*v_ + v
            v = v_
            self.forwardPasser.orthonormal_update(i)
            assert_true(numpy.max(numpy.abs(numpy.dot(B_orth[:,0:i+1].transpose(),B_orth[:,0:i+1]) - numpy.eye(i+1))) < .0000001)
        
    def test_run(self):
        self.forwardPasser.run()
        res = str(self.forwardPasser.get_basis()) + '\n' + str(self.forwardPasser.trace())
#        with open('forward_regress.txt','w') as fl:
#            fl.write(res)
        with open(os.path.join(os.path.dirname(__file__), 'forward_regress.txt'),'r') as fl:
            prev = fl.read()
        assert_equal(res,prev)

if __name__ == '__main__':
    import nose
    nose.run(argv=[__file__, '-s', '-v'])