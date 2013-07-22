'''
This script randomly generates earth-style models, then randomly generates data from those models and 
fits earth models to those data using both the python and R implementations.  It records the sample size,
m, the number of input dimensions, n, the number of forward pass iterations, the runtime, and the r^2 
statistic for each fit and writes the result to a CSV file.
'''

import numpy
import pandas.rpy.common as com
import rpy2.robjects as robjects
import time
import pandas
from sklearn.earth import Earth

class DataGenerator(object):
    def __init__(self):
        pass
    def generate(self, m):
        pass


class NoiseGenerator(DataGenerator):
    def __init__(self, n):
        self.n = n
        
    def generate(self, m):
        X = numpy.random.normal(size=(m,self.n))
        y = numpy.random.normal(size=m)
        return X, y
    
class LinearGenerator(DataGenerator):
    def __init__(self, n):
        self.n = n
        
    def generate(self, m):
        X = numpy.random.normal(size=(m,self.n))
        beta = numpy.random.normal(size=self.n)
        y = numpy.dot(X,beta) + numpy.random.normal(m)
        return X, y

class VFunctionGenerator(DataGenerator):
    def __init__(self, n):
        self.n = n
        
    def generate(self, m):
        X = numpy.random.normal(size=(m,self.n))
        var = numpy.random.randint(self.n)
        y = 10*abs(X[:,var]) + numpy.random.normal(m)
        return X, y
    
class UFunctionGenerator(DataGenerator):
    def __init__(self, n):
        self.n = n
        
    def generate(self, m):
        X = numpy.random.normal(size=(m,self.n))
        var = numpy.random.randint(self.n)
        y = 10*(X[:,var]**2) + numpy.random.normal(m)
        return X, y
        
class RandomComplexityGenerator(DataGenerator):
    def __init__(self, n, max_terms=10, max_degree=2):
        self.n = n
        self.max_terms = max_terms
        self.max_degree = max_degree
        
        
    def generate(self, m):
        X = numpy.random.normal(size=(m,self.n))
        num_terms = numpy.random.randint(2,self.max_terms) #Including the intercept
        coef = 10*numpy.random.normal(size=num_terms)
        B = numpy.ones(shape=(m,num_terms))
        B[:,0] += coef[0]
        for i in range(1,num_terms):
            degree = numpy.random.randint(1,self.max_degree)
            for bf in range(degree):
                knot = numpy.random.normal()
                dir = 1 - 2*numpy.random.binomial(1,.5)
                var = numpy.random.randint(0,self.n)
                B[:,i] *= (dir*(X[:,var] - knot)) * (dir*(X[:,var] - knot) > 0)
        y = numpy.dot(B,coef) + numpy.random.normal(size=m)
        return X, y
        
        
def run_earth(X, y, **kwargs):
    '''Run with the R package earth.  Return prediction value, training time, and number of forward pass iterations.'''
    r = robjects.r
    m,n = X.shape
    data = pandas.DataFrame(X)
    data['y'] = y
    r_data = com.convert_to_r_dataframe(data)
    r('library(earth)')
    r_func = '''
        run <-  function(data, degree=1, fast.k=0, penalty=3.0){
                    time = system.time(model <- earth(y~.,data=data,degree=degree,penalty=penalty))[3]
                    forward_terms = dim(summary(model)$prune.terms)[1]
                    y_pred = predict(model,data)
                    return(list(y_pred, time, forward_terms, model))
                }
        '''
    r(r_func)
    run = r('run')
    r_list = run(**{'data':r_data,'degree':kwargs['max_degree'],'fast.k':0,'penalty':kwargs['penalty']})
    y_pred = numpy.array(r_list[0]).reshape(m)
    time = r_list[1][0]
    forward_terms = r_list[2][0]
    return y_pred, time, (forward_terms - 1) / 2
    
def run_pyearth(X, y, **kwargs):
    '''Run with pyearth.  Return prediction value, training time, and number of forward pass iterations.'''
    model = Earth(**kwargs)
    t0 = time.time()
    model.fit(X,y)
    t1 = time.time()
    y_pred = model.predict(X)
    forward_iterations = len(model.forward_trace()) - 1
    return y_pred, t1-t0, forward_iterations
    
def compare(generator_class, sample_sizes, dimensions, repetitions, **kwargs):
    '''Return a data table that includes m, n, pyearth or earth, training time, and number of forward pass iterations.'''
    header = ['m','n','pyearth','earth','time','forward_iterations','rsq']
    data = []
    for n in dimensions:
        generator = generator_class(n=n)
        for m in sample_sizes:
            for rep in range(repetitions):
                print n, m, rep
                X, y = generator.generate(m=m)
                y_pred_r, time_r, iter_r = run_earth(X,y,**kwargs)
                rsq_r = 1 - (numpy.sum((y-y_pred_r)**2))/(numpy.sum((y-numpy.mean(y))**2))
                data.append([m,n,0,1,time_r,iter_r,rsq_r])
                y_pred_py, time_py, iter_py = run_pyearth(X,y,**kwargs)
                rsq_py = 1 - (numpy.sum((y-y_pred_py)**2))/(numpy.sum((y-numpy.mean(y))**2))
                data.append([m,n,1,0,time_py,iter_py,rsq_py])
    return pandas.DataFrame(data,columns=header)

if __name__ == '__main__':
    sample_sizes = [100, 200, 300, 500]
    dimensions = [10, 20, 30]
    rep = 5
    numpy.random.seed(1)
    data = compare(RandomComplexityGenerator,sample_sizes,dimensions,rep,max_degree=2,penalty=3.0)
    print data
    data.to_csv('comparison.csv')

        