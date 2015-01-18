##############################
# Author: B. Hsu 
# Licence: BSD 3 claused
#
#
###############################
import numpy as np
cimport numpy as np

cimport cython

@cython.cdivision(True)
cpdef np.ndarray[double, ndim=2] wcov(np.ndarray[double, ndim=2] data, np.ndarray[double, ndim=1] weights):
    """ 
    Function to compute weighted covariance matrix correctly. 
    """
    w = np.diagflat( weights/np.sum(weights), k=0)
    data = np.sqrt(w).dot(data - np.sum( w.dot(data),axis=0))

    return 1./(1.-np.sum(w**2)) * data.T.dot(data)

@cython.cdivision(True)
cpdef np.ndarray[double, ndim=2] wcorr(np.ndarray[double, ndim=2] data, np.ndarray[double, ndim=1] weights):
    """
    weighted correlation function
    """
    std = np.diagflat( np.sqrt( wcov(data, weights).diagonal())**-1 )
    return std.dot(wcov(data, weights) ).dot(std)

cdef int factorial(int x):
    if x<=1:
        return 1
    else:
        return x * factorial(x-1)
    
cdef int NVars(int n):
    """
    Count the number of non-zero bits in the integer n. These correspond to the columns/variables
    included in the model.
    """
    cdef int count = 0
    while n > 0:
        count += 1
        n = n & (n-1)
    return count

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[int, ndim=1] asMask(np.ndarray[np.uint8_t,ndim=1] A):
    cdef int length = len(A)
    cdef int i
    cdef list B = []
    
    for i in xrange(length):
        if A[i] != 0: B.append(i)
            
    return np.array(B)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double Rsq(np.ndarray[double, ndim=2] S, int i, int n):
    """
    Computs the coefficient of determinantion (R-squared) from the covariance matrix S for a model 
    containing variables i out of n total covariates. 
    """
    cdef np.ndarray[double, ndim=2] Sxx, Sxx_inv, 
    cdef np.ndarray[double, ndim=1] Syx, Sxy
    cdef np.ndarray[np.uint8_t, ndim=1] cov = np.zeros(n+1, dtype=np.uint8)
    cdef int k 

    if i == 0: return 0.0
    
    for k in xrange(n):
        if i & (1<<k)!=0: cov[k+1] = 1

    Sxx = S[:, asMask(cov)]
    Sxx_inv = np.linalg.inv( Sxx[asMask(cov),:] )
    Syx = S[0, asMask(cov)]; Sxy = S[asMask(cov),0]
    return Syx.dot( Sxx_inv.dot(Sxy) )

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[double, ndim=1] ShapleyValue( np.ndarray[double, ndim=2] S ):
    """
    Computes the Shapley importance of each n covariate in a linear model from the (weighted) covariance
    (n+1) x (n+1) matrix, S. Returns a vector of length n giving the average amount of R2 attributed 
    to the nth variable. 
    """
    cdef unsigned int n_cov, k, ij 
    n_cov = S.shape[1]-1
    cdef np.ndarray[double,ndim=1] model = np.zeros(2**n_cov, dtype = np.double) ### storage space to memoize the models computed
    cdef np.ndarray[double, ndim=1] shapley = np.zeros(n_cov, dtype = np.double)   ### storage space for the shapley importances
    
    
    model[0] = 0.  ### no covariates in model, base case, R2 for no covariates
    for i in xrange(2**n_cov):
        
        k = NVars(i)
        if model[i] == 0: model[i] = Rsq(S,i,n_cov)
        
        for ij in xrange(n_cov):
            ### add the ij variable to the i. if its already in i, continue, else compute new model
            j = (1<<ij)  
            if i == i|j: continue
                
            if model[i|j] == 0:  model[i|j] = Rsq(S,i|j,n_cov) 

            ### compute the improvement in R2 given the addition of the jth variable and average over
            ### permutations possible
            shapley[ ij ] += factorial(k) *1. * factorial(n_cov - k -1)/factorial(n_cov)* (model[i|j]-model[i])/S[0,0]  
    return shapley