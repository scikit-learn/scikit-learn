# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer and Brian Holt
#
# License: BSD Style.


import numpy as np
cimport numpy as np

cimport cython

cdef extern from "math.h":
    cdef extern double log(double x)
    cdef extern double pow(double base, double exponent)

cdef extern from "float.h":
    cdef extern double DBL_MAX

"""
 Classification entropy measures
 
    From Hastie et al. Elements of Statistical Learning, 2009.
         
    If a target is a classification outcome taking on values 0,1,...,K-1
    In node m, representing a region Rm with Nm observations, let
            
       pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)
          
    be the proportion of class k observations in node m   
"""  

cdef class Criterion:
    """Interface for splitting criteria, both regression and classification."""
    
    cdef void init(self, np.ndarray[np.float64_t, ndim=1] y, 
                   int *sorted_X_i):
        """Initialise the criterion class."""
        pass

    cdef int update(self, int a, int b, 
                     np.float64_t *y,
                     int *sorted_X_i):
        """
        Update the criteria for each value in interval [a,b) 
        where a and b are indices in `sorted_X_i`. 
        """
        pass

    cdef double eval(self):
        """Evaluate the criteria (aka the split error). """
        pass

cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification.

    Attributes
    ----------
    n_classes : int
        number of classes
    n_samples : int
        number of samples
    label_count_left : int*
        label counts for samples left of splitting point.
    label_count_right : int*
        label counts for samples right of splitting point.
    n_left : int
        number of samples left of splitting point.
    n_right : int
        number of samples right of splitting point.
    """
    cdef n_classes
    cdef n_samples
    cdef int* label_count_left
    cdef int* label_count_right
    cdef int n_left
    cdef int n_right
    
    def __init__(self, 
                 int n_classes, 
                 np.ndarray[np.int32_t, ndim=1] label_count_left,
                 np.ndarray[np.int32_t, ndim=1] label_count_right):
        self.n_classes = n_classes
        self.n_samples = 0
        self.n_left = 0
        self.n_right = 0
        self.label_count_left = <int*>label_count_left.data
        self.label_count_right = <int*>label_count_right.data

    cdef void init(self, np.ndarray[np.float64_t, ndim=1] y, 
                   int *sorted_X_i):
        """Initialise the criterion class."""
        
        self.n_samples = y.shape[0]
        self.n_left = 0
        self.n_right = 0
        
        cdef int c = 0      
        for c from 0 <= c < self.n_classes:
            self.label_count_left[c] = 0
            self.label_count_right[c] = 0
        
        cdef int j = 0
        for j from 0 <= j < self.n_samples:
            c = <int>(y[j])
            self.label_count_right[c] += 1
            self.n_right += 1
        
        """    
        print "ClassificationCriterion.init: "
        print "    label_count_left [",    
        for c from 0 <= c < self.n_classes:    
            print self.label_count_left[c],  
        print "] ", self.n_left          
        print "    label_count_right [",    
        for c from 0 <= c < self.n_classes:    
            print self.label_count_right[c],
        print "] ", self.n_right        
        """

    cdef int update(self, int a, int b, 
                     np.float64_t *y,
                     int *sorted_X_i):
        """
        Update the criteria for each value in interval [a,b) 
        where a and b are indices in `sorted_X_i`. 
        """
        cdef int c
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            s = sorted_X_i[idx]
            c = <int>(y[s])
            self.label_count_right[c] -= 1
            self.label_count_left[c] += 1
            self.n_right -= 1
            self.n_left += 1

        """  
        print "ClassificationCriterion.update: a = ", a, " b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_samples:    
            print sorted_X_i[c],  
        print "] "                  
        print "    y [",    
        for c from 0 <= c < self.n_samples:    
            print <int>y[sorted_X_i[c]],
        print "] "               

        print "    label_count_left [",    
        for c from 0 <= c < self.n_classes:    
            print self.label_count_left[c],  
        print "] ", self.n_left          
        print "    label_count_right [",    
        for c from 0 <= c < self.n_classes:    
            print self.label_count_right[c],
        print "] ", self.n_right   
        """   
            
        return self.n_left

    cdef double eval(self):
        pass


cdef class Gini(ClassificationCriterion):
    """Gini Index splitting criteria.
    
    Gini index = \sum_{k=0}^{K-1} pmk (1 - pmk)
               = 1 - \sum_{k=0}^{K-1} pmk ** 2
    """

    cdef double eval(self):
        """Returns Gini index of left branch + Gini index of right branch. """
        cdef double H_left = 1.0
        cdef double H_right = 1.0
        cdef int k
        cdef double e1, e2
        cdef double n_left = <double> self.n_left
        cdef double n_right = <double> self.n_right
        
        for k from 0 <= k < self.n_classes:
            if self.label_count_left[k] > 0:
                H_left -= (self.label_count_left[k] / n_left) * (self.label_count_left[k] / n_left)
            if self.label_count_right[k] > 0:
                H_right -= (self.label_count_right[k] / n_right) * (self.label_count_right[k] / n_right)

        e1 = (n_left / self.n_samples) * H_left
        e2 = (n_right / self.n_samples) * H_right
        return e1 + e2


cdef class Entropy(ClassificationCriterion):
    """Entropy splitting criteria.

    Cross Entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
    
    """

    cdef double eval(self):
        """Returns Entropy of left branch + Entropy index of right branch. """
        cdef double H_left = 0.0
        cdef double H_right = 0.0
        cdef int k
        cdef double e1, e2
        cdef double n_left = <double> self.n_left
        cdef double n_right = <double> self.n_right
        
        for k from 0 <= k < self.n_classes:
            if self.label_count_left[k] > 0:
                H_left -= (self.label_count_left[k] / n_left) * log(self.label_count_left[k] / n_left)
            if self.label_count_right[k] > 0:    
                H_right -= (self.label_count_right[k] / n_right) * log(self.label_count_right[k] / n_right)

        e1 = (n_left / self.n_samples) * H_left
        e2 = (n_right / self.n_samples) * H_right
        return e1 + e2


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression.

    Attributes
    ----------
    n_samples : int
        The number of samples
    sum_left : double
        The sum of the samples left of the split point.
    sum_right : double
        The sum of the samples right of the split.
    y : double*
        Pointer to the labels array    
    n_left : int
        number of samples left of split point.
    n_right : int
        number of samples right of split point.
    """

    cdef n_samples
    cdef double sum_left
    cdef double sum_right
    cdef double* y
    cdef int n_left
    cdef n_right
    
    def __init__(self, np.ndarray[np.float64_t, ndim=1] y):
        self.sum_left = 0.0
        self.sum_right = 0.0
        self.n_left = 0
        self.n_right = 0
        self.y = <double*> y.data
        
    cdef void init(self, np.ndarray[np.float64_t, ndim=1] y, 
                   int *sorted_X_i):
        """Initialise the criterion class."""
        
        self.n_samples = y.shape[0]
        self.n_left = 0
        self.n_right = 0
        self.sum_left = 0.0
        self.sum_right = 0.0

        cdef int j = 0
        for j from 0 <= j < self.n_samples:
            self.y[j] = y[sorted_X_i[j]]
            self.sum_right += self.y[j]
            self.n_right += 1 

        """
        print "RegressionCriterion.init: "
        print "    sum_left = ", self.sum_left , self.n_left 
        print "    sum_right = ", self.sum_right , self.n_right     
        """
        
    cdef int update(self, int a, int b, 
                     np.float64_t *y,
                     int *sorted_X_i):
        """
        Update the criteria for each value in interval [a,b) 
        where a and b are indices in `sorted_X_i`. 
        """
        cdef double val = 0.0
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            val = self.y[idx]
            self.sum_right -= val
            self.sum_left += val
            self.n_right -= 1
            self.n_left += 1

        """
        print "RegressionCriterion.update: a = ", a, " b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_samples:    
            print sorted_X_i[c],  
        print "] "                 
        print "    y [",    
        for c from 0 <= c < self.n_samples:    
            print y[sorted_X_i[c]],
        print "] "               

        print "    sum_left = ", self.sum_left , self.n_left 
        print "    sum_right = ", self.sum_right , self.n_right      
        """
         
        return self.n_left

    cdef double eval(self):
        pass


cdef class MSE(RegressionCriterion):

    cdef double eval(self):
        """             
        MSE =  \sum_i (y_i - c0)^2  / N
        
        """
        cdef double mean_left = 0
        cdef double mean_right = 0

        if self.n_left > 0:
            mean_left = self.sum_left / self.n_left
        if self.n_right > 0:
            mean_right = self.sum_right / self.n_right
            
        #print "MSE.eval: mean_left = ", mean_left
        #print "MSE.eval: mean_right = ", mean_right        

        cdef double variance_left = 0.0
        cdef double variance_right = 0.0
        cdef int j
        cdef double e1, e2

        for j from 0 <= j < self.n_samples:
            #print "y[",j,"] = ", self.y[j]
            if j < self.n_left:
                variance_left += (self.y[j] - mean_left) * (self.y[j] - mean_left)
            else: 
                variance_right += (self.y[j] - mean_right) * (self.y[j] - mean_right)

        #print "MSE.eval: variance_left = ", variance_left, 
        #print "MSE.eval: variance_right = ", variance_right   

        if self.n_left > 0:
            variance_left /= self.n_left
        if self.n_right > 0:
            variance_right /= self.n_right

        e1 = ((<double> self.n_left) / self.n_samples) * variance_left
        e2 = ((<double> self.n_right) / self.n_samples) * variance_right

        return e1 + e2


cdef int smallest_sample_larger_than(int sample_idx,
                                     np.float64_t *X_i,
                                     int *sorted_X_i,
                                     int n_samples):
    """Find the largest next sample.
    
    Find the index in the `X_i` array for sample
    who's feature `i` value is just about
    greater than those of the sample `sorted_X_i[sample_idx]`.

    Returns
    -------
    next_sample_idx : int
        The index of the next smallest sample in `sorted_X`
        with different feature value than `sample_idx` .
        I.e. `sorted_X_i[sample_idx] < sorted_X_i[next_sample_idx]`
        -1 if no such element exists.
    """
    
    cdef int idx = 0
    cdef np.float64_t threshold = -DBL_MAX
    if sample_idx > -1:
        threshold = X_i[sorted_X_i[sample_idx]]
    for idx from sample_idx < idx < n_samples:
        if X_i[sorted_X_i[idx]] > threshold:
            return idx
    return -1


def _find_best_split(np.ndarray[np.float64_t, ndim=2, mode="fortran"] X,
                     np.ndarray[np.float64_t, ndim=1, mode="c"] y,
                     Criterion criterion):
    """Find the best dimension and threshold that minimises the error.
    
    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features), dtype=np.float64
        The feature values.
    y : ndarray, shape (n_samples,), dtype=float
        The y.
    criterion : Criterion
        The criterion function to be minimized.

    Returns
    -------
    best_i : int
        The split feature or -1 if criterion not smaller than `parent_split_error`.
    best_t : np.float64_t
        The split threshold
    best_error : np.float64_t
        The error of the split.
    """

    cdef int n_samples = X.shape[0]
    cdef int n_features = X.shape[1]
    cdef int i, a, b, best_i = -1
    cdef np.float64_t t, initial_error, error
    cdef np.float64_t best_error = np.inf, best_t = np.inf

    # pointer access to ndarray data
    cdef double *y_ptr = <double*>y.data 
    cdef np.float64_t *X_i = NULL
    cdef np.ndarray[np.int32_t, ndim=2, mode="fortran"] sorted_X = \
        np.asfortranarray(np.argsort(X, axis=0).astype(np.int32))
    cdef int *sorted_X_i = NULL    
    
    # Compute the column strides (inc in pointer elem to get from col i to i+1)
    # for `X` and `sorted_X`
    cdef int X_elem_stride = X.strides[0]
    cdef int X_col_stride = X.strides[1]
    cdef int X_stride = X_col_stride / X_elem_stride
    cdef int sorted_X_elem_stride = sorted_X.strides[0]
    cdef int sorted_X_col_stride = sorted_X.strides[1]
    cdef int sorted_X_stride = sorted_X_col_stride / sorted_X_elem_stride   
    
    #compute the initial entropy in the node
    sorted_X_i = (<int *>sorted_X.data)
    criterion.init(y, sorted_X_i)
    initial_error = criterion.eval()
    if initial_error == 0: # break early if the node is pure
        return best_i, best_t, best_error, initial_error 
    best_error = initial_error
    #print 'at init, best error = ', best_error
    
    for i from 0 <= i < n_features:
        # get i-th col of X and X_sorted
        X_i = (<np.float64_t *>X.data) + X_stride * i
        sorted_X_i = (<int *>sorted_X.data) + sorted_X_stride * i

        # init the criterion for this feature
        criterion.init(y, sorted_X_i)

        # get sample in mask with smallest value for i-th feature
        a = smallest_sample_larger_than(-1, X_i, sorted_X_i, n_samples)
        
        while True:
            b = smallest_sample_larger_than(a, X_i, sorted_X_i, n_samples)
            
            # if -1 there's none and we are finished
            if b == -1:
                break

            criterion.update(a, b, y_ptr, sorted_X_i)

            # get criterion value
            error = criterion.eval()
            #print 'error = ', error
            
            # check if current error is smaller than previous best
            # if this is never true best_i is -1.
            if error < best_error:
                t = (X_i[sorted_X_i[a]] + X_i[sorted_X_i[b]]) / 2.0
                best_i = i
                best_t = t
                best_error = error

            a = b
                
    return best_i, best_t, best_error, initial_error
