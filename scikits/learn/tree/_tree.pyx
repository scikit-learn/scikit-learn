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
    """Interface for splitting criteria, both regression and classification.
    """
    
    cdef void init(self, np.ndarray[np.float_t, ndim=1] labels, 
                   int *sorted_features_i):
        """Init the criteria for each feature `i`. """
        pass

    cdef int update(self, int a, int b, 
                     np.float64_t *labels,
                     int *sorted_features_i):
        """
        Update the criteria for each value in interval [a,b) 
        where a and b are indices in `sorted_features_i`. 
        """
        pass

    cdef double eval(self):
        """Evaluate the criteria (aka the split error). """
        pass

cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification.

    Attributes
    ----------
    pm_left_ptr : np.ndarray
        label counts for samples left of splitting point.
    pm_right_ptr : int*
        label counts for samples right of splitting point.
    n_left : int
        number of samples left of splitting point.
    n_right : int
        number of samples right of splitting point.
    """

    cdef int* pm_left
    cdef int* pm_right

    cdef int n_left, n_right, K, n_samples
    
    def __init__(self, 
                 int K, 
                 np.ndarray[np.int32_t, ndim=1] pm_left,
                 np.ndarray[np.int32_t, ndim=1] pm_right):
        self.K = K
        self.n_left = 0
        self.n_right = 0
        self.pm_left = <int*>pm_left.data
        self.pm_right = <int*>pm_right.data

    cdef void init(self, np.ndarray[np.float_t, ndim=1] labels, 
                   int *sorted_features_i):
        """
        Initializes the criterion for a new feature (col of `features`).
        """
        
        self.n_samples = labels.shape[0]
        self.n_left = 0
        self.n_right = 0
        
        cdef int c = 0      
        for c from 0 <= c < self.K:
            self.pm_left[c] = 0
            self.pm_right[c] = 0
        
        cdef int j = 0
        for j from 0 <= j < self.n_samples:
            c = <int>(labels[j])
            self.pm_right[c] += 1
            self.n_right += 1
        
        """    
        print "ClassificationCriterion.init: "
        print "    pm_left [",    
        for c from 0 <= c < self.K:    
            print self.pm_left[c],  
        print "] ", self.n_left          
        print "    pm_right [",    
        for c from 0 <= c < self.K:    
            print self.pm_right[c],
        print "] ", self.n_right        
        """

    cdef int update(self, int a, int b, 
                     np.float64_t *labels,
                     int *sorted_features_i):

        cdef int c
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            s = sorted_features_i[idx]
            c = <int>(labels[s])
            self.pm_right[c] -= 1
            self.pm_left[c] += 1
            self.n_right -= 1
            self.n_left += 1

        """  
        print "ClassificationCriterion.update: a = ", a, " b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_samples:    
            print sorted_features_i[c],  
        print "] "                  
        print "    labels [",    
        for c from 0 <= c < self.n_samples:    
            print <int>labels[sorted_features_i[c]],
        print "] "               

        print "    pm_left [",    
        for c from 0 <= c < self.K:    
            print self.pm_left[c],  
        print "] ", self.n_left          
        print "    pm_right [",    
        for c from 0 <= c < self.K:    
            print self.pm_right[c],
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
        
        for k from 0 <= k < self.K:
            if self.pm_left[k] > 0:
                H_left -= (self.pm_left[k] / n_left) * (self.pm_left[k] / n_left)
            if self.pm_right[k] > 0:
                H_right -= (self.pm_right[k] / n_right) * (self.pm_right[k] / n_right)

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
        
        for k from 0 <= k < self.K:
            if self.pm_left[k] > 0:
                H_left -= (self.pm_left[k] / n_left) * log(self.pm_left[k] / n_left)
            if self.pm_right[k] > 0:    
                H_right -= (self.pm_right[k] / n_right) * log(self.pm_right[k] / n_right)

        e1 = (n_left / self.n_samples) * H_left
        e2 = (n_right / self.n_samples) * H_right
        return e1 + e2


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression.

    Attributes
    ----------
    sum_left : double
        The sum of the samples left of the SP.
    sum_right : double
        The sum of the samples right of the SP.
    n_left : int
        number of samples left of splitting point.
    n_right : int
        number of samples right of splitting point.
    """

    cdef double sum_left
    cdef double sum_right

    cdef double* labels

    cdef int n_left, n_right, n_samples
    
    def __init__(self, np.ndarray[np.float_t, ndim=1] labels):
        self.sum_left = 0.0
        self.sum_right = 0.0
        self.n_left = 0
        self.n_right = 0
        self.labels = <double*> labels.data
        
    cdef void init(self, np.ndarray[np.float64_t, ndim=1] labels, 
                   int *sorted_features_i):
        """Initializes the criterion for a new feature (col of `features`)."""
        
        self.n_samples = labels.shape[0]
        self.n_left = 0
        self.n_right = 0
        self.sum_left = 0.0
        self.sum_right = 0.0

        cdef int j = 0
        for j from 0 <= j < self.n_samples:
            self.labels[j] = labels[sorted_features_i[j]]
            self.sum_right += self.labels[j]
            self.n_right += 1 

        """
        print "RegressionCriterion.init: "
        print "    sum_left = ", self.sum_left , self.n_left 
        print "    sum_right = ", self.sum_right , self.n_right     
        """
        
    cdef int update(self, int a, int b, 
                     np.float64_t *labels,
                     int *sorted_features_i):
        cdef double val = 0.0
        # post condition: all samples from [0:b) are on the left side
        for idx from a <= idx < b:
            val = self.labels[idx]
            self.sum_right -= val
            self.sum_left += val
            self.n_right -= 1
            self.n_left += 1

        """
        print "RegressionCriterion.update: a = ", a, " b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_samples:    
            print sorted_features_i[c],  
        print "] "                 
        print "    labels [",    
        for c from 0 <= c < self.n_samples:    
            print labels[sorted_features_i[c]],
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

        cdef double var_left = 0.0
        cdef double var_right = 0.0
        cdef int j
        cdef double e1, e2

        for j from 0 <= j < self.n_samples:
            #print "labels[",j,"] = ", self.labels[j]
            if j < self.n_left:
                var_left += (self.labels[j] - mean_left) * (self.labels[j] - mean_left)
            else: 
                var_right += (self.labels[j] - mean_right) * (self.labels[j] - mean_right)

        #print "MSE.eval: var_left = ", var_left, 
        #print "MSE.eval: var_right = ", var_right   

        if self.n_left > 0:
            var_left /= self.n_left
        if self.n_right > 0:
            var_right /= self.n_right

        e1 = ((<double> self.n_left) / self.n_samples) * var_left
        e2 = ((<double> self.n_right) / self.n_samples) * var_right

        return e1 + e2


cdef int smallest_sample_larger_than(int sample_idx,
                                     np.float64_t *features_i,
                                     int *sorted_features_i,
                                     int n_samples):
    """Find index in the `sorted_features` matrix for sample
    who's feature `i` value is just about
    greater than those of the sample `sorted_features_i[sample_idx]`.

    Returns
    -------
    next_sample_idx : int
        The index of the next smallest sample in `sorted_features`
        with different feature value than `sample_idx` .
        I.e. `sorted_features_i[sample_idx] < sorted_features_i[next_sample_idx]`
        -1 if no such element exists.
    """
    cdef int idx = 0
    cdef np.float64_t threshold = -DBL_MAX
    if sample_idx > -1:
        threshold = features_i[sorted_features_i[sample_idx]]
    for idx from sample_idx < idx < n_samples:
        if features_i[sorted_features_i[idx]] > threshold:
            return idx
    return -1


def _find_best_split(np.ndarray[np.float_t, ndim=2, mode="fortran"] features,
                     np.ndarray[np.float_t, ndim=1, mode="c"] labels,
                     Criterion criterion):
    """
    Parameters
    ----------
    features : ndarray, shape (n_samples, n_features), dtype=np.float
        The feature values.
    labels : ndarray, shape (n_samples,), dtype=float
        The labels.
    criterion : Criterion
        The criterion function to be minimized.

    Returns
    -------
    best_i : int
        The split feature or -1 if criterion not smaller than `parent_split_error`.
    best_t : np.float_t
        The split threshold
    best_error : np.float_t
        The error of the split.
    """
    cdef int n_samples = features.shape[0]
    cdef int n_features = features.shape[1]
    cdef int i, a, b, best_i = -1
    cdef np.float_t t, initial_error, error, best_error, best_t

    # pointer access to ndarray data
    cdef double *labels_ptr = <double*>labels.data 
    cdef np.float_t *features_i = NULL
    cdef np.ndarray[np.int32_t, ndim=2, mode="fortran"] sorted_features = \
        np.asfortranarray(np.argsort(features, axis=0).astype(np.int32))
    cdef int *sorted_features_i = NULL    
    
    # Compute the column strides (inc in pointer elem to get from col i to i+1)
    # for `features` and `sorted_features`
    cdef int features_elem_stride = features.strides[0]
    cdef int features_col_stride = features.strides[1]
    cdef int features_stride = features_col_stride / features_elem_stride
    cdef int sorted_features_elem_stride = sorted_features.strides[0]
    cdef int sorted_features_col_stride = sorted_features.strides[1]
    cdef int sorted_features_stride = sorted_features_col_stride / sorted_features_elem_stride   
    
    #compute the initial entropy in the node
    sorted_features_i = (<int *>sorted_features.data)
    criterion.init(labels, sorted_features_i)
    initial_error = criterion.eval()
    if initial_error == 0: # break early if the node is pure
        return best_i, best_t, best_error, initial_error 
    best_error = initial_error
    #print 'at init, best error = ', best_error
    
    for i from 0 <= i < n_features:
        # get i-th col of features and features_sorted
        features_i = (<np.float_t *>features.data) + features_stride * i
        sorted_features_i = (<int *>sorted_features.data) + sorted_features_stride * i

        # init the criterion for this feature
        criterion.init(labels, sorted_features_i)

        # get sample in mask with smallest value for i-th feature
        a = smallest_sample_larger_than(-1, features_i, sorted_features_i,
                                        n_samples)
        while True:
            b = smallest_sample_larger_than(a, features_i, sorted_features_i,
                                            n_samples)
            # if -1 there's none and we are finished
            if b == -1:
                break

            criterion.update(a, b, labels_ptr, sorted_features_i)

            # get criterion value
            error = criterion.eval()
            #print 'error = ', error
            
            # check if current error is smaller than previous best
            # if this is never true best_i is -1.
            if error < best_error:
                t = (features_i[sorted_features_i[a]] +
                     features_i[sorted_features_i[b]]) / 2.0
                best_i = i
                best_t = t
                best_error = error

            a = b
                
    return best_i, best_t, best_error, initial_error
