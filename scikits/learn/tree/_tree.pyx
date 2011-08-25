# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
#
# Adapted for CART by Brian Holt <bdholt1@gmail.com>
#

import numpy as np
cimport numpy as np

cimport cython

cdef extern from "math.h":
    cdef extern double log(double x)
    cdef extern double pow(double base, double exponent)

cdef extern from "float.h":
    cdef extern double DBL_MIN


cdef class Criterion:
    """Interface for splitting criteria, both regression and classification.
    """
    
    cdef void init(self, np.float64_t *labels, 
                   np.int8_t *sample_mask,
                   int *sorted_features_i,
                   int n_total_samples):
        """Init the criteria for each feature `i`. """
        pass

    cdef int update(self, int a, int b, np.float64_t *labels,
                     np.int8_t *sample_mask, np.float64_t *features_i,
                     int *sorted_features_i):
        """Update the criteria for each interval [a,b) where a and b
        are indices in `sorted_features_i`. """
        pass

    cdef double eval(self):
        """Evaluate the criteria (aka the split error). """
        pass

cdef class ClassificationCriterion(Criterion):
    """Abstract criterion for classification.

    Attributes
    ----------
    pm_left_ptr : int*
        label counts for samples left of splitting point.
    pm_right_ptr : int*
        label counts for samples right of splitting point.
    nml : int
        number of samples left of splitting point.
    nmr : int
        number of samples right of splitting point.
    """

    cdef int *pm_left_ptr
    cdef int *pm_right_ptr

    cdef int nml, nmr, K, n_total_samples, n_samples
    
    def __init__(self, int K, np.ndarray[np.int32_t, ndim=1] pm_left,
                 np.ndarray[np.int32_t, ndim=1] pm_right):
        self.K = K
        self.nml = 0
        self.nmr = 0
        self.pm_left_ptr = <int *>pm_left.data
        self.pm_right_ptr = <int *>pm_right.data

    cdef void init(self, np.float64_t *labels, 
                   np.int8_t *sample_mask,
                   int *sorted_features_i,
                   int n_total_samples):
        """Initializes the criterion for a new feature (col of `features`)."""
        cdef int c = 0, j = 0
        self.nml = 0
        self.nmr = 0
        for c from 0 <= c < self.K:
            self.pm_left_ptr[c] = 0
            self.pm_right_ptr[c] = 0

        for j from 0 <= j < n_total_samples:
            if sample_mask[j] == 0:
                continue
            c = <int>(labels[j])
            self.pm_right_ptr[c] += 1
            self.nmr += 1
            
        self.n_samples = self.nmr
        self.n_total_samples = n_total_samples
        
        """
        print "ClassificationCriterion.init: "
        print "    pm_left [",    
        for c from 0 <= c < self.K:    
            print self.pm_left_ptr[c],  
        print "] ", self.nml          
        print "    pm_right [",    
        for c from 0 <= c < self.K:    
            print self.pm_right_ptr[c],
        print "] ", self.nmr        
        """
        
    cdef int update(self, int a, int b, np.float64_t *labels,
                    np.int8_t *sample_mask, np.float64_t *features_i,
                    int *sorted_features_i):

        cdef int c, idx, j
        # move samples from a to b-1 to the left side
        for idx from a <= idx < b:
            j = sorted_features_i[idx]
            if sample_mask[j] == 0:
                continue
            c = <int>(labels[j])
            self.pm_right_ptr[c] -= 1
            self.pm_left_ptr[c] += 1
            self.nmr -= 1
            self.nml += 1
        
        """    
        print "ClassificationCriterion.update: a = ", a, ", b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_total_samples:
            if sample_mask[sorted_features_i[c]] == 0:
                print " *",   
            print sorted_features_i[c],  
        print "] "                  
        print "    labels [",    
        for c from 0 <= c < self.n_total_samples:
            if sample_mask[sorted_features_i[c]] == 0:
                print " *",   
            print <int>labels[sorted_features_i[c]],
        print "] "               

        print "    pm_left [",    
        for c from 0 <= c < self.K:    
            print self.pm_left_ptr[c],  
        print "] ", self.nml         
        print "    pm_right [",    
        for c from 0 <= c < self.K:    
            print self.pm_right_ptr[c],
        print "] ", self.nmr   
        """
            
        return self.nml

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
        cdef double nml = <double> self.nml
        cdef double nmr = <double> self.nmr
        
        for k from 0 <= k < self.K:
            H_left -= (self.pm_left_ptr[k] / nml) * (self.pm_left_ptr[k] / nml)         
            H_right -= (self.pm_right_ptr[k] / nmr) * (self.pm_right_ptr[k] / nmr)    

        e1 = (nml / self.n_samples) * H_left
        e2 = (nmr / self.n_samples) * H_right
                
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
        cdef double nml = <double> self.nml
        cdef double nmr = <double> self.nmr
        
        for k from 0 <= k < self.K:
            if self.pm_left_ptr[k] > 0:
                H_left += -(self.pm_left_ptr[k] / nml) * log(self.pm_left_ptr[k] / nml)
            if self.pm_right_ptr[k] > 0:
                H_right += -(self.pm_right_ptr[k] / nmr) * log(self.pm_right_ptr[k] / nmr)

        e1 = (nml / self.n_samples) * H_left
        e2 = (nmr / self.n_samples) * H_right
        return e1 + e2


cdef class RegressionCriterion(Criterion):
    """Abstract criterion for regression.

    Attributes
    ----------
    sum_left : double
        The sum of the samples left of the SP.
    sum_right : double
        The sum of the samples right of the SP.
    nml : int
        number of samples left of splitting point.
    nmr : int
        number of samples right of splitting point.
    """

    cdef double sum_left
    cdef double sum_right

    cdef np.float64_t *labels
    cdef np.int8_t *sample_mask

    cdef int nml, nmr, K, n_total_samples, n_samples
    
    def __init__(self, np.ndarray[np.float64_t, ndim=1] labels, 
                       np.ndarray[np.int8_t, ndim=1] sample_mask):
        self.sum_left = 0.0
        self.sum_right = 0.0
        self.nml = 0
        self.nmr = 0
        self.labels = <np.float64_t*> labels.data
        self.sample_mask = <np.int8_t*> sample_mask.data
        
    cdef void init(self, np.float64_t *labels, 
                   np.int8_t *sample_mask,
                   int *sorted_features_i,
                   int n_total_samples):
        """Initializes the criterion for a new feature (col of `features`)."""
        cdef int j = 0
        self.nml = 0
        self.nmr = 0
        self.sum_left = 0.0
        self.sum_right = 0.0

        for j from 0 <= j < n_total_samples:
            self.sample_mask[j] = sample_mask[j]
            self.labels[j] = labels[sorted_features_i[j]]
            if sample_mask[j] == 0:
                continue              
            self.sum_right += self.labels[j]
            self.nmr += 1
            
        self.n_samples = self.nmr
        self.n_total_samples = n_total_samples
        
        """
        print "RegressionCriterion.init: "
        print "    sum_left = ", self.sum_left , self.nml
        print "    sum_right = ", self.sum_right , self.nmr     
        """
        
    cdef int update(self, int a, int b, np.float64_t *labels,
                    np.int8_t *sample_mask, np.float64_t *features_i,
                    int *sorted_features_i):

        cdef int idx, j
        cdef double val
        # move samples from a to b-1 to the left side
        for idx from a <= idx < b:
            j = sorted_features_i[idx]
            if sample_mask[j] == 0:
                continue
            val = self.labels[j]
            self.sum_right -= val
            self.sum_left += val
            self.nmr -= 1
            self.nml += 1
        
        """    
        print "RegressionCriterion.update: a = ", a, ", b = ", b
        print "    sorted [",    
        for c from 0 <= c < self.n_total_samples:
            if sample_mask[sorted_features_i[c]] == 0:
                print " *",               
            print sorted_features_i[c],  
        print "] "                 
        print "    labels [",    
        for c from 0 <= c < self.n_total_samples:
            if sample_mask[sorted_features_i[c]] == 0:
                print " *",            
            print labels[sorted_features_i[c]],
        print "] "               

        print "    sum_left = ", self.sum_left , self.nml 
        print "    sum_right = ", self.sum_right , self.nmr                 
        """
            
        return self.nml

    cdef double eval(self):
        pass


cdef class MSE(RegressionCriterion):

    cdef double eval(self):
        """             
        MSE =  \sum_i (y_i - c0)^2  / N
        
        """
        cdef double mean_left = self.sum_left / self.nml
        cdef double mean_right = self.sum_right / self.nmr

        """
        print "MSE.eval: mean_left = ", mean_left
        print "MSE.eval: mean_right = ", mean_right  
        """
        
        cdef double var_left = 0.0
        cdef double var_right = 0.0
        cdef int j, c
        cdef double e1, e2

        c = 0
        for j from 0 <= j < self.n_total_samples:
            if self.sample_mask[j] == 0:
                continue     
            var_left += (self.labels[j] - mean_left) * (self.labels[j] - mean_left)
            c += 1
            if c >= self.nml:
                break
                
        for j from c <= j < self.n_total_samples:
            if self.sample_mask[j] == 0:
                continue                
            var_right += (self.labels[j] - mean_right) * (self.labels[j] - mean_right)

        """
        print "MSE.eval: var_left = ", var_left, 
        print "MSE.eval: var_right = ", var_right  
        """
                
        var_left /= self.nml
        var_right /= self.nmr
        
        """
        print "num left = ", self.nml
        print "num right = ", self.nmr
        print "num samples = ", self.n_samples           
        """
        
        e1 = ((<double> self.nml) / self.n_samples) * var_left
        e2 = ((<double> self.nmr) / self.n_samples) * var_right

        return e1 + e2


cdef int next_sample(int sample_idx,
                     np.float64_t *features_i,
                     int *sorted_features_i,
                     np.int8_t *sample_mask,
                     int n_total_samples):
    """Find index in the `sorted_features` matrix for sample

    Returns
    -------
    next_sample_idx : int
        The index of the next smallest sample in `sorted_features`
        -1 if no such element exists.
    """
    cdef int idx = 0

    for idx from sample_idx < idx < n_total_samples:
        if sample_mask[sorted_features_i[idx]] != 0:
            return idx         
    return -1


def fill_counts(np.ndarray[np.int32_t, ndim=1, mode="c"] counts,
                np.ndarray[np.float64_t, ndim=1, mode="c"] labels,
                np.ndarray sample_mask):
    """The same as np.bincount but casts elements in `labels` to integers.

    Parameters
    ----------
    counts : ndarray, shape = K
        The bin counts to be filled.
    labels : ndarray, shape = n_total_samples
        The labels.
    sample_mask : ndarray, shape=n_total_samples, dtype=bool
        The samples to be considered.
    """
    cdef int j
    cdef int n_total_samples = labels.shape[0]
    cdef char *sample_mask_ptr = <char *>sample_mask.data
    for j from 0 <= j < n_total_samples:
        if sample_mask_ptr[j] == 0:
            continue
        c = <int>labels[j]
        counts[c] += 1


def _find_best_split(np.ndarray sample_mask,
                     np.ndarray[np.float64_t, ndim=2, mode="fortran"] features,
                     np.ndarray[np.int32_t, ndim=2, mode="fortran"] sorted_features,
                     np.ndarray[np.float64_t, ndim=1, mode="c"] labels,
                     Criterion criterion, int K, int n_samples):
    """
    Parameters
    ----------
    sample_mask : ndarray, shape (n_samples,), dtype=bool
        A mask for the samples to be considered. Only samples `j` for which
        sample_mask[j] != 0 are considered.
    parent_split_error : np.float64
        The split error (aka criterion) of the parent node.
    features : ndarray, shape (n_samples, n_features), dtype=np.float64
        The feature values (aka `X`).
    sorted_features : ndarray, shape (n_samples, n_features)
        Argsort of cols of `features`. `sorted_features[0,j]` gives the example
        index of the smallest value of feature j.
    labels : ndarray, shape (n_samples,), dtype=float64
        The labels.
    criterion : Criterion
        The criterion function to be minimized.
    K : int
        The number of classes - for regression use 0.
    n_samples : int
        The number of samples in the current sample_mask (i.e. `sample_mask.sum()`).

    Returns
    -------
    best_i : int
        The split feature or -1 if criterion not smaller than `parent_split_error`.
    best_t : np.float64_t
        The split threshold
    best_error : np.float64_t
        The error (criterion) of the split.
    """
    cdef int n_total_samples = features.shape[0]
    cdef int n_features = features.shape[1]
    cdef int i, j, best_i = -1, best_nml, nml = 0
    cdef np.float64_t t, error, best_error, best_t

    # pointer access to ndarray data
    cdef np.float64_t *labels_ptr = <np.float64_t *>labels.data
    cdef np.float64_t *features_i = NULL
    cdef int *sorted_features_i = NULL
    cdef np.int8_t *sample_mask_ptr = <np.int8_t *>sample_mask.data

    # Compute the column strides (inc in pointer elem to get from col i to i+1)
    # for `features` and `sorted_features`
    cdef int features_elem_stride = features.strides[0]
    cdef int features_col_stride = features.strides[1]
    cdef int features_stride = features_col_stride / features_elem_stride
    cdef int sorted_features_elem_stride = sorted_features.strides[0]
    cdef int sorted_features_col_stride = sorted_features.strides[1]
    cdef int sorted_features_stride = sorted_features_col_stride / sorted_features_elem_stride   
    
    cdef double parent_split_error = np.inf
    
    for i from 0 <= i < n_features:
        # get i-th col of features and features_sorted
        features_i = (<np.float64_t *>features.data) + features_stride * i
        sorted_features_i = (<int *>sorted_features.data) + sorted_features_stride * i

        # init the criterion for this feature
        criterion.init(labels_ptr, 
                       sample_mask_ptr, 
                       sorted_features_i, 
                       n_total_samples)

        # get sample in mask with smallest value for i-th feature
        a = next_sample(-1, features_i, sorted_features_i,
                        sample_mask_ptr, n_total_samples)
        while True:
            # get sample in mask with val for i-th feature just larger than `a`
            b = next_sample(a, features_i, sorted_features_i,
                            sample_mask_ptr, n_total_samples)
            # if -1 there's none and we are fin
            if b == -1:
                break

            # update criterion for interval [a, b)
            #print "updating: a = ", a, " b = ", b
            nml = criterion.update(a, b, labels_ptr, 
                                   sample_mask_ptr, 
                                   features_i,
                                   sorted_features_i)

            # get criterion value
            error = criterion.eval()
            
            # check if current criterion smaller than parent criterion
            # if this is never true best_i is -1.
            if error < parent_split_error:
                # compute split point
                t = (features_i[sorted_features_i[a]] +
                     features_i[sorted_features_i[b]]) / 2.0
                parent_split_error = error
                best_i = i
                best_t = t
                best_error = error
                best_nml = nml
                
            a = b
    return best_i, best_t, best_error, best_nml
