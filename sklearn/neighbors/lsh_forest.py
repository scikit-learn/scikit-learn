"""
Locality Sensitive Hashing Forest for Approximate Nearest Neighbor Search
-------------------------------------------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.utils.validation import safe_asarray

__all__ = ["LSHForest"]


def _bisect_left(a, x):
    """Private function to perform bisect left operation"""
    return np.searchsorted(a, x)

def _bisect_right(a, x):
    """Private function to perform bisect right operation""" 
    lo = 0
    hi = a.shape[0]
    while lo < hi:
        mid = (lo+hi)//2
        #Only string length of x is considered when comparing.
        if x < a[mid] and not a[mid][:len(x)]==x:
            hi = mid
        else:
            lo = mid + 1
    return lo

def _find_matching_indices(sorted_array, item, h):
    """
    Finds indices in sorted array of strings where their first
    h elements match the items' first h elements
    """
    left_index = bisect_left(sorted_array, item[:h])
    right_index = bisect_right(sorted_array, item[:h])
    return np.arange(left_index, right_index)

def _find_longest_prefix_match(bit_string_array, query):
    """
    Private function to find the longest prefix match for query
    in the bit_string_array
    """    
    hi = len(query)
    lo = 0
    
    if _find_matching_indices(bit_string_array, query, hi).shape[0] > 0:
        return hi
    
    while lo < hi:
        mid = (lo+hi)//2        
        k = _find_matching_indices(bit_string_array, query, mid).shape[0]
        if k > 0:
            lo = mid + 1
            res = mid
        else:
            hi = mid     
    return res

def _simple_euclidean_distance(query, candidates):
    """
    Private function to calculate Euclidean distances between each
    point in candidates and query
    """
    distances = np.zeros(candidates.shape[0])    
    for i in range(candidates.shape[0]):        
        distances[i] = np.linalg.norm(candidates[i]-query)        
    return distances

class LSH_forest(BaseEstimator):
    """
    Performs approximate nearest neighbor search using LSH forest.
    
    LSH Forest: Locality Sensitive Hashing forest [1] is an alternative 
    method for vanilla approximate nearest neighbor search methods.
    LSH forest data structure has been implemented using sorted 
    arrays and binary search. 
    
    Parameters
    ----------
    max_label_length: int, optional (default = 32)
        Maximum length of a binary hash
    
    n_trees: int, optional (default = 10)
        Number of trees in the LSH Forest.
        
    hashing_algorithm: {'random_projections'}, 
        optional (default = 'random_projections')
        Algorithm of LSH family by which the hashing is performed on
        the data vectors.
        
        -'random_projections': hash using :class:`RandomProjections`
        
    Attributes
    ----------
    max_label_lengh: maximum length of hash 
    number_of_trees: number of trees build in indexing
    """
    def __init__(self, max_label_length = 32, number_of_trees = 10, 
                 hashing_algorithm = 'random_projections'):
        self.max_label_length = max_label_length
        self.number_of_trees = number_of_trees
        self.hashing_algorithm = hashing_algorithm
        self.random_state = np.random.RandomState(seed=1)
        
    def _select_hashing_algorithm():
        """ Selectes the LSH algorithm """
        if self.hashing_algorithm=='random_projections':
            return RandomProjections
        else:
            raise ValueError("Unknown hashing algorithm: %s" 
                             %(self.hashing_algorithm))
        
    def _get_random_hyperplanes(self, hash_size = None, dim = None):
        """ 
        Generates hyperplanes from standard normal distribution  and return 
        it as a 2D numpy array. This is g(p,x) for a particular tree.
        """
        if hash_size == None or dim == None:
            raise ValueError("hash_size or dim(number of dimensions) cannot be None.")        
        
        return self.random_state.randn(hash_size, dim) 
    
    def _hash(self, input_point = None, hash_function = None):
        """
        Does hash on the data point with the provided hash_function: g(p,x).
        """
        if input_point == None or hash_function == None:
            raise ValueError("input_point or hash_function cannot be None.")
            
        projections = np.dot(hash_function, input_point) 
            
        return "".join(['1' if i > 0 else '0' for i in projections])
    
    def _create_tree(self, hash_function = None):
        """
        Builds a single tree (in this case creates a sorted array of 
        binary hashes).
        """
        if hash_function == None:
            raise ValueError("hash_funciton cannot be None.")
            
        number_of_points = self.input_array.shape[0]
        binary_hashes = []
        for i in range(number_of_points):
            binary_hashes.append(self._hash(self.input_array[i], hash_function))
        
        binary_hashes = np.array(binary_hashes)
        return np.argsort(binary_hashes), np.sort(binary_hashes)
    
    def _compute_distances(self, query, candidates):
        #distances = euclidean_distances(query, self.input_array[candidates])
        distances = simple_euclidean_distance(query, self.input_array[candidates])
        return np.argsort(distances), distances
        
        
    def fit(self, X = None):
        """
        Fit the LSH forest on the data. 
        
        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.       
       """ 
        if input_array == None:
            raise ValueError("input_array cannot be None")
            
        self._input_array = safe_asarray(X)
        number_of_points = input_array.shape[0]
        n_dim = input_array.shape[1]
        
        self.hash_generator = self._select_hashing_algorithm(n_dim=n_dim, 
                                                             hash_size=self.max_label_length)
        
        #Creates a g(p,x) for each tree
        self.hash_functions = []
        self.trees = []
        self.original_indices = []
        for i in range(self.number_of_trees):
            #This is g(p,x) for a particular tree.
            hash_function = self.hash_generator.generate_hash_function()
            original_index, bin_hashes = self._create_tree(hash_function)
            self.original_indices.append(original_index)
            self.trees.append(bin_hashes)
            self.hash_functions.append(hash_function)
        
        self.hash_functions = np.array(self.hash_functions)
        self.trees = np.array(self.trees)
        self.original_indices = np.array(self.original_indices)
        
    def _query(self, query = None, c = 1, m = 10, lower_bound = 4):
        """
        returns the number of neighbors for a given query.
        """
        if query == None:
            raise ValueError("query cannot be None.")
        query = np.array(query)
        
        #descend phase
        max_depth = 0
        for i in range(len(self.trees)):
            bin_query = self._hash(query, self.hash_functions[i])
            k = get_longest_prefix_length(self.trees[i], bin_query)
            if k > max_depth:
                max_depth = k
        
        bin_queries = []
        for i in range(len(self.trees)):
            bin_queries.append(self._hash(query, self.hash_functions[i]))

        #Synchronous ascend phase
        candidates = []
        number_of_candidates = c*len(self.trees)
        while max_depth > lower_bound and (len(candidates) < number_of_candidates or len(set(candidates)) < m):
            for i in range(len(self.trees)):                
                candidates.extend(self.original_indices[i,simpleFunctionBisectReImplemented(self.trees[i], 
                                                                                            bin_queries[i], max_depth)].tolist())
                #candidates = list(OrderedSet(candidates)) #this keeps the order inserted into the list 
            max_depth = max_depth - 1
            #print max_depth, len(candidates) ,len(set(candidates))
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)
        #print ranks[0,:m]
        print candidates.shape
        return candidates[ranks[:m]]
    
    
    def query_num_candidates(self, query = None, c = 1, m = 10, lower_bound = 4):
        """
        returns the nearest neighbors for a given query the number of required 
        candidates.
        """
        if query == None:
            raise ValueError("query cannot be None.")
        query = np.array(query)
        
        #descend phase
        max_depth = 0
        for i in range(len(self.trees)):
            bin_query = self._hash(query, self.hash_functions[i])
            k = get_longest_prefix_length(self.trees[i], bin_query)
            if k > max_depth:
                max_depth = k
                
        bin_queries = []
        for i in range(len(self.trees)):
            bin_queries.append(self._hash(query, self.hash_functions[i]))
                
        #Synchronous ascend phase
        candidates = []
        number_of_candidates = c*len(self.trees)
        while max_depth > lower_bound and (len(candidates) < number_of_candidates or len(set(candidates)) < m):
            for i in range(len(self.trees)):
                candidates.extend(self.original_indices[i,simpleFunctionBisectReImplemented(self.trees[i], 
                                                                                            bin_queries[i], max_depth)].tolist())
                #candidates = list(OrderedSet(candidates)) #this keeps the order inserted into the list 
            max_depth = max_depth - 1
            #print max_depth, len(candidates) ,len(set(candidates))
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)
        #print ranks[0,:m]        
        return candidates[ranks[:m]], candidates.shape[0]


    def query_candidates(self, query = None, c = 1, m = 10, lower_bound = 5):
        """
        returns the nearest neighbors for a given query the number of required 
        candidates.
        """
        if query == None:
            raise ValueError("query cannot be None.")
        query = np.array(query)
        
        #descend phase
        max_depth = 0
        for i in range(len(self.trees)):
            bin_query = self._hash(query, self.hash_functions[i])
            k = get_longest_prefix_length(self.trees[i], bin_query)
            if k > max_depth:
                max_depth = k

        bin_queries = []
        for i in range(len(self.trees)):
            bin_queries.append(self._hash(query, self.hash_functions[i]))  

        #Synchronous ascend phase
        candidates = []
        number_of_candidates = c*len(self.trees)
        while max_depth > lower_bound and (len(candidates) < number_of_candidates or len(set(candidates)) < m):
            for i in range(len(self.trees)):            
                candidates.extend(self.original_indices[i,simpleFunctionBisectReImplemented(self.trees[i], 
                                                                                            bin_queries[i], max_depth)].tolist())
                #candidates = list(OrderedSet(candidates)) #this keeps the order inserted into the list 
            max_depth = max_depth - 1
            #print max_depth, len(candidates) ,len(set(candidates))
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)
        #print ranks[0,:m]        
        return candidates[ranks[:m]], candidates


    def get_candidates_for_hash_length(self, query, hash_length):
        candidates = []        
        for i in range(len(self.trees)):
            bin_query = self._hash(query, self.hash_functions[i])
            candidates.extend(self.original_indices[i,simpleFunctionBisectReImplemented(self.trees[i], 
                                                                                            bin_query, hash_length)].tolist())
            
        return np.unique(candidates)

