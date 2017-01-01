# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 12:58:49 2016

@author: anki08
"""

import numpy as np

from sklearn.metrics.cluster.bicluster import _jaccard
from sklearn.metrics.cluster import consensus_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics.cluster import entropy
from sklearn.metrics.cluster import expected_mutual_information
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import v_measure_score
from sklearn.metrics.cluster import silhouette_score
from sklearn.metrics.cluster import silhouette_samples

from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_almost_equal


# Dictionaries of metrics
# ------------------------
# The goal of having those dictionaries is to have an easy way to call a
# particular metric and associate a name to each function:
#   - BICLUSTER_METRICS: all biclusters - (clusters formed from matrices)
#   - SUPERVISED_METRICS: all supervised cluster metrics - (when given a ground truth value)
#   - UNSUPERVISED_METRICS: all unsupervised cluster metrics
#
# Those dictionaries will be used to test systematically some invariance
# properties, e.g. invariance toward several input layout.
#

#metrics used to test similarity between bicluster
BICLUSTER_METRICS = {
    "consensus_score":consensus_score    
    }
      
SUPERVISED_METRICS = {
    "adjusted_mutual_info_score":adjusted_mutual_info_score ,
    "adjusted_rand_score":adjusted_rand_score,
    "completeness_score":completeness_score,
    "contingency_matrix":contingency_matrix,
    "homogeneity_score":homogeneity_score,
    "expected_mutual_information":expected_mutual_information,
    "homogeneity_completeness_v_measure":homogeneity_completeness_v_measure,
    "mutual_info_score":mutual_info_score,
    "normalized_mutual_info_score":normalized_mutual_info_score,
    "v_measure_score":v_measure_score
    }
        
UNSUPERVISED_METRICS = {
    "silhouette_score":silhouette_score
    }
        
ALL_METRICS = dict()
ALL_METRICS.update(BICLUSTER_METRICS)
ALL_METRICS.update(SUPERVISED_METRICS)
ALL_METRICS.update(UNSUPERVISED_METRICS)

# Lists of metrics with common properties
# ---------------------------------------
# Lists of metrics with common properties are used to test systematically some
# functionalities and invariance, e.g. SYMMETRIC_METRICS lists all metrics that
# are symmetric with respect to their input argument y_true and y_pred.
#
#---------------------------------------------------------------------
# Symmetric with respect to their input arguments y_true and y_pred.
# Symmetric metrics only imply to supervised clusters.
SYMMETRIC_METRICS = [
    "adjusted_rand_score","v_measure_score",
    "mutual_info_score","adjusted_mutual_info_score",
    "normalized_mutual_info_score",
    ]
                    
NON_SYMMETRIC_METRICS = ["homogeneity_score" , "completeness_score"]

# When the information is zero these metrics output zero.
METRICS_ZERO_INFO = [
    "normalized_mutual_info_score","v_measure_score",
    "adjusted_mutual_info_score"
    ]
 
# Metrics with output between 0 and 1                                      
METRICS_NORMALIZED_OUTPUT = [
    "adjusted_rand_score","homogeneity_score","completeness_score",
    "v_measure_score","adjusted_mutual_info_score",
    "normalized_mutual_info_score"
    ]                  
                      
# Metrics where permutations of labels dont change score( 0 and 1 exchchanged)
METRICS_PERMUTE_LABELS = ["homogeneity_score","v_measure_score",
    "completeness_score","mutual_info_score","adjusted_mutual_info_score",
    "normalized_mutual_info_score","adjusted_rand_score"
    ]

# Input parameters can be both in the form of arrays and lists
METRICS_WITH_FORMAT_INVARIANCE = ["homogeneity_score","v_measure_score",
    "completeness_score","mutual_info_score","adjusted_mutual_info_score",
    "normalized_mutual_info_score","adjusted_rand_score"
    ] 
    
# If classes members are completely split across different clusters,
# the assignment is totally in-complete, hence the score of these metrics is 0
# they are perfect when the clusters are both homoneneous and complete.
#
CLASS_BASED_METRICS = [
    "adjusted_mutual_info_score","normalized_mutual_info_score",
    "v_measure_score"
    ]

                                                
def assert_between(var,score_1,score_2):
    """ Returns a boolean value
    Helper function to check if score lies in between two values
    """
    if assert_greater(var,score_1) and assert_less(var,score_2):
        return True
    else:
        return False

        
def test_symmetry():
    rng = np.random.RandomState(0)
    y1 = rng.randint(3, size=30)
    y2 = rng.randint(3, size=30)
    for name in SYMMETRIC_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric(y1,y2),
                            metric(y2,y1))
            
    for name in NON_SYMMETRIC_METRICS:
        metric = ALL_METRICS[name]
        assert_almost_equal(metric([0, 1, 2, 5, 4, 9],[0, 1, 9, 4, 3, 5]),
                            metric([0, 1, 9, 4, 3, 5],[0, 1, 2, 5, 4, 9]))

                                                
def test_exactly_zero_info_score():
    # Check numerical stability when information is exactly zero
    for i in np.logspace(1, 4, 4).astype(np.int):
        labels_a, labels_b=(np.ones(i, dtype=np.int),
                              np.arange(i, dtype=np.int))
        for name in METRICS_ZERO_INFO:
            metric=ALL_METRICS[name]
            assert_almost_equal(metric(labels_a,labels_b), 0.0)

                                    
def test_normalized_output():
    upper_bound_1 = [0, 0, 0, 1, 1, 1]
    upper_bound_2 = [0, 0, 0, 1, 1, 1]
    for name in METRICS_NORMALIZED_OUTPUT:
        metric=ALL_METRICS[name]
        assert_between(metric([0, 0, 0, 1, 1, 1],[0, 0, 0, 1, 2, 2]), 0.0, 1.0)
        assert_between(metric([0, 0, 1, 1, 2, 2],[0, 0, 1, 1, 1, 1]), 0.0, 1.0)
        assert_equal(metric(upper_bound_1,upper_bound_2), 1.0)

                         
def test_permute_labels():
    for name in METRICS_PERMUTE_LABELS:
        metric = ALL_METRICS[name]
        y_label = np.array([0, 0, 0, 1, 1, 0, 1])
        y_pred=np.array([1, 0, 1, 0, 1, 1, 0])
        assert_almost_equal(metric(y_pred, y_label),metric(y_pred, 1-y_label))
        assert_almost_equal(metric(y_pred, y_label),metric(1-y_pred, y_label))
        assert_almost_equal(metric(y_pred, y_label),metric(1-y_pred, 1-y_label))

        
def test_class_based_clusters():
    for name in CLASS_BASED_METRICS:
        metric = ALL_METRICS[name]
        assert_equal(metric([0, 0, 0, 0],[0, 1, 2, 3]),0.0)
        assert_equal(metric([0, 0, 1, 1],[0, 0, 1, 1]),1.0)
        assert_equal(metric([0, 0, 1, 1],[0, 0, 1, 1]),1.0)   
        

def test_format_invariance():
    for name in METRICS_WITH_FORMAT_INVARIANCE:
        metric = ALL_METRICS[name]
        list_a = [0, 0, 0, 1, 1, 1]
        list_b = [0, 1, 2, 3, 4, 5]
        arr_a = np.array([0, 0, 0, 1, 1, 1])
        arr_b = np.array([0, 1, 2, 3, 4, 5])
        score_list = metric(list_a,list_b)
        score_array = metric(arr_a,arr_b)
        assert_equal(score_list,score_array)
          