# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 12:58:49 2016

@author: anki08
"""

import numpy as np
from sklearn.utils.testing import assert_equal, assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.metrics.cluster.bicluster import _jaccard
from sklearn.metrics import consensus_score

from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics.cluster import entropy
from sklearn.metrics.cluster import expected_mutual_information
#from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_completeness_v_measure
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import v_measure_score

import scipy.sparse as sp
from scipy.sparse import csr_matrix
from sklearn import datasets
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_raise_message
from sklearn.metrics.cluster import silhouette_score
from sklearn.metrics.cluster import silhouette_samples
from sklearn.metrics import pairwise_distances
#from sklearn.metrics.cluster import calinski_harabaz_score

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
        "jaccard":_jaccard ,
        "consensus_score":consensus_score
        }

      
SUPERVISED_METRICS = {
        "adjusted_mutual_info_score":adjusted_mutual_info_score ,
        "adjusted_rand_score":adjusted_rand_score,
        "completeness_score":completeness_score,
        "contingency_matrix":contingency_matrix,
        "homogeneity_score":homogeneity_score,
        "entropy":entropy,
        "expected_mutual_information":expected_mutual_information,
        "fowlkes_mallows_score":fowlkes_mallows_score,
        "homogeneity_completeness_v_measure":homogeneity_completeness_v_measure,
        "mutual_info_score":mutual_info_score,
        "normalized_mutual_info_score":normalized_mutual_info_score,
        "v_measure_score":v_measure_score
        }
        
UNSUPERVISED_METRICS={
        "silhouette_score":silhouette_score,
        "silhouette_samples":silhouette_samples,
        "calinski_harabaz_score":calinski_harabaz_score
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

# Symmetric with respect to their input arguments y_true and y_pred
# metric(y_true, y_pred) == metric(y_pred, y_true).
SYMMETRIC_METRICS = [ "adjusted_rand_score" , "v-measure_score" ,
                      "mutual_info_score","adjusted_mutual_info_score","normalized_mutual_score",
                    ]
                    
NON_SYMMETRIC_METRICS = ["homogeneity_score","completeness_score"]

#METIRICS with output between 0 and 1 and belong to BICLUSTER METRICS
METRICS_NORMALIZED_OUTPUT_BICLUSTER = [ "consensus_score","jaccard"
                                      ]
 
#metrics with output between 0 and 1                                      
METRICS_NORMALIZED_OUTPUT = ["adjusted_rand_score","homogeneity_score",
                            "completeness_score","v-measure_score",
                            "adjusted_mutual_info_score","normalized_mutual_info_score",
                            "fowlkes_mallows_score"
                            ]
                      
                                
#when information is zero these metrics output zero
METRICS_ZERO_INFO = ["normalized_mutual_info_score","v-measure_score",
                     "adjusted_mutual_info_score"]
                      
#METRICS where permutations oflabels dont change score
METRICS_PERMUTE_LABELS = ["homogeneity_score","v-measure_score",
                         "completeness_score","mutual_info_score","adjusted_mutual_info_score",
                         "normalized_mutual_info_score"
                         ]
                        
#metrics which result in 0 when a class is split across different clusters
ClASS_BASED_METRICS = ["adjusted_mutual_info","normalized_mutual_info",
                      "fowlkes_mallows_score","v-measure_score"
                      ]
                        

#@ignore_warnings
def test_symmetry():
    random_state = check_random_state(0)
    y_true=[0, 1, 0]
    y_pred=[42,7,42]
    
    #symmetric_metrics
    for name in SYMMETRIC_METRICS:
        metric=ALL_METRICS[name]
        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_pred, y_true),
                            err_msg="%s is not symmetric" % name)
                            
    # Not symmetric metrics
    for name in NON_SYMMETRIC_METRICS:
        metric = ALL_METRICS[name]
        assert_true(np.any(metric(y_true, y_pred) != metric(y_pred, y_true)),
                    msg="%s seems to be symmetric" % name)   
                    
#test function for metrics whose output in range 0 to 1
#bicluster metrics have different input format so they are handled differently                  
def test_normalized_output_bicluster():
    
    # test for jaccard
    a1 = np.array([True, True, False, False])
    a2 = np.array([True, True, True, True])
    a3 = np.array([False, True, True, False])
    a4 = np.array([False, False, True, True])
    assert_equal(_jaccard(a1, a1, a1, a1), 1)
    assert_greater(_jaccard(a1, a1, a2, a2), 0)
    assert_less(_jaccard(a1, a1, a2, a2), 1)
    assert_greater(_jaccard(a1, a1, a3, a3),0)
    assert_less(_jaccard(a1, a1, a3, a3), 1)
    assert_equal(_jaccard(a1, a1, a4, a4), 0)
    
    #test for consensus score 
    a = [[True, True, False, False],
         [False, False, True, True]]
    b = a[::-1]
    assert_equal(consensus_score((a, a), (a, a)), 1)
    assert_equal(consensus_score((a, a), (b, b)), 1)
    assert_equal(consensus_score((a, b), (a, b)), 1)
    assert_equal(consensus_score((a, b), (b, a)), 1)

    assert_equal(consensus_score((a, a), (b, a)), 0)
    assert_equal(consensus_score((a, a), (a, b)), 0)
    assert_equal(consensus_score((b, b), (a, b)), 0)
    assert_equal(consensus_score((b, b), (b, a)), 0)
    
    b=a[:,2]
    assert_greater(consensus_score((a , a), (b , b)), 0)
    assert_less(consensus_score((a , a), (b , b)), 1)
    
    
#test function for metrics whose output in range 0 to 1
def test_normalized_output():
    for name in METRICS_NORMALIZED_OUTPUT:
        metric=ALL_METRICS[name]
        assert_greater(metric([0, 0, 0, 1, 1, 1],[0, 0, 0, 1, 2, 2]), 0.0)
        assert_less(metric([0, 0, 0, 1, 1, 1],[0, 1, 0, 1, 2, 2]), 1.0)
        assert_less(metric([0, 0, 1, 1, 2, 2],[0, 0, 1, 1, 1, 1]), 1.0)
        assert_greater(metric([0, 0, 1, 1, 2, 2],[0, 0, 1, 1, 1, 1]), 0.0)
        assert_less(metric([0, 0, 0, 2, 2, 2],[0, 1, 0, 1, 2, 2]),1.0)
        assert_greater(metric([0, 0, 0, 2, 2, 2],[0, 1, 0, 1, 2, 2]), 0.0)

                            
#when information is zero these metrics output zero                     
def test_exactly_zero_info_score():
    # Check numerical stability when information is exactly zero
    for i in np.logspace(1, 4, 4).astype(np.int):
        labels_a, labels_b = (np.ones(i, dtype=np.int),
                              np.arange(i, dtype=np.int))
        assert_equal(normalized_mutual_info_score(labels_a, labels_b), 0.0)
        assert_equal(v_measure_score(labels_a, labels_b), 0.0)
        assert_equal(adjusted_mutual_info_score(labels_a, labels_b), 0.0)
        assert_equal(normalized_mutual_info_score(labels_a, labels_b), 0.0)
    
#test function for mtericshaving the property of not changing the score when 
#the labels are permuted.
                       
def permute_labels():
    for name in METRICS_NORMALIZED_OUTPUT:
        metric=ALL_METRICS[name]
        var_1 = metric([0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 2, 2])
        var_2 = metric([0, 1, 0, 1, 0, 1], [2, 0, 1, 1, 0, 2])
        assert_equal(var_1,var_2)
        
# If classes members are completely split across different clusters,
#the assignment is totally in-complete, hence the score of these metrics is 0
#they are perfect when the clusters are both homoneneous and complete
        
def class_based_clusters():
    for name in METRICS_NORMALIZED_OUTPUT:
        metric=ALL_METRICS[name]
        assert_equal(metric([0, 0, 0, 0], [0, 1, 2, 3]), 0.0)
        assert_equal(metric([0, 0, 1, 1], [0, 0, 1, 1]),1.0)
        assert_equal(metric([0, 0, 1, 1], [0, 0, 1, 1]),1.0)

    