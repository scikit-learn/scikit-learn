"""
A framework for comparison between different clustering methods. 
"""

print __doc__

import numpy as np
from scikits.learn.cluster import * 
import pylab as pl
from itertools import cycle
import time


def main():
    "main function"

################################################################################
    

    ################################################################################
    # Generate sample data
    np.random.seed(0)

    n_points_per_cluster = 250
    n_clusters = 3
    n_points = n_points_per_cluster*n_clusters
    means = np.array([[1,1],[-1,-1],[1,-1]])
    std = .6
    clustMed = []
    
    X = np.empty((0, 2))
    labels = np.empty(0) # actual labels (ground truth)
    
    
    
    for i in range(n_clusters):
        X = np.r_[X, means[i] + std * np.random.randn(n_points_per_cluster, 2)]
        labels = np.r_[labels, i * np.ones(n_points_per_cluster)] 
        

        
    clust_labels = np.zeros((n_points,0)); # clust  labels[:,i] = labels given by the ith clustering method
    times = [];
    names = [];
        
    ################################################################################
    # Compute clustering with MeanShift

    names.append('Mean Shift');
    start = time.clock()

    bandwidth = estimate_bandwidth(X, quantile=0.3)
    ms = MeanShift(bandwidth=bandwidth)
    ms.fit(X)

    elapsed = (time.clock() - start)





    clust_labels = np.c_[clust_labels, ms.labels_]

    times.append(elapsed)


  

    ################################################################################
    # Compute clustering with Kmeans
    names.append('k-means');

    start = time.clock()
    km = KMeans(init='k-means++', k=3, n_init=1)
    km.fit(X);
    elapsed = (time.clock() - start)
    
                  
    clust_labels = np.c_[clust_labels, km.labels_]
    
    times.append(elapsed)
    
    
    
    
    ################################################################################
    # An example of how to add a new clustering methods (here named X-Cluster)
    ################################################################################
    # Compute clustering with (METHOD X-CLUSTER)
    
    # names.append('X-CLUSTER');

    # start = time.clock()
    # xc = XCluster(parameters)
    # xc.fit(X);
    # elapsed = (time.clock() - start)
    
                  
    # clust_labels = np.c_[clust_labels, xc.labels_]
    # times.append(elapsed)


    
    ################################################################################

    n_clust_methods = clust_labels.shape[1]
    errors = []
    ################################################################################
    for i in range(n_clust_methods):
        #plot_clustering(clust_labels[:,i],X);
        errors.append(calc_error(labels,clust_labels[:,i]));
        
        
    
    
    
    for t,err, name in zip(times, errors, names):
        print "%s: Error=%f, Time=%f"%(name, err, t)



################################################################################
# Plot function


def calc_error(labels1, labels2):
    labels1 = np.copy(labels1)
    labels2 = np.copy(labels2)

    
    labels1,labels2 = unify(labels1,labels2);


    
    #print labels1, labels2

    

    return np.sum(np.abs(labels1-labels2))/labels1.size


def unify(labels1,labels2):
    '''

    '''
    labels1 = np.copy(labels1)
    labels2 = np.copy(labels2)


    labels_unique = np.unique(labels1);
    n_clust = labels_unique.size;


    M = np.zeros((n_clust,n_clust));
    
    for i,j in zip(labels1,labels2):
        M[i,j] += 1;

    #print '****************'
    #print labels1,labels2

    
    for k in range(n_clust):
        amax = M[k,:].argmax();
        labels1[labels1 == k] = -amax
        M[k+1:,amax] = 0
    

    labels1 = -labels1; 


    return labels1,labels2

    

    
    

        




    


def plot_clustering(labels,X):
    pl.figure()
    
    pl.clf()

    n_clusters = X.shape[0]

    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters), colors):
        my_members = labels == k
        pl.plot(X[my_members,0], X[my_members,1], col+'.')
    
    
    pl.title('Estimated number of clusters: %d' % n_clusters)
    pl.show()






main()







