"""
=========================================================================
Clustering on RandomTreesEmbedding Outputs with Modified Distance Metrics
=========================================================================
This demonstration is inspired by the Sparse Projection Oblique
Randomer Forest(SPORF) package. 

Distance is an important concept in clustering and other machine
learning tasks that gives us an intuitive sense of the similarity
between samples. This example shows how can we generate different 
distance metrics from the tree structure info of 
Random Tree Embedding and use these distance metrics to complete 
various tasks. In this example, we will use them to do clustering. 

The distances between two samples i and j can be defined as: 
1) the depth of the nearest common ancestor of the leaf nodes
sample i and j land into (take the reciprocal to make it a 
distance); 
2) the length of the shortest path of the leaf nodes
sample i and j land into; 
3) the probability of sample i and j landing into the same
lead node(take 1 minus the value to make it a distance). 

The example will be implemented on the Iris dataset as a
simple demonstration. And the performance will be evaluated
using ARI. 

"""

from sklearn import tree
from sklearn import datasets
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

# the function here computes the first distance metrics, 
# which is one over the depth of the nearest common ancestor
# averaged over all estimators

# the input are the estimators and the data

def NCA(forest, samples):
    
    n_sample = samples.shape[0]
    d = np.zeros([n_sample, n_sample])
    n_estimator = len(forest)

    for k in range(n_estimator):
    
        tree = forest[k]
    
        path = tree.decision_path(samples).todense()
    
        for i in range(n_sample):
            for j in range(n_sample):
                sample_ids = [i, j]
                d[i, j] = d[i, j]+1/(path[sample_ids].sum(axis=0) == len(sample_ids)).sum()
            

    d = d / n_estimator
    d_Nearest_Common_Ancestor = [1/x for x in d]
    
    return d_Nearest_Common_Ancestor

# the function here computes the second distance metrics, 
# which is length of shortest path between leaf nodes
# averaged over all estimators

# the input are the estimators and the data

def SP(forest, samples): 
    
    n_sample = samples.shape[0]
    d = np.zeros([n_sample, n_sample])
    n_estimator = len(forest)
    
    for k in range(n_estimator):
    
        tree = forest[k]
    
        path = tree.decision_path(samples).todense()
        
        for i in range(n_sample):
            for j in range(n_sample):
                sample_ids = [i, j]
                splitting_depth = (path[sample_ids].sum(axis=0) == len(sample_ids)).sum()
                depth_i = (path[[i, i]].sum(axis=0) == len(sample_ids)).sum()
                depth_j = (path[[j, j]].sum(axis=0) == len(sample_ids)).sum()
                d[i, j] = d[i, j]+depth_i+depth_j-2*splitting_depth

    d_shortest_path = d / n_estimator
    
    return d_shortest_path

# this function compute the third distance metrics, 
# which is the probability of two samples landing into
# the same leaf node in an estimator

# the inputs are the classifier and the data

def PM(clf, samples): 
    
    n_sample = samples.shape[0]
    Y = clf.transform(samples)
    prob = np.zeros([n_sample, n_sample])
    
    for i in range(n_sample):
        for j in range(n_sample):
            leaf_i = Y[i, :]
            leaf_j = Y[j, :]
            p = leaf_i.dot(leaf_j.transpose()).todense()
            prob[i, j] = p[0, 0]

    prob = 1 - prob / Y.shape[1]
    
    return prob

# use the Iris dataset as a simple example

iris = datasets.load_iris()

samples = iris.data
labels = iris.target
n_estimator = 200

plt.figure(figsize=(15, 10))

# train a random trees embedding model of 200 estimators

clf = RandomTreesEmbedding(n_estimators=n_estimator, max_depth = 10)
clf = clf.fit(samples)
forest = clf.estimators_

# compute the distances

d_Nearest_Common_Ancestor = NCA(forest, samples) # the depth of the nearest common ancestor
d_shortest_path = SP(forest, samples) # the shortest path 
prob = PM(clf, samples) # the proximity matrix

# plot the heatmap

plt.subplot(1, 3, 1)
plt.imshow(d_Nearest_Common_Ancestor)
plt.subplot(1, 3, 2)
plt.imshow(d_shortest_path)
plt.subplot(1, 3, 3)
plt.imshow(prob)

# given the nature of the task, here Agglomerative Clustering is selected
# as the method to perform the clustering task. 
# since the affinity metric is the distances we computed beforehand, so 
# usually we choose "precomputed" as the affinity metric. 
# as for the linkage, it should be chosen from "complete", "average" and "single". 
# the parameters used in this example are chosen because they give the best result. 
# in practice, the parameters should be chosen according to the performance. 

# use agglomerative clustering on the nearest common ancestor

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="complete")
predict_labels = cluster.fit_predict(d_Nearest_Common_Ancestor)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)

# use agglomerative clustering on the shortest path

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="complete")
predict_labels = cluster.fit_predict(d_shortest_path)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)

# use agglomerative clustering on the proximity matrix

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="average")
predict_labels = cluster.fit_predict(1-prob)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)