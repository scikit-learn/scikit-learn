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

# use the Iris dataset as a simple example

iris = datasets.load_iris()

samples = iris.data
labels = iris.target
n_sample = samples.shape[0]
n_feature = samples.shape[1]
n_estimator = 200

plt.figure(figsize=(15, 10))

# train a random trees embedding model of 200 estimators

clf = RandomTreesEmbedding(n_estimators=n_estimator, max_depth = 10)
clf = clf.fit(samples)

# compute the first distance metrics, which is 
# one over the depth of the nearest common ancestor
# averaged over all estimators

forest = clf.estimators_
d = np.zeros([n_sample, n_sample])

for k in range(n_estimator):
    
    tree = forest[k]
    
    path = tree.decision_path(samples).todense()
    n_nodes = tree.tree_.node_count
    
    for i in range(n_sample):
        for j in range(n_sample):
            sample_ids = [i, j]
            d[i, j] = d[i, j]+(path[sample_ids].sum(axis=0) == len(sample_ids)).sum()
            

d = d / n_estimator
d_Nearest_Common_Ancestor = [1/x for x in d]

# plot the heatmap 

plt.subplot(1, 3, 1)
plt.imshow(d_Nearest_Common_Ancestor)

# use agglomerative clustering to classify the output

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="complete")
predict_labels = cluster.fit_predict(d_Nearest_Common_Ancestor)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)

# compute the second distance metrics, which is 
# length of shortest path between leaf nodes
# averaged over all estimators

d = np.zeros([n_sample, n_sample])

for k in range(n_estimator):
    
    tree = forest[k]
    
    path = tree.decision_path(samples).todense()
    n_nodes = tree.tree_.node_count
    l_children = tree.tree_.children_left
    r_children = tree.tree_.children_right
    
    node_depth = np.zeros(n_nodes)
    is_leaves = np.zeros(n_nodes)
    stack = [(0, -1)]
    
    while len(stack) > 0:
        node_id, parent_depth = stack.pop()
        node_depth[node_id] = parent_depth + 1

        if (l_children[node_id] != r_children[node_id]):
            stack.append((l_children[node_id], parent_depth + 1))
            stack.append((r_children[node_id], parent_depth + 1))
        else:
            is_leaves[node_id] = 1
    
    for i in range(n_sample):
        for j in range(n_sample):
            sample_ids = [i, j]
            splitting_depth = (path[sample_ids].sum(axis=0) == len(sample_ids)).sum()
            leaf_i = tree.apply(samples[i].reshape(1, -1))
            leaf_j = tree.apply(samples[j].reshape(1, -1))
            d[i, j] = d[i, j]+node_depth[leaf_i]+node_depth[leaf_j]-2*splitting_depth
            
d_shortest_path = d / n_estimator

# plot the heatmap

plt.subplot(1, 3, 2)
plt.imshow(d_shortest_path)

# use agglomerative clustering to classify the output

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="complete")
predict_labels = cluster.fit_predict(d_shortest_path)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)

# compute the third distance metrics, which is 
# the probability of two samples landing into
# the same leaf node in an estimator

Y = clf.transform(samples)
print(Y.shape)
prob = np.zeros([n_sample, n_sample])
for i in range(n_sample):
    for j in range(n_sample):
        leaf_i = Y[i, :]
        leaf_j = Y[j, :]
        p = leaf_i.dot(leaf_j.transpose()).todense()
        prob[i, j] = p[0, 0]

prob = 1 - prob / Y.shape[1]
        
# plot the heatmap

plt.subplot(1, 3, 3)
plt.imshow(prob)

# use agglomerative clustering to classify the output

cluster = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="average")
predict_labels = cluster.fit_predict(1-prob)
score = adjusted_rand_score(labels, predict_labels)

print("labels\n", labels)
print("predict\n", predict_labels)
print("Adjusted Rand Score:", score)