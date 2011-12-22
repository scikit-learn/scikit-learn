"""
================================
Hierarchical Online Clustering
================================

This example shows how one can "augment" an existing clustering with new data.
This can be useful in situations where not all data is available immediately 
but becomes available gradually. One desirable property is that clusterings
in subsequent steps remain consistent, i.e. that two datapoints that have
been in the same cluster originally remain in the same cluster and that two 
datapoints that have been in different clusters remain in different clusters.
This is achieved in this example by defining constraints that are passed to
create_dendrogram.
"""

import numpy as np
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
from sklearn.cluster.hierarchical import create_dendrogram
from sklearn.datasets.samples_generator import make_swiss_roll
from sklearn.neighbors import kneighbors_graph

np.random.seed(0)

###############################################################################
# Generate data (swiss roll dataset)
n_samples = 1000
noise = 0.05
X, _ = make_swiss_roll(n_samples, noise)
# Make it thinner
X[:, 1] *= .5

X1 = X[:500,:]
X2 = X[500:,:]
X = np.vstack((X1,X2))

# Unconstrained clustering of first 500 datapoints
connectivity = kneighbors_graph(X1, n_neighbors=10)

dendrogram = create_dendrogram(X1, connectivity,
                               linkage_criterion="complete")
cluster_roots = dendrogram.cut_height(max_height=15.0)
labeling = dendrogram.get_labeling(cluster_roots)

# Plot result
fig = pl.figure(0)
ax = p3.Axes3D(fig)
ax.view_init(7, -80)
for l in np.unique(labeling):
    ax.plot3D(X1[labeling == l, 0], X1[labeling == l, 1], X1[labeling == l, 2],
              'o', color=pl.cm.jet(float(l) / np.max(labeling + 1)))
pl.title("Datapoints 0...500 (unconstrained)")

# Unconstrained clustering of whole dataset
connectivity = kneighbors_graph(X, n_neighbors=10)

dendrogram = create_dendrogram(X, connectivity,
                               linkage_criterion="complete")
cluster_roots = dendrogram.cut_height(max_height=15.0, prune=True)
labeling1 = dendrogram.get_labeling(cluster_roots)

# Plot result
fig = pl.figure(1)
ax = p3.Axes3D(fig)
ax.view_init(7, -80)
for l in np.unique(labeling1):
    ax.plot3D(X[labeling1 == l, 0], X[labeling1 == l, 1], X[labeling1 == l, 2],
              'o', color=pl.cm.jet(float(l) / np.max(labeling1 + 1)))
pl.title('Datapoints 0...1000 (unconstrained)')

# Create constraints: 
# Two datapoints that have been in different clusters must not be merged.
def constraintGenerator(s1, s2):
    def constraint(t1, t2):
        t = set(t1 + t2)
        return t.isdisjoint(s1) or t.isdisjoint(s2) \
                    or (t.issuperset(s1) and t.issuperset(s2)) 
    return constraint

constraints = []
for l1 in np.unique(labeling):
    for l2 in np.unique(labeling):
        if l1 >= l2: continue
        constraints.append(constraintGenerator(set(np.where(labeling==l1)[0]),
                                               set(np.where(labeling==l2)[0])))

        
# Constrained clustering of whole dataset
connectivity = kneighbors_graph(X, n_neighbors=10)
dendrogram = create_dendrogram(X, connectivity,
                               linkage_criterion="complete",
                               constraints=constraints)
cluster_roots = dendrogram.cut_height(max_height=15.0)
labeling2 = dendrogram.get_labeling(cluster_roots)

# Check that the constrained clustering is consistent with the original 
# clustering, i. e. that two datapoints that have been in the same cluster 
# originally are still in the same and that two datapoints that have been in 
# different clusters are still in different clusters.
for i in range(len(labeling)):
    for j in range(i, len(labeling)):
        if (labeling[i] != labeling[j]):
            assert labeling2[i] != labeling2[j]
    
print "Constrained clustering consistent with original clustering."

# Plot result
fig = pl.figure(2)
ax = p3.Axes3D(fig)
ax.view_init(7, -80)
for l in np.unique(labeling2):
    ax.plot3D(X[labeling2 == l, 0], X[labeling2 == l, 1], X[labeling2 == l, 2],
              'o', color=pl.cm.jet(float(l) / np.max(labeling2 + 1)))
pl.title("Datapoints 0...1000 constrained to clustering of datapoints 0...500")

pl.show()
