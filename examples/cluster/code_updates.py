import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.datasets import make_blobs
from hdbscan import HDBSCAN

# Generate sample data
data, _ = make_blobs(n_samples=300, centers=3, random_state=42)

# Fit HDBSCAN
hdbscan = HDBSCAN(min_cluster_size=10)
labels = hdbscan.fit_predict(data)

# Extract hierarchy information
hierarchy = hdbscan.condensed_tree_.to_pandas()

# Create a directed graph
G = nx.DiGraph()

# Add nodes to the graph
for idx, row in hierarchy.iterrows():
    G.add_node(idx, size=row['child_size'])

# Add edges based on parent-child relationships
for idx, row in hierarchy.iterrows():
    parent = row['parent']
    if parent >= 0:
        G.add_edge(parent, idx)

# Create a plot
pos = nx.spring_layout(G)
sizes = np.array([data for _, data in G.nodes(data='size')])
nx.draw(G, pos, node_size=sizes * 100, with_labels=False)
plt.show()
