"""Hierarchical Agglomerative Clustering

These routines perform some hierachical agglomerative clustering of some
input data. Currently, only Ward's algorithm is implemented.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux
License: BSD 3 clause
"""
import heapq
from collections import defaultdict
import itertools
import warnings

import numpy as np
from scipy import sparse
from scipy.cluster import hierarchy

from ..base import BaseEstimator
from ..utils._csgraph import cs_graph_components
from ..externals.joblib import Memory

from . import _inertia
from ._feature_agglomeration import AgglomerationTransform


###############################################################################

class WardDistance(object):
    
    def __init__(self, X, connectivity, n_nodes, n_features, n_samples):
        # build moments as a list
        self.moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features))]
        self.moments[0][:n_samples] = 1
        self.moments[1][:n_samples] = X
    
        self.A = []
        # create distances matrix
        coord_row = []
        coord_col = []
        for ind, row in enumerate(connectivity.rows):
            self.A.append(row)
            # We keep only the upper triangular for the moments
            # Generator expressions are faster than arrays on the following
            row = [i for i in row if i < ind]
            coord_row.extend(len(row) * [ind, ])
            coord_col.extend(row)
    
        coord_row = np.array(coord_row, dtype=np.int)
        coord_col = np.array(coord_col, dtype=np.int)
        
        distances = np.empty(len(coord_row), dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   coord_row, coord_col, distances)
    
        self.distances = zip(distances, coord_row, coord_col)
        heapq.heapify(self.distances)

    def update(self, child_node1, child_node2, parent_node, parent):
        # update the moments
        for p in range(2):
            self.moments[p][parent_node] = \
                        self.moments[p][child_node1] + self.moments[p][child_node2]
        # update the structure matrix A and the inertia matrix
        coord_col = []
        visited = set([parent_node])      
        for l in set(self.A[child_node1]).union(self.A[child_node2]):
            while parent[l] != l:
               l = parent[l]
            if l not in visited:
               visited.add(l)
               coord_col.append(l)
               self.A[l].append(parent_node)
        self.A.append(coord_col)
        coord_col = np.array(coord_col, dtype=np.int)
        coord_row = np.empty_like(coord_col)
        coord_row.fill(parent_node)
        
        distances = np.empty(len(coord_row), dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   coord_row, coord_col, distances)
        
        for tupl in itertools.izip(distances, coord_row, coord_col):
            heapq.heappush(self.distances, tupl)
            
    def hasMoreCandidates(self):
        return len(self.distances) > 0
    
    def fetchCandidate(self):
        return heapq.heappop(self.distances)
    
    def computeDistance(self, i, j):
        # Compute inertia of the cluster obtained when merging i and j
        distances = np.empty(1, dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   np.array([i]), np.array([j]), distances)
        return distances[0]


class CompleteLinkageDistance(object):
    
    def __init__(self, X, connectivity, n_nodes, n_features, n_samples):
        self.X = X
        self.A = [] 
        # Heap of possible cluster merges, sorted by heir distances
        self.distances = [] 
        # Distances between two clusters (represented by their root node indices) 
        self.distanceDict = {} 
        # Mapping from a pair of clusters to the set of node pairs which
        # connect these two clusters
        self.clusterConnectingNodes = defaultdict(set)
        # Mapping from a a cluster to a mapping from nodes in this cluster
        # to their distance to the most distant node in the cluster
        self.maximalDistanceFromNodeInCluster = defaultdict(dict)  
        # Determine distances between all connected nodes 
        for ind1, row in enumerate(connectivity.rows):
            self.A.append(row)
            self.maximalDistanceFromNodeInCluster[ind1][ind1] = 0.0
            for ind2 in row:
                self.clusterConnectingNodes[(ind1, ind2)] = set([(ind1, ind2)])
                self.clusterConnectingNodes[(ind2, ind1)] = set([(ind1, ind2)])
                # Compute distance between two connected nodes
                dist = np.linalg.norm(X[ind1] - X[ind2])
                self.distances.append((dist, ind1, ind2))
                self.distanceDict[(ind1, ind2)] = dist
                self.distanceDict[(ind2, ind1)] = dist
                
        # Enforce symmetry of A
        for ind1, row in enumerate(connectivity.rows):
            for ind2 in row:
                if ind1 not in self.A[ind2]:
                    self.A[ind2].append(ind1)
        
        heapq.heapify(self.distances)

    def update(self, child_node1, child_node2, parent_node, parent):
        # Update maximalDistanceFromNodeInCluster for new cluster "parent_node" 
        for node in self.maximalDistanceFromNodeInCluster[child_node1].keys():
            minInterClusterDist = np.inf
            for (connecting_node1, connecting_node2) in self.clusterConnectingNodes[(child_node1, child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if connecting_node1 not in self.maximalDistanceFromNodeInCluster[child_node1]:
                    connecting_node1, connecting_node2 = connecting_node2, connecting_node1
                    
                interClusterDist = \
                    self.maximalDistanceFromNodeInCluster[child_node2][connecting_node2] \
                        + self.distanceDict[(connecting_node1, connecting_node2)] \
                        + np.linalg.norm(self.X[node] - self.X[connecting_node1]) # TODO: graph distance instead of euclidean
                minInterClusterDist = min(minInterClusterDist, interClusterDist)
                      
            self.maximalDistanceFromNodeInCluster[parent_node][node] = \
                max(self.maximalDistanceFromNodeInCluster[child_node1][node], # intra cluster dist
                    minInterClusterDist)
        for node in self.maximalDistanceFromNodeInCluster[child_node2].keys():
            minInterClusterDist = np.inf
            for (connecting_node1, connecting_node2) in self.clusterConnectingNodes[(child_node1, child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if connecting_node1 not in self.maximalDistanceFromNodeInCluster[child_node1]:
                    connecting_node1, connecting_node2 = connecting_node2, connecting_node1
                    
                interClusterDist = \
                    self.maximalDistanceFromNodeInCluster[child_node1][connecting_node1] \
                        + self.distanceDict[(connecting_node1, connecting_node2)] \
                        + np.linalg.norm(self.X[node] - self.X[connecting_node2]) # TODO: graph distance instead of euclidean
                minInterClusterDist = min(minInterClusterDist, interClusterDist)
                      
            self.maximalDistanceFromNodeInCluster[parent_node][node] = \
                max(self.maximalDistanceFromNodeInCluster[child_node2][node], # intra cluster dist
                    minInterClusterDist)
        
        # Determine all other clusters that are connected to one of the child
        # clusters. These cluster will also be connected to the parent cluster
        coord_col = []
        visited = set([parent_node])
        for l in set(self.A[child_node1]).union(self.A[child_node2]):
            while parent[l] != l:
                l = parent[l]
            if l not in visited:
                visited.add(l)
                coord_col.append(l)
                self.A[l].append(parent_node)
        self.A.append(coord_col)
        
        # Determine for all connected clusters the distance to the newly formed
        # cluster
        for l in coord_col:
            self.clusterConnectingNodes[(parent_node, l)] = \
                 self.clusterConnectingNodes[(child_node1, l)].union(self.clusterConnectingNodes[(child_node2, l)])
            self.clusterConnectingNodes[(l, parent_node)] = \
                            self.clusterConnectingNodes[(parent_node, l)]
            # Find the distance between pair of nodes that are most distant in l and parent_node
            interClusterDist = np.inf
            for (connecting_node1, connecting_node2) in self.clusterConnectingNodes[(l, parent_node)]:
                # Canonical ordering (connecting node 1 in l)
                if connecting_node1 not in self.maximalDistanceFromNodeInCluster[l]:
                    connecting_node1, connecting_node2 = connecting_node2, connecting_node1
                dist = self.maximalDistanceFromNodeInCluster[l][connecting_node1] \
                        + self.distanceDict[(connecting_node1, connecting_node2)] \
                        + self.maximalDistanceFromNodeInCluster[parent_node][connecting_node2]
                interClusterDist = min(interClusterDist, dist)
            
            self.distanceDict[(l, parent_node)] = interClusterDist
            self.distanceDict[(parent_node, l)] = interClusterDist
            heapq.heappush(self.distances, (interClusterDist, parent_node, l))
            
    def hasMoreCandidates(self):
        return len(self.distances) > 0
    
    def fetchCandidate(self):
        return heapq.heappop(self.distances)
    
    def computeDistance(self, i, j):
        return self.distanceDict[(i, j)]

            
def ward_tree(X, connectivity=None, n_components=None, return_inertias=False, 
              merge_replay=[], DistanceClass='CompleteLinkageDistance', copy=True):
    """Hierarchical clustering based on a Feature matrix.

    The inertia matrix uses a Heapq-based representation.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neigbhoring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.
        Defaut is None, i.e, the ward algorithm is unstructured.

    n_components : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the connectivity matrix.
        
    return_inertias : bool (optional)
        If true, the inertias are returned additionally

    copy : bool (optional)
        Make a copy of connectivity or work inplace. If connectivity
        is not of LIL type there will be a copy in any case.

    Returns
    -------
    children : list of pairs. Length of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    n_components : sparse matrix.
        The number of connected components in the graph.

    n_leaves : int
        The number of leaves in the tree
        
    heights :  list of floats. Length n_nodes
        The inertias associated to the tree's nodes. Only returned, if 
        return_inertias is True
    """
    X = np.asanyarray(X)
    n_samples, n_features = X.shape
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))

    try:
        DistanceClass = eval(DistanceClass)
    except:
        raise Exception("Unknown distance class %s" % DistanceClass)

    # Compute the number of nodes
    if connectivity is not None:
        if n_components is None:
            n_components, _ = cs_graph_components(connectivity)
        if n_components > 1:
            warnings.warn("the number of connected components of the"
            " connectivity matrix is %d > 1. The tree will be stopped early."
            % n_components)
    else:
        out = hierarchy.ward(X)
        children_ = out[:, :2].astype(np.int)
        return children_, 1, n_samples, None

    n_nodes = 2 * n_samples - 1

    if (connectivity.shape[0] != n_samples or
        connectivity.shape[1] != n_samples):
        raise ValueError('Wrong shape for connectivity matrix: %s '
                         'when X is %s' % (connectivity.shape, X.shape))
    # convert connectivity matrix to LIL eventually with a copy
    if sparse.isspmatrix_lil(connectivity) and copy:
        connectivity = connectivity.copy()
    else:
        connectivity = connectivity.tolil()

    # Remove diagonal from connectivity matrix
    connectivity.setdiag(np.zeros(connectivity.shape[0]))
        
    # Compute distances between connected nodes
    distance_class = DistanceClass(X, connectivity, n_nodes, n_features, n_samples)

    # prepare the main fields
    parent = np.arange(n_nodes, dtype=np.int)
    heights = np.zeros(n_nodes)
    open_nodes = np.ones(n_nodes, dtype=bool)
    children = []
    merges = []
    
    reserved_indices = set([]) # The indices which are contained in the replayed merges
    for (i_, j_), k_ in merge_replay:
        reserved_indices.add(i_)
        reserved_indices.add(j_)
    for (i_, j_), k_ in merge_replay:
        if k_ in reserved_indices:
            reserved_indices.remove(k_) # Add this later since we don't know the new index yet
    
    index_mapping = dict(zip(reserved_indices, reserved_indices))    
    # recursive merge loop
    for k in range(n_samples, n_nodes):
        # Fetch merge that will be reapplied next (if any)
        if len(merge_replay) > 0:
            (i_, j_), k_ = merge_replay[0]
            merge_replay = merge_replay[1:]
             # Associate indices
            i = index_mapping[i_]
            j = index_mapping[j_]
            # Compute merge distance
            merge_distance = distance_class.computeDistance(i,j)
            index_mapping[k_] = k
#            print "Reapply", merge_distance, i, j, k
        else: # No merges to be reapplied left
            # Identify the merge that will be applied next
            merge_distance = np.inf
            while distance_class.hasMoreCandidates():
                merge_distance, i, j = distance_class.fetchCandidate()
                if open_nodes[i] and open_nodes[j]:
                    break
                if not distance_class.hasMoreCandidates():
                    merge_distance = np.inf
                    break
            if not distance_class.hasMoreCandidates() and merge_distance == np.inf: 
                # Merge unconnected components with height infinity
                i = k - 1
                desc = set(_hc_get_descendent([i], np.array(children), n_samples, 
                                              add_intermediate_nodes=True))
                for j in range(k-2, 0, -1):
                    if j not in desc:
                        break           
            # Update index mapping
            if i in index_mapping.values():
                i_ = index_mapping.keys()[index_mapping.values().index(i)]
                index_mapping[i_] = k
            if j in index_mapping.values():
                j_ = index_mapping.keys()[index_mapping.values().index(j)]
                index_mapping[j_] = k  
#            print "Novel", merge_distance, i, j, k

            
        parent[i], parent[j], heights[k] = k, k, merge_distance
        merges.append(((i, j), k))
        children.append([i, j])
        open_nodes[i], open_nodes[j] = False, False
        
        # Add new possible merges and their distances to heap
        distance_class.update(i, j, k, parent)
                    

    # Separate leaves in children (empty lists up to now)
    n_leaves = n_samples
    children = np.array(children)  # return numpy array for efficient caching

    if return_inertias:
        return children, n_components, n_leaves, merges, heights
    else:
        return children, n_components, n_leaves, merges


###############################################################################
# Functions for cutting  hierarchical clustering tree

def _hc_get_descendent(ind, children, n_leaves, add_intermediate_nodes=False):
    """Function returning all the descendent leaves of a set of nodes.

    Parameters
    ----------
    ind : list of int
        A list that indicates the nodes for which we want the descendents.

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        This is not defined for leaves.

    n_leaves : int
        Number of leaves.

    Return
    ------
    descendent : list of int
    """
    descendent = []
    while len(ind) != 0:
        i = ind.pop()
        if i < n_leaves:
            descendent.append(i)
        else:
            if add_intermediate_nodes:
                descendent.append(i)
            ci = children[i - n_leaves]
            ind.extend((ci[0], ci[1]))
    return descendent


def _hc_cut(n_clusters, children, n_leaves):
    """Function cutting the ward tree for a given number of clusters.

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters to form.

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        Leaves have empty list of children and are not stored.

    n_leaves : int
        Number of leaves of the tree.

    Return
    ------
    labels : array [n_points]
        cluster labels for each point

    """
    nodes = [np.max(children[-1]) + 1]
    for i in range(n_clusters - 1):
        nodes.extend(children[np.max(nodes) - n_leaves])
        nodes.remove(np.max(nodes))
    labels = np.zeros(n_leaves, dtype=np.int)
    for i, node in enumerate(nodes):
        labels[_hc_get_descendent([node], children, n_leaves)] = i
    return labels

def _hc_cut_inertia(max_inertia, children, n_leaves, inertias):
    """Function cutting the ward tree for a given maximal inertia.

    Parameters
    ----------
    max_inertia : float
        The maximal inertia a cluster is allowed to have. Determines indirectly
        how many clusters are formed. 

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        Leaves have empty list of children and are not stored.

    n_leaves : int
        Number of leaves of the tree.
        
    inertias :  list of floats. Length of n_nodes
        The inertias associated to the tree's nodes.

    Return
    ------
    labels : array
        cluster labels for each point

    """
    open_nodes = [len(inertias)-1] # root of the tree
    cluster_roots = [] 
    while open_nodes != []:
        node = open_nodes[0]
        open_nodes = open_nodes[1:]
        if inertias[node] <= max_inertia: 
            # This tree node is the root of a cluster
            cluster_roots.append(node)
        else:
            # Tree node induces subtree with too large inertia; split it
            open_nodes.extend(children[node - n_leaves])
            
    labels = np.zeros(n_leaves, dtype=np.int)
    for i, node in enumerate(cluster_roots):
        labels[_hc_get_descendent([node], children, n_leaves)] = i
    return labels


###############################################################################
# Class for Ward hierarchical clustering

class Ward(BaseEstimator):
    """Ward hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray or None
        The number of clusters to find.
        
    max_inertia : float or None
        If this value is not None, the max_inertia is used to determine the
        number of returned clusters. In this case, the n_clusters parameter
        is agnored and may be None.

    connectivity : sparse matrix.
        Connectivity matrix. Defines for each sample the neigbhoring
        samples following a given structure of the data.
        Defaut is None, i.e, the hiearchical clustering algorithm is
        unstructured.

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    copy : bool
        Copy the connectivity matrix or work inplace.

    n_components : int (optional)
        The number of connected components in the graph defined by the
        connectivity matrix. If not set, it is estimated.

    Methods
    -------
    fit:
        Compute the clustering

    Attributes
    ----------
    children_ : array-like, shape = [n_nodes, 2]
        List of the children of each nodes.
        Leaves of the tree do not appear.

    labels_ : array [n_points]
        cluster labels for each point

    n_leaves_ : int
        Number of leaves in the hiearchical tree.

    """

    def __init__(self, n_clusters=2, max_inertia=None,
                 memory=Memory(cachedir=None, verbose=0), connectivity=None,
                 copy=True, n_components=None):
        self.n_clusters = n_clusters
        self.max_inertia = max_inertia
        self.memory = memory
        self.copy = copy
        self.n_components = n_components
        self.connectivity = connectivity

    def fit(self, X, merge_replay=[]):
        """Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The samples a.k.a. observations.

        Returns
        -------
        self
        """
        memory = self.memory
        if isinstance(memory, basestring):
            memory = Memory(cachedir=memory)

        # Construct the tree
        if self.max_inertia is None:
            children, n_components, n_leaves, self.merges = \
                    memory.cache(ward_tree)(X, self.connectivity,
                                            n_components=self.n_components,
                                            merge_replay=merge_replay,
                                            copy=self.copy)
            # Cut the tree based on number of desired clusters
            self.labels_ = _hc_cut(self.n_clusters, children, n_leaves)
        else:
            children, n_components, n_leaves, self.merges, inertias = \
                    memory.cache(ward_tree)(X, self.connectivity,
                                            n_components=self.n_components,
                                            return_inertias=True, 
                                            merge_replay=merge_replay,
                                            copy=self.copy)
            # Cut the tree based on maximally allowed inertia
            self.labels_ = _hc_cut_inertia(self.max_inertia, children, n_leaves,
                                           inertias)
            
        return self

###############################################################################
# Ward-based feature agglomeration

class WardAgglomeration(AgglomerationTransform, Ward):
    """Feature agglomeration based on Ward hierarchical clustering

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters.

    connectivity : sparse matrix
        connectivity matrix. Defines for each feature the neigbhoring
        features following a given structure of the data.
        Defaut is None, i.e, the hiearchical agglomeration algorithm is
        unstructured.

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    copy : bool
        Copy the connectivity matrix or work inplace.

    n_components : int (optional)
        The number of connected components in the graph defined by the
        connectivity matrix. If not set, it is estimated.

    Methods
    -------
    fit:
        Compute the clustering of features

    Attributes
    ----------
    children_ : array-like, shape = [n_nodes, 2]
        List of the children of each nodes.
        Leaves of the tree do not appear.

    labels_ : array [n_points]
        cluster labels for each point

    n_leaves_ : int
        Number of leaves in the hiearchical tree.

    """

    def fit(self, X, y=None, **params):
        """Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The data

        Returns
        -------
        self
        """
        return Ward.fit(self, X.T, **params)
