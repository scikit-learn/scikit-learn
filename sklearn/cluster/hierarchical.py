"""Hierarchical Agglomerative Clustering

These routines perform some hierarchical agglomerative clustering of some
input data. Currently, only complete-linkage and ward's criterion are
implemented as linkage criteria.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux, Jan Hendrik Metzen
License: BSD 3 clause

.. TODO:: Consider nearest-neighbor chain algorithm for creating dendrogram?
 (https://secure.wikimedia.org/wikipedia/en/wiki/Nearest-neighbor_chain_algorithm)
"""

import warnings

import numpy as np
from scipy import sparse
from scipy.cluster import hierarchy

from ..base import BaseEstimator
from ..utils._csgraph import cs_graph_components
from ..externals.joblib import Memory

from .linkage import WardsLinkage
from ._feature_agglomeration import AgglomerationTransform

class Dendrogram(object):
    """
    
    heights...
    
    Attributes
    ----------
    heights :  list of floats. Length of n_nodes
            The heights associated to the tree's nodes.
    
    children : list of pairs. Length of n_nodes
        List of the children of each nodes. This is not defined for leaves.
        
    n_leaves : int
            Number of leaves.
    """

    def __init__(self, n_samples, n_components):
        self.n_samples = n_samples 
        self.n_components = n_components
        
        n_nodes = 2 * self.n_samples - 1
        self.parent = np.arange(n_nodes, dtype=np.int)
        self.heights = np.zeros(n_nodes)
        self.children = []
        
        self.merges = [] # TODO: Get rid of this
        
     
    def merge(self, i, j, k, merge_distance):
        self.parent[i] = self.parent[j] = k
        self.heights[k] = merge_distance
        self.children.append([i, j])
        self.merges.append(((i, j), k))
        
        
    # Functions for cutting  hierarchical clustering tree
    def _get_descendent(self, ind, add_intermediate_nodes=False):
        """Function returning all the descendent leaves of a set of nodes.
    
        Parameters
        ----------
        ind : list of int
            A list that indicates the nodes for which we want the descendents.
            
        add_intermediate_nodes : bool
            If true, leaves and inenr nodes in the subtree are returnd, 
            otherwise only the leaves. Defaults to False.
    
        Return
        ------
        descendent : list of int
        """
        descendent = []
        while len(ind) != 0:
            i = ind.pop()
            if i < self.n_samples: #its a leaf
                descendent.append(i)
            else: # inner node, go to children
                if add_intermediate_nodes:
                    descendent.append(i)
                ci = self.children[i - self.n_samples]
                ind.extend((ci[0], ci[1]))
        return descendent
    
    
    def _cut(self, n_clusters):
        """Function cutting the dendrogram for a given number of clusters.
    
        Parameters
        ----------
        n_clusters : int or ndarray
            The number of clusters to form.
    
        Return
        ------
        labels : array [n_points]
            cluster labels for each point
    
        """
        nodes = [np.max(self.children[-1]) + 1]
        for i in xrange(n_clusters - 1):
            nodes.extend(self.children[np.max(nodes) - self.n_samples])
            nodes.remove(np.max(nodes))
        labels = np.zeros(self.n_samples, dtype=np.int)
        for i, node in enumerate(nodes):
            labels[self._get_descendent([node])] = i
        return labels
    
    
    def _cut_height(self, max_height):
        """ Function cutting the dendrogram for a given maximal height.
        
        Parameters
        ----------
        max_height : float
            The maximal height a tree node is allowed to have to form a single
            cluster. Determines indirectly how many clusters are formed. 
    
        Return
        ------
        labels : array
            cluster labels for each point
    
        """
        open_nodes = [len(self.heights) - 1]  # root of the tree
        cluster_roots = [] 
        while open_nodes != []:
            node = open_nodes[0]
            open_nodes = open_nodes[1:]
            if self.heights[node] <= max_height: 
                # This tree node is the root of a cluster
                cluster_roots.append(node)
            else:
                # Tree node induces subtree with too large height; split it
                open_nodes.extend(self.children[node - self.n_samples])
                
        labels = np.zeros(self.n_samples, dtype=np.int)
        for i, node in enumerate(cluster_roots):
            labels[self._get_descendent([node])] = i
        return labels

###############################################################################
              
def create_dendrogram(X, connectivity, n_components=None, merge_replay=[], 
                      linkage_criterion=WardsLinkage, linkage_kwargs={}, 
                      copy=True):
    """Hierarchical clustering based on a Feature matrix.

    The height matrix uses a Heapq-based representation.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neighboring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.
        Default is None, i.e, full connectivity is assumed

    n_components : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the connectivity matrix.

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
        The heights associated to the tree's nodes.
    """
    X = np.asarray(X)
    n_samples, n_features = X.shape
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))

    if n_components is None:
        n_components, _ = cs_graph_components(connectivity)
    if n_components > 1:
        warnings.warn("the number of connected components of the"
        " connectivity matrix is %d > 1. The tree will be stopped early."
        % n_components)

    # Compute the number of nodes
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
    linkage = linkage_criterion(X, connectivity, **linkage_kwargs)

    # Create the dendrogram
    dendrogram = Dendrogram(n_samples, n_components)
    
    open_nodes = np.ones(n_nodes, dtype=bool)
    
    # The indices which are contained in the replayed merges
    reserved_indices = set([])
    for (i_, j_), k_ in merge_replay:
        reserved_indices.add(i_)
        reserved_indices.add(j_)
    for (i_, j_), k_ in merge_replay:
        if k_ in reserved_indices:
             # Add k_ later since we don't know the new index yet
            reserved_indices.remove(k_)
    
    index_mapping = dict(zip(reserved_indices, reserved_indices))
    
    augmentation = None
    # recursive merge loop
    replaying = len(merge_replay) > 0
    for k in xrange(n_samples, n_nodes):
        # Fetch merge that will be reapplied next (if any)
        if len(merge_replay) > 0:
            (i_, j_), k_ = merge_replay[0]
            merge_replay = merge_replay[1:]
             # Associate indices
            i = index_mapping[i_]
            j = index_mapping[j_]
            # Compute merge distance
            merge_distance = linkage.get_cluster_distance(i, j)
            index_mapping[k_] = k
#            print "Reapply", merge_distance, i, j, k
        else:  # No merges to be reapplied left
            if augmentation: 
                # There are nodes that can be added to cluster created during
                # replaying 
                merge_distance, i, j = augmentation
                augmentation = None
#                print "Augmentation", merge_distance, i, j, k
            else:
                # Identify the merge that will be applied next using the 
                # standard method
                merge_distance = np.inf
                while linkage.has_more_candidates():
                    merge_distance, i, j = linkage.fetch_candidate()
                    if open_nodes[i] and open_nodes[j]:
                        break
                    if not linkage.has_more_candidates():
                        merge_distance = np.inf
                        break
                if not linkage.has_more_candidates() \
                                        and merge_distance == np.inf: 
                    # Merge unconnected components with height infinity
                    i = k - 1
                    desc = set(_hc_get_descendent([i], np.array(children), 
                                                  n_samples, 
                                                  add_intermediate_nodes=True))
                    for j in xrange(k - 2, 0, -1):
                        if j not in desc:
                            break      
#                print "Novel", merge_distance, i, j, k     
            # Update index mapping
            if i in index_mapping.values():
                i_ = index_mapping.keys()[index_mapping.values().index(i)]
                index_mapping[i_] = k
            if j in index_mapping.values():
                j_ = index_mapping.keys()[index_mapping.values().index(j)]
                index_mapping[j_] = k  

        dendrogram.merge(i, j, k, merge_distance)
        
        open_nodes[i], open_nodes[j] = False, False
        
        # Add new possible merges and their distances to heap
        linkage.update(i, j, k, dendrogram.parent)
        
        # If we are finished with replaying
        if replaying and len(merge_replay) == 0:
            # Check for possible augmentations of the clusters created during 
            # replaying
            possible_augmentations = []
            for cluster_root in linkage.get_nontrivial_clusters():
                for node in linkage.get_nodes_connected_to_clusters(
                                                            [cluster_root]):
                    if not open_nodes[node]:
                        continue
                    dist = linkage.get_cluster_distance(node, cluster_root)
                    if dist < heights[cluster_root]:
                        possible_augmentations.append(
                                                (dist - heights[cluster_root],
                                                 heights[cluster_root],
                                                 node, cluster_root))
            if len(possible_augmentations) > 0:
                augmentation = min(possible_augmentations)[1:]
            else:
                replaying = False
                
    return dendrogram

def ward_tree(X, connectivity=None, n_components=None, copy=True):
    if connectivity is not None:    
        dendrogram = create_dendrogram(X, connectivity, n_components, 
                                       merge_replay=[], 
                                       linkage_criterion=WardsLinkage, copy=True) 
        return dendrogram.children, dendrogram.n_components, dendrogram.n_samples
    else:
        out = hierarchy.ward(X)
        children_ = out[:, :2].astype(np.int)
        return children_, 1, X.shape[0]

###############################################################################

class HierarchicalClustering(BaseEstimator):
    """Hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray or None
        The number of clusters to find.
        
    max_height : float or None
        If this value is not None, the max_height is used to determine the
        number of returned clusters. In this case, the n_clusters parameter
        is ignored and may be None.

    connectivity : sparse matrix.
        Connectivity matrix. Defines for each sample the neigbhoring
        samples following a given structure of the data.
        Default is None, i.e, the hiearchical clustering algorithm is
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

    def __init__(self, n_clusters=None, max_height=None,
                 memory=Memory(cachedir=None, verbose=0), connectivity=None,
                 copy=True, n_components=None, linkage_criterion=WardsLinkage,
                 linkage_kwargs={}):
        self.n_clusters = n_clusters
        self.n_components = n_components
        self.max_height = max_height
                
        self.connectivity = connectivity
        
        self.linkage_criterion = linkage_criterion
        self.linkage_kwargs = linkage_kwargs
        
        self.memory = memory
        self.copy = copy
        
        self.merges = []
        
    def fit(self, X):
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
        dendrogram = \
            memory.cache(create_dendrogram)(X, self.connectivity,
                                            n_components=self.n_components,
                                            merge_replay=self.merges,
                                            linkage_criterion=self.linkage_criterion,
                                            linkage_kwargs=self.linkage_kwargs,
                                            copy=self.copy)
        # Cut the tree ... 
        if self.max_height is None:
            # based on number of desired clusters
            self.labels_ = dendrogram._cut(self.n_clusters)
        else:
            # based on maximally allowed height
            self.labels_ = dendrogram._cut_height(self.max_height)

        self.merges = dendrogram.merges           
        # Undo merges in the uppermost part of the tree that have been 
        # "cut off" by _hc_cut_height or _hc_cut
        if len(np.unique(self.labels_)) > 1:
            self.merges = self.merges[:-len(np.unique(self.labels_)) + 1]
            
        return self


###############################################################################
class Ward(HierarchicalClustering):
    
    def __init__(self, *args, **kwargs):
        super(Ward, self).__init__(*args, linkage_criterion=WardsLinkage,
                                   **kwargs)


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
        Default is None, i.e, the hiearchical agglomeration algorithm is
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
