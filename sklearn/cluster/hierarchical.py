"""Hierarchical Agglomerative Clustering

These routines perform hierarchical agglomerative clustering of some
input data. 

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux, Jan Hendrik Metzen
License: BSD 3 clause

.. TODO:: Consider nearest-neighbor chain algorithm for creating dendrogram?
 (https://secure.wikimedia.org/wikipedia/en/wiki/
                 Nearest-neighbor_chain_algorithm)
"""

import numpy as np
from scipy import sparse
from scipy.cluster import hierarchy

from ..base import BaseEstimator
from ..utils._csgraph import cs_graph_components
from ..externals.joblib import Memory

from .linkage import WardsLinkage
from ._feature_agglomeration import AgglomerationTransform


class Dendrogram(object):
    """ Dendrogram that is constructed during hierarchical clustering.
    
    The dendrogram encodes the agglomerative clustering process. The dendrogram
    is a binary tree in which each datapoint of the dataset corresponds to
    one leaf and each inner node corresponds to one cluster that was formed 
    during the agglomerative clustering. The cluster associated to a tree node
    contains exactly those datapoints that are the leaves of the subtree 
    induced by this node.
    
    Each node of the tree has an additional "height" attribute (not to be
    mistaken as the node's height in the tree) which corresponds to the 
    distance of its two child nodes. This distance is computed using one of
    the linkage criteria defined in linkage.py.   
    
    Parameters
    ----------
    n_samples : int
        The number of samples in the dataset to be clustered (i.e. the number
        of leaves of the dendrogram tree)
    
    n_components : int
        The number of connected components in the dataset

    Methods
    ----------
    merge:
        Adds a new tree node corresponding to merging two clusters.
        
    _get_descendent:
        Return all the descendent leaves of a set of nodes.
    
    _cut:
       Cut the dendrogram for a given number of clusters.
       
    _cut_height:
        Cut the dendrogram for a given maximal height.
        
    
    Attributes
    ----------
    n_samples : int
        The number of samples in the dataset to be clustered (i.e. the number
        of leaves of the dendrogram tree)
    
    n_components : int
        The number of connected components in the dataset

    parent : array of floats. Length of n_nodes
        For tree node with index i, this array contains as i-th entry the index
        of its parent. Nodes without parents contain themselves as parents. 
        
    heights :  array of floats. Length of n_nodes
            The heights associated to the tree's nodes.
    
    children : list of pairs. Length of n_nodes
        List of the children of node. This is not defined for leaves.
        
    next_index : int
        The index that a new tree node (resulting from merging) will obtain.
    """

    def __init__(self, n_samples, n_components):
        self.n_samples = n_samples 
        self.n_components = n_components
        
        n_nodes = 2 * self.n_samples - 1
        self.parent = np.arange(n_nodes, dtype=np.int)
        self.heights = np.zeros(n_nodes)
        self.children = []
        
        self.next_index = self.n_samples
        
    def merge(self, child_node_1, child_node_2, merge_distance):
        """ Adds a new tree node corresponding to merging two clusters.
        
        Create tree node *parent_node* that is the parent of *child_node_1* and 
        *child_node_2* and gets as its height the *merge_distance*.
        
        Parameters
        ----------
        child_node_1 : int
            The index of the first child node
            
        child_node_2 : int
            The index of the second child node
        
        merge_distance: float
            The distance of the two child nodes corresponding to the two 
            clusters.
            
        Return
        ------    
        next_index : int
            The index of the new parent node that connects the two child nodes
            
        """
        self.parent[child_node_1] = self.parent[child_node_1] = self.next_index
        self.heights[self.next_index] = merge_distance
        self.children.append([child_node_1, child_node_2])
        self.next_index += 1
        
        return self.next_index
        
    def _get_descendent(self, ind, add_intermediate_nodes=False):
        """ Return all the descendent leaves of a set of nodes.
    
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
            if i < self.n_samples:  # its a leaf
                descendent.append(i)
            else:  # inner node, go to children
                if add_intermediate_nodes:
                    descendent.append(i)
                ci = self.children[i - self.n_samples]
                ind.extend((ci[0], ci[1]))
        return descendent
    
    def _cut(self, n_clusters):
        """ Cut the dendrogram for a given number of clusters.
    
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
        """ Cut the dendrogram for a given maximal height.
        
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


def create_dendrogram(X, connectivity, n_components=None, 
                      linkage_criterion=WardsLinkage, linkage_kwargs={}, 
                      copy=True):
    """ Hierarchical clustering algorithm that creates a dendrogram.

    This is the structured version, that takes into account the topological
    structure between samples. The linkage is Ward's linkage per default;
    however one can use any linkage object that implements the Linkage
    interface (see interface.py)

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neighboring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.

    n_components : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the connectivity matrix.

    linkage_criterion : object implementing the Linkage interface
        The linkage criterion used to determine the distances of two clusters.
        Defaults to Ward's linkage; however one can pass any object 
        implementing the Linkage interface.
        
    linkage_kwargs : dict
        Additional keyword arguments that are directly passed to the __init__
        method of the linkage object.

    copy : bool (optional)
        Make a copy of connectivity or work inplace. If connectivity
        is not of LIL type there will be a copy in any case.

    Returns
    -------
    dendrogram : instance of Dendrogram
        The dendrogram that was formed during hierarchical clustering.
    """
    X = np.asarray(X)
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))
    n_samples, n_features = X.shape
    
    # Sanity checks
    if (connectivity.shape[0] != n_samples or
        connectivity.shape[1] != n_samples):
        raise ValueError('Wrong shape for connectivity matrix: %s '
                         'when X is %s' % (connectivity.shape, X.shape))
        
    # Compute the number of nodes of the tree
    # Binary tree with n leaves has 2n-1 nodes...
    n_nodes = 2 * n_samples - 1

    # convert connectivity matrix to LIL, possibly with a copy
    if sparse.isspmatrix_lil(connectivity) and copy:
        connectivity = connectivity.copy()
    else:
        connectivity = connectivity.tolil()
    
    if n_components is None:
        n_components, _ = cs_graph_components(connectivity)

    # Remove diagonal from connectivity matrix
    connectivity.setdiag(np.zeros(connectivity.shape[0]))        
    
    # Create the dendrogram
    dendrogram = Dendrogram(n_samples, n_components)
    
    # Linkage object that manages the distance computations between clusters
    # distances are updated incrementally during merging of clusters...
    linkage = linkage_criterion(X, connectivity, **linkage_kwargs)
    
    # Array in which open_nodes[i] indicate whether the cluster with index i
    # has nor yet been merged into a larger cluster
    open_nodes = np.ones(n_nodes, dtype=bool)
    
    # Recursive merge loop
    k = n_samples
    while k < n_nodes:
        # Identify the merge that will be applied next
        # This is the merge with minimal distance of two cluster that haven't
        # been merged into a larger cluster yet.
        merge_distance = np.inf
        while linkage.has_more_candidates():
            merge_distance, i, j = linkage.fetch_candidate()
            if open_nodes[i] and open_nodes[j]:
                break
            if not linkage.has_more_candidates():
                merge_distance = np.inf
                break
            
        # Check if we have fully merged all connected components
        if merge_distance == np.inf and not linkage.has_more_candidates(): 
            # Merge unconnected components with height infinity
            i = k - 1
            desc = set(dendrogram._get_descendent([i], 
                                                  add_intermediate_nodes=True))
            for j in xrange(k - 2, 0, -1):
                if j not in desc:
                    break      

        # Add one node to dendrogram tree that is the parent of the two
        # tree nodes i and j. Store the corresponding merge distance as the 
        # nodes' height 
        k = dendrogram.merge(i, j, merge_distance)
        
        # Update linkage object
        linkage.update(i, j, k - 1, dendrogram.parent)
        
        open_nodes[i], open_nodes[j] = False, False
                
    return dendrogram


def ward_tree(X, connectivity=None, n_components=None, copy=True):
    """Hierarchical clustering based on Ward's criterion.

    If connectivity is None, the scipy implementation of Ward's algorithm
    is used. Otherwise, the structured version, that takes into account 
    the topological structure between samples is used.

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
    """
    # TODO: Contained only for backward compatibility. Add a deprecation 
    #       warning?
    if connectivity is not None:    
        dendrogram = create_dendrogram(X, connectivity, n_components, 
                                       linkage_criterion=WardsLinkage, 
                                       copy=True) 
        return dendrogram.children, dendrogram.n_components, \
                dendrogram.n_samples
    else:
        out = hierarchy.ward(X)
        children_ = out[:, :2].astype(np.int)
        return children_, 1, X.shape[0]

###############################################################################


class HierarchicalClustering(BaseEstimator):
    """Hierarchical clustering: constructs a dendrogram and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray or None
        The number of clusters to find. If None, max_height must be specified 
        instead.
        
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
        
    linkage_criterion : object implementing the Linkage interface
        The linkage criterion used to determine the distances of two clusters.
        Defaults to Ward's linkage; however one can pass any object 
        implementing the Linkage interface.
        
    linkage_kwargs : dict
        Additional keyword arguments that are directly passed to the __init__
        method of the linkage object.

    Methods
    -------
    fit:
        Compute the clustering

    Attributes
    ----------
    labels_ : array [n_points]
        cluster labels for each point

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
        
        if self.connectivity is not None:
            # Construct the tree
            dendrogram = \
                memory.cache(create_dendrogram)(X, self.connectivity, 
                                                n_components=self.n_components,
                                                linkage_criterion=\
                                                        self.linkage_criterion,
                                                linkage_kwargs=\
                                                        self.linkage_kwargs,
                                                copy=self.copy)
            # Cut the tree ... 
            if self.max_height is None:
                # based on number of desired clusters
                self.labels_ = dendrogram._cut(self.n_clusters)
            else:
                # based on maximally allowed height
                self.labels_ = dendrogram._cut_height(self.max_height)
        else:  # Fall back to scipy
            assert self.n_clusters is not None, \
                "Unstructured clustering requires the number of clusters to "\
                "be specified explicitly." 
            
            out = hierarchy.ward(X)
            # Put result into a dendrogram and cut it to get labeling
            dendrogram = Dendrogram(X.shape[0], 1)
            dendrogram.children = out[:, :2].astype(np.int)
            self.labels_ = dendrogram._cut(self.n_clusters)
            
        return self


###############################################################################

class Ward(HierarchicalClustering):
    # TODO: Contained only for backward compatibility. Add a deprecation 
    #       warning?
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
    labels_ : array [n_points]
        cluster labels for each point

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
