"""Hierarchical Agglomerative Clustering

These routines perform some hierachical agglomerative clustering of some
input data. Currently, only Ward's algorithm is implemented.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux
License: BSD 3 clause
"""
import heapq
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
# Ward's algorithm

def ward_tree(X, connectivity=None, n_components=None, return_inertias=False, 
              merge_replay=[], copy=True):
    """Ward clustering based on a Feature matrix.

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
        return children_, 1, n_samples

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
    
    # create inertia matrix
    coord_row = []
    coord_col = []
    A = []
    for ind, row in enumerate(connectivity.rows):
        A.append(row)
        # We keep only the upper triangular for the moments
        # Generator expressions are faster than arrays on the following
        row = [i for i in row if i < ind]
        coord_row.extend(len(row) * [ind, ])
        coord_col.extend(row)

    coord_row = np.array(coord_row, dtype=np.int)
    coord_col = np.array(coord_col, dtype=np.int)

    # build moments as a list
    moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features))]
    moments[0][:n_samples] = 1
    moments[1][:n_samples] = X
    inertia = np.empty(len(coord_row), dtype=np.float)
    _inertia.compute_ward_dist(moments[0], moments[1],
                               coord_row, coord_col, inertia)
    inertia = zip(inertia, coord_row, coord_col)
    heapq.heapify(inertia)    

    # prepare the main fields
    parent = np.arange(n_nodes, dtype=np.int)
    heights = np.zeros(n_nodes)
    used_node = np.ones(n_nodes, dtype=bool)
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
            # Compute inertia of the cluster obtained when merging i and j
            ini = np.empty(len(coord_row), dtype=np.float)
            _inertia.compute_ward_dist(moments[0], moments[1],
                                       np.array([i]), np.array([j]), ini)
            height = ini[0]
            index_mapping[k_] = k
#            print "Reapply", height, i, j, k
        else: # No merges to be reapplied left
            # Identify the merge that will be applied next
            height = np.inf
            while len(inertia) > 0:
                height, i, j = heapq.heappop(inertia)
                if used_node[i] and used_node[j]:
                    break
                if len(inertia) == 0:
                    height = np.inf
                    break
            if len(inertia) == 0 and height == np.inf: 
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
#            print "Novel", height, i, j, k
            
        parent[i], parent[j], heights[k] = k, k, height
        merges.append(((i,j),k))
        children.append([i, j])
        used_node[i], used_node[j] = False, False

        # update the moments
        for p in range(2):
            moments[p][k] = moments[p][i] + moments[p][j]

        # update the structure matrix A and the inertia matrix
        coord_col = []
        visited = set([k])        
        for l in set(A[i]).union(A[j]):
            while parent[l] != l:
               l = parent[l]
            if l not in visited:
               visited.add(l)
               coord_col.append(l)
               A[l].append(k)
        A.append(coord_col)
        coord_col = np.array(coord_col, dtype=np.int)
        coord_row = np.empty_like(coord_col)
        coord_row.fill(k)
        ini = np.empty(len(coord_row), dtype=np.float)

        _inertia.compute_ward_dist(moments[0], moments[1],
                                   coord_row, coord_col, ini)
        
        for tupl in itertools.izip(ini, coord_row, coord_col):
            heapq.heappush(inertia, tupl)

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
    nodes = [np.max(children[-1]) + 1]
    for i in range(len(inertias)):
        if inertias[-(i+1)] <= max_inertia: break
        nodes.extend(children[np.max(nodes) - n_leaves])
        nodes.remove(np.max(nodes))
    labels = np.zeros(n_leaves, dtype=np.int)
    for i, node in enumerate(nodes):
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
