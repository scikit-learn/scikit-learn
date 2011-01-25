"""
These routines perform some hierachical agglomerative clustering of some input
data. Currently, only Ward's algorithm is implemented.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux
License: BSD 3 clause
"""
import heapq as heapq
import warnings

import numpy as np
from scipy import sparse

from ..base import BaseEstimator
from ..utils._csgraph import cs_graph_components
from ..externals.joblib import Memory

from . import _inertia
from ._feature_agglomeration import AgglomerationTransform

###############################################################################
# Ward's algorithm

def ward_tree(X, connectivity=None, n_comp=None, copy=True):
    """Ward clustering based on a Feature matrix. Heapq-based representation
    of the inertia matrix.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neigbhoring samples
        following a given structure of the data.
        Defaut is None, i.e, the ward algorithm is unstructured.

    n_comp : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the connectivity matrix.

    copy : bool (optional)
        Make a copy of connectivity or work inplace. If connectivity
        is not of LIL type there will be a copy in any case.

    Returns
    -------
    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

    children : list of pairs. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    heights : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    n_comp : sparse matrix.
        The number of connected components in the graph.

    """
    X = np.asanyarray(X)
    n_samples, n_features = X.shape
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))

    # Compute the number of nodes
    if connectivity is not None:
        if n_comp is None:
            n_comp, _ = cs_graph_components(connectivity)
    else:
        n_comp = 1

    if n_comp > 1:
        warnings.warn("the number of connected components of the"
        " connectivity matrix is %d > 1. The tree will be stopped early."
        % n_comp)

    n_nodes = 2 * n_samples - n_comp

    # convert connectivity matrix to LIL eventually with a copy
    if connectivity is None:
        connectivity = np.ones([n_samples, n_samples])
        connectivity.flat[::n_samples+1] = 0 # set diagonal to 0
        connectivity = sparse.lil_matrix(connectivity)
        n_nodes = 2 * n_samples - 1
    else:
        if sparse.isspmatrix_lil(connectivity) and copy:
            connectivity = connectivity.copy()
        else:
            connectivity = connectivity.tolil()

    # Remove diagonal from connectivity matrix
    connectivity.setdiag(np.zeros(connectivity.shape[0]))

    # build moments as a list
    moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features)),
                np.zeros((n_nodes, n_features))]
    moments[0][:n_samples] = 1
    moments[1][:n_samples] = X
    moments[2][:n_samples] = X ** 2

    # create inertia matrix
    cord_row = []
    cord_col = []
    B = []
    for ind, row in enumerate(connectivity.rows):
        cord_row.extend(list(ind * np.ones(len(row), dtype=int)))
        cord_col.extend(row)
        B.append(row)
    A = B
    inertia = np.zeros(len(cord_row), dtype=np.float)
    _inertia.compute_inertia(moments[0][cord_row], moments[0][cord_col], \
                             moments[1][cord_row], moments[1][cord_col], \
                             moments[2][cord_row], moments[2][cord_col], \
                             inertia)
    inertia = zip(inertia, cord_row, cord_col)
    heapq.heapify(inertia)

    # prepare the main fields
    parent = np.arange(n_nodes, dtype=np.int)
    heights = np.zeros(n_nodes)
    used_node = np.ones(n_nodes, dtype=bool)
    children = []
    for k in range(n_samples):
        children.append([])

    # recursive merge loop
    for k in range(n_samples, n_nodes):

        # identify the merge
        while True:
            node = heapq.heappop(inertia)
            i, j = node[1], node[2]
            if used_node[i] and used_node[j]:
                break
        parent[i], parent[j], heights[k] = k, k, node[0]
        children.append([i, j])
        used_node[i], used_node[j] = False, False

        # update the moments
        for p in range(3):
            moments[p][k] = moments[p][i] + moments[p][j]

        # update the structure matrix A and the inertia matrix
        cord_col = []
        for l in set(A[i]).union(A[j]):
            if parent[l] == l:
                cord_col.append(l)
                A[l].append(k)
        A.append(cord_col)
        cord_row = len(cord_col) * [k]
        ini = np.zeros(len(cord_row), dtype=np.float)
        _inertia.compute_inertia(moments[0][cord_row], moments[0][cord_col], \
                             moments[1][cord_row], moments[1][cord_col], \
                             moments[2][cord_row], moments[2][cord_col], \
                             ini)
        ini = zip(ini, cord_row, cord_col)
        for tupl in ini:
            heapq.heappush(inertia, tupl)

    # Separate leaves in children (empty lists up to now)
    n_leaves = 0
    while len(children) > 0 and len(children[0]) == 0:
        n_leaves += 1
        children.pop(0)
    children = np.array(children) # return as numpy array for efficient caching

    return parent, children, heights, n_comp, n_leaves


###############################################################################
# Functions for cutting  hierarchical clustering tree

def _hc_get_descendent(ind, children, n_leaves):
    """
    Function returning all the descendent leaves of a set of nodes in the tree.

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
            ci = children[i - n_leaves]
            ind.extend((ci[0], ci[1]))
    return descendent


def _hc_cut(n_clusters, parent, children, n_leaves, heights):
    """
    Function cutting the ward tree for a given number of clusters.

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters to form.

    parent : array-like, shape = [n_nodes]
        Int. Gives the parent node for each node, i.e. parent[i] is the
        parent node of the node i. The last value of parent is the
        root node, that is its self parent, so the last value is taken
        3 times in the array.
        The n_nodes is equal at  (2*n_samples - 1), and takes into
        account the nb_samples leaves, and the unique root.

    children : list of pairs. Lenght of n_nodes
        List of the children of each nodes.
        Leaves have empty list of children and are not stored.

    n_leaves : int
        Number of leaves of the tree.

    heights : array-like, shape = [n_nodes]
        Float. Gives the inertia of the created nodes. The n_samples first
        values of the array are 0, and thus the values are positive (or
        null) and are ranked in an increasing order.

    Return
    ------
    labels_ : array [n_points]
        cluster labels for each point

    active_nodes : list of int
                index of the nodes kept for the labeling
    """
    parent = parent[:-1]
    heights = heights[:-1]
    active_nodes = [len(parent)]
    node_to_cut = active_nodes[0]
    for i in range(n_clusters - 1):
        if np.sum(parent == node_to_cut) != 0:
            active_nodes.append(np.where(parent == node_to_cut)[0][0])
            active_nodes.append(np.where(parent == node_to_cut)[0][1])
            active_nodes.remove(node_to_cut)
        else:
            active_nodes.append(node_to_cut)
        node_to_cut = active_nodes[np.argmax(heights[active_nodes])]
    label = np.zeros(n_leaves)
    for j in active_nodes[:n_clusters - 1]:
        ind = [j]
        label[_hc_get_descendent(ind, children, n_leaves)] = np.max(label) + 1
    return label, active_nodes

###############################################################################
# Class for Ward hierarchical clustering

class Ward(BaseEstimator):
    """
    Class for Ward hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters.

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    Methods
    -------
    fit:
        Compute the clustering

    Attributes
    ----------
    parent_ : array-like, shape = [n_nodes]
        Int. Gives the parent node for each node, i.e. parent[i] is the
        parent node of the node i. The last value of parent is the
        root node, that is its self parent, so the last value is taken
        3 times in the array.
        The n_nodes is equal at  (2*n_samples - 1), and takes into
        account the nb_samples leaves, and the unique root.

    children_ : list of pairs. Length of n_nodes
        List of the children of each nodes.
        Leaves of the tree have empty list of children.

    heights_ : array-like, shape = [n_nodes]
        Gives the inertia of the created nodes. The n_samples first
        values of the array are 0, and thus the values are positive (or
        null) and are ranked in an increasing order.

    labels_ : array [n_points]
        cluster labels for each point

    n_leaves_ : int
        Number of leaves in the hiearchical tree.

    Return
    ------
    self
    """

    def __init__(self, n_clusters=2, memory=Memory(cachedir=None, verbose=0),
                 connectivity=None, copy=True, n_comp=None):
        """
        connectivity : sparse matrix.
            connectivity matrix. Defines for each sample the neigbhoring
            samples following a given structure of the data.
            Defaut is None, i.e, the hiearchical clustering algorithm is
            unstructured.
        """
        self.n_clusters = n_clusters
        self.memory = memory
        self.copy = copy
        self.n_comp = n_comp
        self.connectivity = connectivity

    def fit(self, X, **params):
        """
        Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            A M by N array of M observations in N dimensions or a length
            M array of M one-dimensional observations.

        Returns
        -------
        self
        """
        self._set_params(**params)

        memory = self.memory
        if isinstance(memory, basestring):
            memory = Memory(cachedir=memory)

        # Construct the tree
        self.parent_, self.children_, self.heights_, self.n_comp, \
            self.n_leaves_ = memory.cache(ward_tree)(X, self.connectivity,
                                          n_comp=self.n_comp, copy=self.copy)

        # Cut the tree
        self.labels_, self.active_nodes_ = _hc_cut(self.n_clusters,
                                self.parent_, self.children_, self.n_leaves_,
                                self.heights_)
        return self

###############################################################################
# Ward-based feature agglomeration

class WardAgglomeration(AgglomerationTransform, Ward):
    """Feature agglomeration based on Ward hierarchical clustering

    XXX
    """

    def fit(self, X, y=None, **params):
        return Ward.fit(self, X.T, n_clusters=self.n_clusters, **params)
