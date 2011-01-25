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

def ward_tree(X, adjacency_matrix=None, n_comp=None, copy=True):
    """Ward clustering based on a Feature matrix. Heapq-based representation
    of the inertia matrix.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    adjacency_matrix : sparse matrix.
        adjacency matrix. Defines for each sample the neigbhoring samples
        following a given structure of the data.
        Defaut is None, i.e, the ward algorithm is unstructured.

    n_comp : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the adjacency matrix.

    copy : bool (optional)
        Make a copy of adjacency_matrix or work inplace. If adjacency_matrix
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
    if adjacency_matrix is not None:
        if n_comp is None:
            n_comp, _ = cs_graph_components(adjacency_matrix)
    else:
        n_comp = 1

    if n_comp > 1:
        warnings.warn("the number of connected components of the"
        " adjacency matrix is %d > 1. The tree will be stopped early."
        % n_comp)

    n_nodes = 2 * n_samples - n_comp

    # convert adjacency matrix to LIL eventually with a copy
    if adjacency_matrix is None:
        adjacency_matrix = np.ones([n_samples, n_samples])
        adjacency_matrix.flat[::n_samples+1] = 0 # set diagonal to 0
        adjacency_matrix = sparse.lil_matrix(adjacency_matrix)
        n_nodes = 2 * n_samples - 1
    else:
        if sparse.isspmatrix_lil(adjacency_matrix) and copy:
            adjacency_matrix = adjacency_matrix.copy()
        else:
            adjacency_matrix = adjacency_matrix.tolil()

    # Remove diagonal from adjacency matrix
    adjacency_matrix.setdiag(np.zeros(adjacency_matrix.shape[0]))

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
    for ind, row in enumerate(adjacency_matrix.rows):
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

    return parent, children, heights, n_comp


###############################################################################
# Functions for cutting  hierarchical clustering tree

def _hc_get_descendent(ind, children):
    """
    Function returning all the descendent leaves of a set of nodes in the tree.

    Parameters
    ----------
    ind : list of int
          A list that indicates the nodes for which we want the descendents.

    children : list of pairs. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    Return
    ------
    descendent : list of int
    """
    descendent = []
    while len(ind) != 0:
        i = ind.pop()
        ci = children[i]
        if len(ci) == 0:
            descendent.append(i)
        else:
            ind.extend(ci)
    return descendent


def _hc_cut(k, parent, children, heights):
    """
    Function cutting the ward tree for a given number of clusters.

    Parameters
    ----------
    k : int or ndarray
        The number of clusters to form.

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
    for i in range(k - 1):
        if np.sum(parent == node_to_cut) != 0:
            active_nodes.append(np.where(parent == node_to_cut)[0][0])
            active_nodes.append(np.where(parent == node_to_cut)[0][1])
            active_nodes.remove(node_to_cut)
        else:
            active_nodes.append(node_to_cut)
        node_to_cut = active_nodes[np.argmax(heights[active_nodes])]
    label = np.zeros(children.count([]))
    for j in active_nodes[:k - 1]:
        ind = [j]
        label[_hc_get_descendent(ind, children)] = np.max(label) + 1
    return label, active_nodes


###############################################################################
# Display functions for hierarchical clustering

def plot_dendrogram(children, parent, heights, ax=None, active_nodes=None,
           cmap_nodes=None, weights_nodes=None, **kwargs):
    """
    Plot the dendrogram

    Parameters
    ----------
    children : list of pairs. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.


    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

    heights : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    ax : a pylab.axes instance (defaut is None).
        If None, create a new instance.

    active_nodes : list of int (defaut is None).
                   If active_nodes is not None, use it to add color to the
                   tree. List of nodes use for labeling.

    nodes_cmap : pylab color map (defaut is None and thus setted to pl.cm.jet)
                 Color maps used for plotting the active nodes.

    weights_nodes : list of float (defaut is None)
                   If active_nodes and weights_nodes are not None, use
                    weights_nodes to plot the branches corresponding to the
                    actives_nodes.

    Return
    ------
    a pylab.scatter instance
    """
    lx_, ly_, colors_ = _mk_dendogram(children, parent, heights, cmap_nodes,
                        active_nodes, weights_nodes)
    return _plot_graph(lx_, ly_, colors_, ax=ax, **kwargs)


def _mk_dendogram(children, parent, heights, cmap_node,
                       active_nodes=None, weights_nodes=None):
    """
    Function for computing the dendrogram

    Parameters
    ----------
    children : list of pairs. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

    heights : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    active_nodes : list of int (defaut is None).
                   If active_nodes is not None, use it to add color to the
                   tree. List of nodes use for labeling.

    nodes_cmap : pylab color map (defaut is None and thus setted to pl.cm.jet)
                 Color maps used for plotting the active nodes.

    weights_nodes : list of float (defaut is None)
                   If active_nodes and weights_nodes are not None, use
                    weights_nodes to plot the branches corresponding to the
                    actives_nodes.

    Return
    ------
    A tuple (lx_, ly_, color_) that gives the x coordinates (lx_), y
    coordinates (ly_) and color (color_) for the different nodes of the
    dendrogram.
    """

    # Compute the coordinates of the tree
    dax = np.zeros(len(children))
    dax[-1] = 0.5 * n_leaves
    x_, y_, lx_, ly_, nodes_ = [], [], [], [], []
    x_.append(dax[-1])
    y_.append(heights[-1])
    nodes_.append(len(heights))
    for i in range(len(heights), n_leaves, -1):
        ci = children[i - 1]
        dx = 0.5 * dax[-1]
        if (i - 1) != parent[i - 1]:
            dx = 0.5 * (dax[i - 1] - dax[parent[i - 1]])
        dax[ci[0]] = dax[i - 1] + dx
        dax[ci[1]] = dax[i - 1] - dx
        for j in [0, 1]:
            x_.append(dax[ci[j]])
            y_.append(heights[ci[j]])
            nodes_.append(ci[j])
            lx_.append([dax[i - 1], heights[i - 1]])
            ly_.append([dax[ci[j]], heights[ci[j]]])

    # Compute the colors of the tree
    if nodes_cmap is None:
        import pylab as pl
        nodes_cmap = pl.cm.jet
    if active_nodes is not None:
        used_scores_nodes = np.zeros(len(active_nodes), dtype=float)
        if weights_nodes is None:
            for i in range(len(active_nodes)):
                used_scores_nodes[i] = (1. + np.float(i))\
                                    / (len(active_nodes) + 1.)
        if weights_nodes is not None:
            for i in range(len(active_nodes)):
                used_scores_nodes = weights_nodes - np.min(weights_nodes)
                used_scores_nodes /= np.max(used_scores_nodes)
    colorx = np.zeros(len(children), dtype=float)
    color_ = []
    color_.append(cmap_nodes(0))
    for i in range(len(heights), n_leaves, -1):
        ci = children[i - 1]
        if active_nodes is not None:
            for j in [0, 1]:
                colorx[ci[j]] = colorx[i - 1]
                if ci[j] in active_nodes:
                    colorx[ci[j]] = \
                                used_scores_nodes[active_nodes.index(ci[j])]
                color_.append(cmap_nodes(colorx[ci[j]]))
        else:
            color_.append(cmap_nodes(0.5))
            color_.append(cmap_nodes(0.5))

    return (lx_, ly_, color_)


def _plot_graph(lx_, ly_, colors_, ax=None, **kwargs):
    """
    Function that plots a dendrogram.

    Parameters
    ----------
    lx_ : lists
    Defines the x coordinates for the different nodes of the dendrogram.

    ly_ : lists
    Defines the y coordinates for the different nodes of the dendrogram.

    color_ : lists
    Defines the color for the different nodes of the dendrogram.

    ax : a pylab.axes instance (defaut is None).
    If None, create a new instance.

    Return
    ------
    a pylab.scatter instance
    """

    import pylab as pl
    if ax is None:
        ax = pl.gca()
    line_segments = pl.matplotlib.collections.LineCollection(zip(lx_, ly_),
                                          color=color_[1:], **kwargs)
    ax.add_collection(line_segments)
    scatter = ax.scatter(x_, y_, c=color_[1:], **kwargs)
    return scatter


###############################################################################
# Class for Ward hierarchical clustering

class Ward(BaseEstimator):
    """
    Class for Ward hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    k : int or ndarray
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

    children_ : list of pairs. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    heights_ : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    labels_ : array [n_points]
        cluster labels for each point

    Return
    ------
    self
    """

    def __init__(self, k, memory=Memory(cachedir=None, verbose=0),
                 adjacency_matrix=None, copy=True, n_comp=None):
        """
        adjacency_matrix : sparse matrix.
            adjacency matrix. Defines for each sample the neigbhoring
            samples following a given structure of the data.
            Defaut is None, i.e, the hiearchical clustering algorithm is
            unstructured.
        """
        self.k = k
        self.memory = memory
        self.copy = copy
        self.n_comp = n_comp
        self.adjacency_matrix = adjacency_matrix

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
        self.parent_, self.children_, self.heights_, self.n_comp = \
                    memory.cache(ward_tree)(X, self.adjacency_matrix,
                                        n_comp=self.n_comp, copy=self.copy)

        # Cut the tree
        # self.labels_, self.active_nodes_ = memory.cache(_hc_cut)(self.k,
        self.labels_, self.active_nodes_ = _hc_cut(self.k,
                                self.parent_, self.children_, self.heights_)
        return self

###############################################################################
# Ward-based feature agglomeration

class WardAgglomeration(AgglomerationTransform, Ward):
    """Feature agglomeration based on Ward hierarchical clustering

    XXX
    """

    def fit(self, X, y=None, **params):
        return Ward.fit(self, X.T, k=self.k, **params)

