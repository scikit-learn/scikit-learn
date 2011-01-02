"""
These routines perform some hierachical agglomerative clustering of some input
data. Currently, only Ward's algorithm is implemented.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort
License: BSD 3 clause
"""
import numpy as np
import heapq as heapq
from scipy import sparse
import _inertia
from scikits.learn.base import BaseEstimator
from scikits.learn.utils._csgraph import cs_graph_components

###########################################################################
# Ward's algorithm

def ward_tree(X, adjacency_matrix=None):
    """Ward clustering based on a Feature matrix. Heapq-based representation
    of the inertia matrix.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X:  array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    adjacency_matrix : sparse matrix.
        adjacency matrix. Defines for each sample the neigbhoring samples
        following a given structure of the data.
        Defaut is None, i.e, the ward algorithm is unstructured.

    Returns
    -------
    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

    children : list of two-elements lists. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    height : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    adjacency_matrix : sparse matrix.
        The update version of adjacency matrix. Defines for each node the
        neigbhoring nodes following a given structure of the data.

    """
    X = np.asanyarray(X)
    n_samples, n_features = X.shape
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))

    # Adjacency matrix
    if adjacency_matrix is None:
        adjacency_matrix = np.ones([X.shape[0], X.shape[0]])
        adjacency_matrix[range(X.shape[0]), range(X.shape[0])] = 0
        adjacency_matrix = sparse.lil_matrix(adjacency_matrix)
        n_nodes = 2 * n_samples - 1
    else:
        adjacency_matrix = adjacency_matrix.tolil()

    # Remove diagonal from adjacency matrix
    adjacency_matrix.setdiag(np.zeros(adjacency_matrix.shape[0]))

    # Compute the number of nodes
    n_comp, label = cs_graph_components(adjacency_matrix)
    n_nodes = 2 * n_samples - n_comp
    if n_comp > 1:
        print "Warning: the number of connected compoments of the" + \
    " adjacency matrix is ", n_comp, " > 1. The tree will be stopped early."

    # build moments as a list
    moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features)),
                np.zeros((n_nodes, n_features))]
    moments[0][:n_samples] = 1
    moments[1][:n_samples] = X
    moments[2][:n_samples] = X ** 2

    # create a inertia matrix
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
    parent = np.arange(n_nodes).astype(np.int)
    height = np.zeros(n_nodes)
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
        parent[i], parent[j], height[k] = k, k, node[0]
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

    return parent, children, height, A


###########################################################################
# Functions for cutting  hierarchical clustering tree

def _hc_get_descendent(ind, children):
    """
    Function returning all the descendent leaves of a set of nodes in the tree.

    Parameters
    ----------
    ind : list of int
          A list that indicates the nodes for which we want the descendents.

    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

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


def _hc_cut(k, parent, children, height):
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

    height : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    Return
    ------
    label_ : array [n_points]
        cluster labels for each point

    active_nodes : list of int
                index of the nodes kept for the labeling
    """
    parent = parent[:-1]
    height = height[:-1]
    active_nodes = [len(parent)]
    node_to_cut = active_nodes[0]
    for i in range(k - 1):
        if np.sum(parent == node_to_cut) != 0:
            active_nodes.append(np.where(parent == node_to_cut)[0][0])
            active_nodes.append(np.where(parent == node_to_cut)[0][1])
            active_nodes.remove(node_to_cut)
        else:
            active_nodes.append(node_to_cut)
        node_to_cut = active_nodes[np.argmax(height[active_nodes])]
    label = np.zeros(children.count([]))
    for j in active_nodes[:k - 1]:
        ind = [j]
        label[_hc_get_descendent(ind, children)] = np.max(label) + 1
    return label, active_nodes


###########################################################################
# Display functions for hierarchical clustering

def plot_dendrogram(axe, parent, children, height, active_nodes=None,
           weights_nodes=None, cmap_nodes=None, **kwargs):
    """
    Plot the dendrogram

    Parameters
    ----------
    axe : a pylab.axes instance

    parent : array-like, shape = [n_nodes]
            Int. Gives the parent node for each node, i.e. parent[i] is the
            parent node of the node i. The last value of parent is the
            root node, that is its self parent, so the last value is taken
            3 times in the array.
            The n_nodes is equal at  (2*n_samples - 1), and takes into
            account the nb_samples leaves, and the unique root.

    children : list of two-elements lists. Lenght of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    height : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    active_nodes : list of int (defaut is None).
                   If active_nodes is not None, use it to add color to the
                   tree. List of nodes use for labeling.

    weights_nodes : list of float (defaut is None)
                   If active_nodes and weights_nodes are not None, use
                    weights_nodes to plot the branches corresponding to the
                    actives_nodes.

    cmap_nodes : pylab color map (defaut is None and thus setted to pl.cm.jet)
                 Color maps used for plotting the active nodes.

    Return
    ------
    axe : a pylab axe instance
    """
    import pylab as pl
    n_leaves = np.min(parent)

    # Compute the coordinates of the tree
    dax = np.zeros(len(children))
    dax[-1] = 0.5 * n_leaves
    x_, y_, lx_, ly_, nodes_ = [], [], [], [], []
    x_.append(dax[-1])
    y_.append(height[-1])
    nodes_.append(len(height))
    for i in range(len(height), n_leaves, -1):
        ci = children[i - 1]
        dx = 0.5 * dax[-1]
        if (i - 1) != parent[i - 1]:
            dx = 0.5 * (dax[i - 1] - dax[parent[i - 1]])
        dax[ci[0]] = dax[i - 1] + dx
        dax[ci[1]] = dax[i - 1] - dx
        for j in [0, 1]:
            x_.append(dax[ci[j]])
            y_.append(height[ci[j]])
            nodes_.append(ci[j])
            lx_.append([dax[i - 1], height[i - 1]])
            ly_.append([dax[ci[j]], height[ci[j]]])

    # Compute the colors of the tree
    if cmap_nodes is None:
            cmap_nodes = pl.cm.jet
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
    for i in range(len(height), n_leaves, -1):
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

    # Plot the tree in axe
    color_ = np.array(color_)
    line_segments = pl.matplotlib.collections.LineCollection(zip(lx_, ly_),
                                            color=color_[1:], **kwargs)
    axe.add_collection(line_segments)
    axe.scatter(x_, y_, c=color_[1:], **kwargs)
    return axe


##############################################################################
# Classes for hiearchical clustering

class HierarchicalClustering(BaseEstimator):
    """
    General class for hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray
                 The number of clusters.

    tree_func: function that takes a matrix X of data, and possibly an
        adjacency matrix.
        Defaut is ward_tree (construction of the tree by ward algorithm).

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

    height_ : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    label_ : array [n_points]
        cluster labels for each point

    """

    def __init__(self, n_clusters, tree_func=ward_tree):
        self.tree_func = tree_func
        self.n_clusters = n_clusters

    def fit(self, X, adjacency_matrix=None, copy=True, **params):
        """
        Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            A M by N array of M observations in N dimensions or a length
            M array of M one-dimensional observations.

        adjacency_matrix : sparse matrix.
            adjacency matrix. Defines for each sample the neigbhoring
            samples following a given structure of the data.
            Defaut is None, i.e, the hiearchical clustering algorithm is
            unstructured.

        Returns
        -------
        self
        """
        self._set_params(**params)

        # If necessary, copy the adjacency matrix
        if copy and adjacency_matrix is not None:
            self.adjacency_matrix = adjacency_matrix.copy()
        else:
            self.adjacency_matrix = adjacency_matrix

        # Check if the adjacency matrix is well-connected
        if self.adjacency_matrix is not None:
            self.adjacency_matrix = self.adjacency_matrix.tolil()
            self.adjacency_matrix.setdiag(
                                    np.zeros(self.adjacency_matrix.shape[0]))
            n_comp, label = cs_graph_components(self.adjacency_matrix)
            if n_comp > 1:
                print "Warning: the number of connected compoments of the" + \
                " adjacency matrix is > 1. The tree will be stopped early," + \
                " and the maximal number of clusters will be ", n_comp
            self.n_clusters = np.max([self.n_clusters, n_comp])

        # Construct the tree
        self.parent_, self.children_, self.height_, self.adjacency_matrix = \
                                    self.tree_func(X, self.adjacency_matrix)

        # Cut the tree
        self.label_, self.active_nodes_ = _hc_cut(self.n_clusters,
                                self.parent_, self.children_, self.height_)
        return self


class Ward(HierarchicalClustering):
    """
    Class for ward clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray
                 The number of clusters.

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

    height_ : array-like, shape = [n_nodes]
            Float. Gives the inertia of the created nodes. The n_samples first
            values of the array are 0, and thus the values are positive (or
            null) and are ranked in an increasing order.

    label_ : array [n_points]
        cluster labels for each point

    """

    def __init__(self, n_clusters):
        self.tree_func = ward_tree
        self.n_clusters = n_clusters
