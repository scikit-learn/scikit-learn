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

from .linkage import Linkage, WardLinkage, CompleteLinkage
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

    propagate_heights : bool
        For certain linkages, it is not guaranteed that a tree node gets
        assigned a larger height than its child nodes. This may be critical
        in combination with cutting the dendrogram using cut_height. If
        propagate_heights is set to True, the height of a tree node is set
        to the maximum of the height provided by the linkage and the maximal
        height of any of its childs. Defaults to False.

    Methods
    ----------
    merge:
        Adds a new tree node corresponding to merging two clusters.

    get_descendent:
        Return all the descendent leaves of a set of nodes.

    cut:
       Cut the dendrogram for a given number of clusters.

    cut_height:
        Cut the dendrogram for a given maximal height.


    Attributes
    ----------
    n_samples : int
        The number of samples in the dataset to be clustered (i.e. the number
        of leaves of the dendrogram tree)

    n_components : int
        The number of connected components in the dataset

    parent : array, shape = [n_nodes]
        For tree node with index i, this array contains as i-th entry the index
        of its parent. Nodes without parents contain themselves as parents.

    heights : array, shape = [n_nodes]
        The heights associated to the tree's nodes.

    children : array, shape = [n_samples - 1, 2]
        Array containing the children of a node. This is not defined for leaves.

    next_index : int
        The index that a new tree node (resulting from merging) will obtain.
    """

    def __init__(self, n_samples, n_components, propagate_heights=False):
        self.n_samples = n_samples
        self.n_components = n_components
        self.propagate_heights = propagate_heights

        n_nodes = 2 * self.n_samples - 1
        self.parent = np.arange(n_nodes, dtype=np.int)
        self.heights = -np.inf * np.ones(n_nodes)
        self.children = np.zeros((self.n_samples - 1, 2), dtype=np.int) - 1

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
        self.parent[child_node_1] = self.parent[child_node_2] = self.next_index
        if self.propagate_heights:
            self.heights[self.next_index] = \
                        max(merge_distance, self.heights[child_node_1],
                            self.heights[child_node_2])
        else:
            self.heights[self.next_index] = merge_distance
        self.children[self.next_index - self.n_samples, 0] = child_node_1
        self.children[self.next_index - self.n_samples, 1] = child_node_2
        self.next_index += 1

        return self.next_index

    def get_descendent(self, ind, add_intermediate_nodes=False):
        """ Return all the descendent leaves of a set of nodes.

        Parameters
        ----------
        ind : list of int
            A list that indicates the nodes for which we want the descendents.

        add_intermediate_nodes : bool
            If true, leaves and inner nodes in the subtree are returned,
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
                ind.append(self.children[i - self.n_samples, 0])
                ind.append(self.children[i - self.n_samples, 1])
        return descendent

    def cut(self, n_clusters):
        """ Cut the dendrogram for a given number of clusters.

        Parameters
        ----------
        n_clusters : int or ndarray
            The number of clusters to form.

        Return
        ------
        cluster_roots : list of int
            A list of tree nodes. Each tree node has a height less than
            max_height and is maximal with regard to that criterion (i.e.
            its parent node in the tree has a height larger than max_height).
            The induced subtrees of these nodes are disjoint and cover the
            whole dataset.

        """
        cluster_roots = [np.max(self.children[-1]) + 1]
        for i in xrange(n_clusters - 1):
            cluster_roots.extend(self.children[np.max(cluster_roots)
                                                        - self.n_samples])
            cluster_roots.remove(np.max(cluster_roots))

        return cluster_roots

    def cut_height(self, max_height, prune=False):
        """ Cut the dendrogram for a given maximal height.

        Parameters
        ----------
        max_height : float
            The maximal height a tree node is allowed to have to form a single
            cluster. Determines indirectly how many clusters are formed.

        prune : bool
            If True, all tree nodes with an height above max_height are
            actually removed. Otherwise, the tree itself is unchanged.
            Defaults to False.

        Return
        ------
        cluster_roots : list of int
            A list of tree nodes. Each tree node has a height less than
            max_height and is maximal with regard to that criterion (i.e.
            its parent node in the tree has a height larger than max_height).
            The induced subtrees of these nodes are disjoint and cover the
            whole dataset.
        """
        open_nodes = [len(self.heights) - 1]  # root of the tree
        cluster_roots = []
        while open_nodes != []:
            node = open_nodes[0]
            open_nodes = open_nodes[1:]
            if self.heights[node] <= max_height or node <= self.n_samples:
                # This tree node is the root of a cluster
                cluster_roots.append(node)
            else:
                # Tree node induces subtree with too large height; split it
                child_node1 = self.children[node - self.n_samples, 0]
                child_node2 = self.children[node - self.n_samples, 1]
                open_nodes.append(child_node1)
                open_nodes.append(child_node2)
                open_nodes.sort(reverse=True)
                # If pruning is enabled, we also remove the parts of the tree
                # that are "too high"
                if prune:
                    self.parent[child_node1] = child_node1
                    self.parent[child_node2] = child_node2
                    self.children[node - self.n_samples] = -1

        return cluster_roots

    def get_labeling(self, cluster_roots):
        """Return labeling based on the induced subtrees of cluster root nodes

        Parameters
        ----------
        cluster_roots : list of int
            A list of tree nodes. All datapoints that are leaves of the
            induced subtree of one of these nodes gets the same label while
            leaves in different subtree get different labels. For getting a
            valid subtree, the induced subtrees of the all passed nodes
            must be disjoint and cover all datapoints.

        Return
        ------
        labels : array
            cluster labels for each point
        """

        labels = np.zeros(self.n_samples, dtype=np.int)
        for i, node in enumerate(cluster_roots):
            labels[self.get_descendent([node])] = i
        return labels


###############################################################################


def create_dendrogram(X, connectivity, n_components=None,
                      linkage_criterion="ward", metric="euclidean",
                      linkage_kwargs={}, propagate_heights=False,
                      constraints=[], copy=True):
    """Hierarchical clustering algorithm that creates a dendrogram.

    This is the structured version, that takes into account the topological
    structure between samples. The linkage is Ward's linkage per default;
    however one can use any linkage object that implements the Linkage
    interface (see linkage.py)

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

    linkage_criterion : str or subclass of Linkage
        The linkage criterion used to determine the distances of two clusters.
        This can be either "ward" or "complete" or an any subclass of the
        Linkage class.

    linkage_kwargs : dict
        Additional keyword arguments that are directly passed to the __init__
        method of the linkage object.

    propagate_heights : bool
        For certain linkages, it is not guaranteed that a tree node gets
        assigned a larger height than its child nodes. This could be critical
        in combination with cutting the dendrogram using cut_height. If
        propagate_heights is set to True, the height of a tree node is set
        to the maximum of the height provided by the linkage and the maximal
        height of any of its childs. Defaults to False.

    constraints : list of functions
        A list of constraints that are checked when merging two clusters.
        For a merge candidate consisting of two clusters with datapoints with
        indices c1 and c2, each constraint is checked. If constraint(c1, c2)
        is not True, the two clusters are not merged.
        Defaults to empty list.

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

    # Determine linkage class for known linkage names
    if linkage_criterion == "ward":
        linkage_criterion = WardLinkage
    elif linkage_criterion == "complete":
        linkage_criterion = CompleteLinkage

    # Linkage criterion must be a subclass of Linkage.
    assert (issubclass(linkage_criterion, Linkage)), \
        "linkage_criterion must be a subclass of Linkage or a known " \
        "linkage name."

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
    dendrogram = Dendrogram(n_samples, n_components, propagate_heights)

    # Linkage object that manages the distance computations between clusters
    # distances are updated incrementally during merging of clusters...
    linkage = linkage_criterion(X, connectivity, metric=metric,
                                **linkage_kwargs)

    # Array in which open_nodes[i] indicate whether the cluster with index i
    # has not yet been merged into a larger cluster
    open_nodes = np.ones(n_nodes, dtype=bool)

    # Recursive merge loop
    k = n_samples
    while k < n_nodes:
        # Identify the merge that will be applied next
        # This is the merge with minimal distance of two cluster that fulfill 
        # all constraints and haven't been merged into a larger cluster yet. 
        merge_distance = np.inf
        while linkage.has_more_candidates():
            merge_distance, i, j = linkage.fetch_candidate()
            # Check if all constraints are fulfilled
            constraints_satisfied = True
            for constraint in constraints:
                if not constraint(linkage.get_nodes_of_cluster(i),
                                  linkage.get_nodes_of_cluster(j)):
                    constraints_satisfied = False
                    break
            if not constraints_satisfied:
                # i and j have merge distance infinity because they violate 
                # at least one of the constraints
                merge_distance = np.inf
                continue
            # Check if nodes haven't been merged already
            if open_nodes[i] and open_nodes[j]:
                break  # i and j can be merged!
            # Check if only unconnected components are left
            if not linkage.has_more_candidates():
                merge_distance = np.inf
                break

        # Check if we have fully merged all connected components
        if merge_distance == np.inf and not linkage.has_more_candidates():
            # Merge unconnected components with height infinity
            i = k - 1
            desc = set(dendrogram.get_descendent([i],
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
                                       linkage_criterion="ward",
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
        Default is None, i.e, the hierarchical clustering algorithm is
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

    linkage_criterion : str or subclass of Linkage
        The linkage criterion used to determine the distances of two clusters.
        In the structured case (connectivity not None), this can be either
        "ward" or "complete" or an any subclass of the Linkage class.
        In the unstructured case (connectivity is None), this can be any of
        "average", "centroid", "complete", "median", "single", "ward",
        or "weighted". Defaults to "ward".

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

    # Mapping from name to scipy implementation of hierarchical cluster
    # algorithm
    unstructured_cluster_algorithms = {"average": hierarchy.average,
                                       "centroid": hierarchy.centroid,
                                       "complete": hierarchy.complete,
                                       "median": hierarchy.median,
                                       "single": hierarchy.single,
                                       "ward": hierarchy.ward,
                                       "weighted": hierarchy.weighted}

    def __init__(self, n_clusters=None, max_height=None, connectivity=None,
                 n_components=None, linkage_criterion="ward",
                 metric="euclidean", linkage_kwargs={}, 
                 memory=Memory(cachedir=None, verbose=0), copy=True):
        self.n_clusters = n_clusters
        self.max_height = max_height
        self.connectivity = connectivity
        self.n_components = n_components
        self.linkage_criterion = linkage_criterion
        self.metric = metric
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
            self.dendrogram_ = memory.cache(create_dendrogram)(X,
                                    self.connectivity,
                                    n_components=self.n_components,
                                    linkage_criterion=self.linkage_criterion,
                                    metric=self.metric,
                                    linkage_kwargs=self.linkage_kwargs,
                                    copy=self.copy)
            # Cut the tree ...
            if self.max_height is None:
                # based on number of desired clusters
                cluster_roots_ = self.dendrogram_.cut(self.n_clusters)
            else:
                # based on maximally allowed height
                cluster_roots_ = self.dendrogram_.cut_height(self.max_height)
        else:  # Fall back to scipy
            assert self.n_clusters is not None, \
                "Unstructured clustering requires the number of clusters to "\
                "be specified explicitly."

            assert self.linkage_criterion \
                        in self.unstructured_cluster_algorithms, \
                    "Unknown linkage criterion %s. Must be one of %s." \
                        % (self.linkage_criterion,
                           self.unstructured_cluster_algorithms.keys())
            # Invoke clustering algorithm
            clustering_algorithm = \
                self.unstructured_cluster_algorithms[self.linkage_criterion]
            out = clustering_algorithm(X)
            # Put result into a dendrogram and cut it to get labeling
            self.dendrogram_ = Dendrogram(X.shape[0], 1)
            self.dendrogram_.children = out[:, :2].astype(np.int)
            cluster_roots_ = self.dendrogram_.cut(self.n_clusters)

        # Determine labeling of datapoints based on the induced subtrees
        # of the cut of the dendrogram
        self.labels_ = self.dendrogram_.get_labeling(cluster_roots_)

        return self


###############################################################################

class Ward(HierarchicalClustering):
    # TODO: Contained only for backward compatibility. Add a deprecation
    #       warning?
    pass


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
