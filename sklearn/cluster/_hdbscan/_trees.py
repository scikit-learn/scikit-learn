# -*- coding: utf-8 -*-
# Author: Leland McInnes <leland.mcinnes@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from ._hdbscan_tree import compute_stability, labelling_at_cut, recurse_leaf_dfs

CB_LEFT = 0
CB_RIGHT = 1
CB_BOTTOM = 2
CB_TOP = 3


def _bfs_from_cluster_tree(tree, bfs_root):
    """
    Perform a breadth first search on a tree in condensed tree format
    """

    result = []
    to_process = [bfs_root]

    while to_process:
        result.extend(to_process)
        to_process = tree["child"][np.in1d(tree["parent"], to_process)].tolist()

    return result


def _recurse_leaf_dfs(cluster_tree, current_node):
    children = cluster_tree[cluster_tree["parent"] == current_node]["child"]
    if len(children) == 0:
        return [
            current_node,
        ]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


def _get_leaves(condensed_tree):
    cluster_tree = condensed_tree[condensed_tree["child_size"] > 1]
    if cluster_tree.shape[0] == 0:
        # Return the only cluster, the root
        return [condensed_tree["parent"].min()]

    root = cluster_tree["parent"].min()
    return _recurse_leaf_dfs(cluster_tree, root)


class CondensedTree(object):
    """The condensed tree structure, which provides a simplified or smoothed version
    of the :class:`~hdbscan.plots.SingleLinkageTree`.

    Parameters
    ----------
    condensed_tree_array : numpy recarray from :class:`~hdbscan.HDBSCAN`
        The raw numpy rec array version of the condensed tree as produced
        internally by hdbscan.

    cluster_selection_method : string, optional (default 'eom')
        The method of selecting clusters. One of 'eom' or 'leaf'

    allow_single_cluster : Boolean, optional (default False)
        Whether to allow the root cluster as the only selected cluster

    """

    def __init__(
        self,
        condensed_tree_array,
        cluster_selection_method="eom",
        allow_single_cluster=False,
    ):
        self._raw_tree = condensed_tree_array
        self.cluster_selection_method = cluster_selection_method
        self.allow_single_cluster = allow_single_cluster

    def _select_clusters(self):
        if self.cluster_selection_method == "eom":
            stability = compute_stability(self._raw_tree)
            if self.allow_single_cluster:
                node_list = sorted(stability.keys(), reverse=True)
            else:
                node_list = sorted(stability.keys(), reverse=True)[:-1]
            cluster_tree = self._raw_tree[self._raw_tree["child_size"] > 1]
            is_cluster = {cluster: True for cluster in node_list}

            for node in node_list:
                child_selection = cluster_tree["parent"] == node
                subtree_stability = np.sum(
                    [
                        stability[child]
                        for child in cluster_tree["child"][child_selection]
                    ]
                )

                if subtree_stability > stability[node]:
                    is_cluster[node] = False
                    stability[node] = subtree_stability
                else:
                    for sub_node in _bfs_from_cluster_tree(cluster_tree, node):
                        if sub_node != node:
                            is_cluster[sub_node] = False

            return sorted([cluster for cluster in is_cluster if is_cluster[cluster]])

        elif self.cluster_selection_method == "leaf":
            return _get_leaves(self._raw_tree)
        else:
            raise ValueError(
                "Invalid Cluster Selection Method: %s\n"
                'Should be one of: "eom", "leaf"\n'
            )

    def to_numpy(self):
        """Return a numpy structured array representation of the condensed tree."""
        return self._raw_tree.copy()


def _get_dendrogram_ordering(parent, linkage, root):

    if parent < root:
        return []

    return (
        _get_dendrogram_ordering(int(linkage[parent - root][0]), linkage, root)
        + _get_dendrogram_ordering(int(linkage[parent - root][1]), linkage, root)
        + [parent]
    )


class SingleLinkageTree(object):
    """A single linkage format dendrogram tree, with plotting functionality
    and networkX support.

    Parameters
    ----------
    linkage : ndarray (n_samples, 4)
        The numpy array that holds the tree structure. As output by
        scipy.cluster.hierarchy, hdbscan, of fastcluster.

    """

    def __init__(self, linkage):
        self._linkage = linkage

    def to_numpy(self):
        """Return a numpy array representation of the single linkage tree.

        This representation conforms to the scipy.cluster.hierarchy notion
        of a single linkage tree, and can be used with all the associated
        scipy tools. Please see the scipy documentation for more details
        on the format.
        """
        return self._linkage.copy()

    def get_clusters(self, cut_distance, min_cluster_size=5):
        """Return a flat clustering from the single linkage hierarchy.

        This represents the result of selecting a cut value for robust single linkage
        clustering. The `min_cluster_size` allows the flat clustering to declare noise
        points (and cluster smaller than `min_cluster_size`).

        Parameters
        ----------

        cut_distance : float
            The mutual reachability distance cut value to use to generate a
            flat clustering.

        min_cluster_size : int, default=5
            Clusters smaller than this value with be called 'noise' and remain
            unclustered in the resulting flat clustering.

        Returns
        -------

        labels : array (n_samples,)
            An array of cluster labels, one per datapoint. Unclustered points
            are assigned the label -1.
        """
        return labelling_at_cut(self._linkage, cut_distance, min_cluster_size)


class MinimumSpanningTree(object):
    def __init__(self, mst, data):
        self._mst = mst
        self._data = data

    def to_numpy(self):
        """Return a numpy array of weighted edges in the minimum spanning tree"""
        return self._mst.copy()
