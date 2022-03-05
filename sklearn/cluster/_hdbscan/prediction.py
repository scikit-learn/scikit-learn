# Support various prediction methods for predicting cluster membership
# of new or unseen points. There are several ways to interpret how
# to do this correctly, so we provide several methods for
# the different use cases that may arise.

import numpy as np

from sklearn.neighbors import KDTree, BallTree
from .dist_metrics import DistanceMetric
from ._hdbscan_tree import recurse_leaf_dfs
from ._prediction_utils import (
    get_tree_row_with_child,
    dist_membership_vector,
    outlier_membership_vector,
    prob_in_some_cluster,
    all_points_dist_membership_vector,
    all_points_outlier_membership_vector,
    all_points_prob_in_some_cluster,
)
from warnings import warn


class PredictionData(object):
    """
    Extra data that allows for faster prediction if cached.

    Parameters
    ----------

    data : array (n_samples, n_features)
        The original data set that was clustered.

    condensed_tree : CondensedTree
        The condensed tree object created by a clustering.

    min_samples : int
        The min_samples value used in clustering.

    tree_type : str, default="kdtree"
        Which type of space tree to use for core distance computation.
        One of:
            * ``kdtree``
            * ``balltree``

    metric : str, default="euclidean"
        The metric used to determine distance for the clustering.
        This is the metric that will be used for the space tree to determine
        core distances etc.

    **kwargs :
        Any further arguments to the metric.

    Attributes
    ----------

    raw_data : array (n_samples, n_features)
        The original data set that was clustered

    tree : KDTree or BallTree
        A space partitioning tree that can be queried for nearest neighbors.

    core_distances : array (n_samples,)
        The core distances for every point in the original data set.

    cluster_map : dict
        A dictionary mapping cluster numbers in the condensed tree to labels
        in the final selected clustering.

    cluster_tree : structured array
        A version of the condensed tree that only contains clusters, not
        individual points.

    max_lambdas : dict
        A dictionary mapping cluster numbers in the condensed tree to the
        maximum lambda value seen in that cluster.
    """

    _tree_type_map = {"kdtree": KDTree, "balltree": BallTree}

    def _clusters_below(self, cluster):
        result = []
        to_process = [cluster]

        while to_process:
            result.extend(to_process)
            to_process = self.cluster_tree["child"][
                np.in1d(self.cluster_tree["parent"], to_process)
            ]
            to_process = to_process.tolist()

        return result

    def _recurse_leaf_dfs(self, current_node):
        children = self.cluster_tree[self.cluster_tree["parent"] == current_node][
            "child"
        ]
        if len(children) == 0:
            return [
                current_node,
            ]
        else:
            return sum(
                [recurse_leaf_dfs(self.cluster_tree, child) for child in children], []
            )

    def __init__(
        self,
        data,
        condensed_tree,
        min_samples,
        tree_type="kdtree",
        metric="euclidean",
        **kwargs,
    ):
        self.raw_data = data.astype(np.float64)
        self.tree = self._tree_type_map[tree_type](
            self.raw_data, metric=metric, **kwargs
        )
        self.core_distances = self.tree.query(data, k=min_samples)[0][:, -1]
        self.dist_metric = DistanceMetric.get_metric(metric, **kwargs)

        selected_clusters = sorted(condensed_tree._select_clusters())
        # raw_condensed_tree = condensed_tree.to_numpy()
        raw_condensed_tree = condensed_tree._raw_tree

        self.cluster_map = {c: n for n, c in enumerate(sorted(list(selected_clusters)))}
        self.reverse_cluster_map = {n: c for c, n in self.cluster_map.items()}

        self.cluster_tree = raw_condensed_tree[raw_condensed_tree["child_size"] > 1]
        self.max_lambdas = {}
        self.leaf_max_lambdas = {}
        self.exemplars = []

        all_clusters = set(
            np.hstack([self.cluster_tree["parent"], self.cluster_tree["child"]])
        )

        for cluster in all_clusters:
            self.leaf_max_lambdas[cluster] = raw_condensed_tree["lambda_val"][
                raw_condensed_tree["parent"] == cluster
            ].max()

        for cluster in selected_clusters:
            self.max_lambdas[cluster] = raw_condensed_tree["lambda_val"][
                raw_condensed_tree["parent"] == cluster
            ].max()

            for sub_cluster in self._clusters_below(cluster):
                self.cluster_map[sub_cluster] = self.cluster_map[cluster]
                self.max_lambdas[sub_cluster] = self.max_lambdas[cluster]

            cluster_exemplars = np.array([], dtype=np.int64)
            for leaf in self._recurse_leaf_dfs(cluster):
                leaf_max_lambda = raw_condensed_tree["lambda_val"][
                    raw_condensed_tree["parent"] == leaf
                ].max()
                points = raw_condensed_tree["child"][
                    (raw_condensed_tree["parent"] == leaf)
                    & (raw_condensed_tree["lambda_val"] == leaf_max_lambda)
                ]
                cluster_exemplars = np.hstack([cluster_exemplars, points])

            self.exemplars.append(self.raw_data[cluster_exemplars])


def _find_neighbor_and_lambda(
    neighbor_indices, neighbor_distances, core_distances, min_samples
):
    """
    Find the nearest mutual reachability neighbor of a point, and  compute
    the associated lambda value for the point, given the mutual reachability
    distance to a nearest neighbor.

    Parameters
    ----------
    neighbor_indices : array (2 * min_samples, )
        An array of raw distance based nearest neighbor indices.

    neighbor_distances : array (2 * min_samples, )
        An array of raw distances to the nearest neighbors.

    core_distances : array (n_samples, )
        An array of core distances for all points

    min_samples : int
        The min_samples value used to generate core distances.

    Returns
    -------
    neighbor : int
        The index into the full raw data set of the nearest mutual reachability
        distance neighbor of the point.

    lambda_ : float
        The lambda value at which this point joins/merges with `neighbor`.
    """
    neighbor_core_distances = core_distances[neighbor_indices]
    point_core_distances = neighbor_distances[min_samples] * np.ones(
        neighbor_indices.shape[0]
    )
    mr_distances = np.vstack(
        (neighbor_core_distances, point_core_distances, neighbor_distances)
    ).max(axis=0)

    nn_index = mr_distances.argmin()

    nearest_neighbor = neighbor_indices[nn_index]
    if mr_distances[nn_index] > 0.0:
        lambda_ = 1.0 / mr_distances[nn_index]
    else:
        lambda_ = np.finfo(np.double).max

    return nearest_neighbor, lambda_


def _extend_condensed_tree(
    tree, neighbor_indices, neighbor_distances, core_distances, min_samples
):
    """
    Create a new condensed tree with an additional point added, allowing for
    computations as if this point had been part of the original tree. Note
    that this makes as little change to the tree as possible, with no
    re-optimizing/re-condensing so that the selected clusters remain
    effectively unchanged.

    Parameters
    ----------
    tree : structured array
        The raw format condensed tree to update.

    neighbor_indices : array (2 * min_samples, )
        An array of raw distance based nearest neighbor indices.

    neighbor_distances : array (2 * min_samples, )
        An array of raw distances to the nearest neighbors.

    core_distances : array (n_samples, )
        An array of core distances for all points

    min_samples : int
        The min_samples value used to generate core distances.

    Returns
    -------
    new_tree : structured array
        The original tree with an extra row providing the parent cluster
        and lambda information for a new point given index -1.
    """
    tree_root = tree["parent"].min()

    nearest_neighbor, lambda_ = _find_neighbor_and_lambda(
        neighbor_indices, neighbor_distances, core_distances, min_samples
    )

    neighbor_tree_row = get_tree_row_with_child(tree, nearest_neighbor)
    potential_cluster = neighbor_tree_row["parent"]

    if neighbor_tree_row["lambda_val"] <= lambda_:
        # New point departs with the old
        new_tree_row = (potential_cluster, -1, 1, neighbor_tree_row["lambda_val"])
    else:
        # Find appropriate cluster based on lambda of new point
        while (
            potential_cluster > tree_root
            and tree[tree["child"] == potential_cluster]["lambda_val"] >= lambda_
        ):
            potential_cluster = tree["parent"][tree["child"] == potential_cluster][0]

        new_tree_row = (potential_cluster, -1, 1, lambda_)

    return np.append(tree, new_tree_row)


def _find_cluster_and_probability(
    tree,
    cluster_tree,
    neighbor_indices,
    neighbor_distances,
    core_distances,
    cluster_map,
    max_lambdas,
    min_samples,
):
    """
    Return the cluster label (of the original clustering) and membership
    probability of a new data point.

    Parameters
    ----------
    tree : CondensedTree
        The condensed tree associated with the clustering.

    cluster_tree : structured_array
        The raw form of the condensed tree with only cluster information (no
        data on individual points). This is significantly more compact.

    neighbor_indices : array (2 * min_samples, )
        An array of raw distance based nearest neighbor indices.

    neighbor_distances : array (2 * min_samples, )
        An array of raw distances to the nearest neighbors.

    core_distances : array (n_samples, )
        An array of core distances for all points

    cluster_map : dict
        A dictionary mapping cluster numbers in the condensed tree to labels
        in the final selected clustering.

    max_lambdas : dict
        A dictionary mapping cluster numbers in the condensed tree to the
        maximum lambda value seen in that cluster.

    min_samples : int
        The min_samples value used to generate core distances.
    """
    raw_tree = tree._raw_tree
    tree_root = cluster_tree["parent"].min()

    nearest_neighbor, lambda_ = _find_neighbor_and_lambda(
        neighbor_indices, neighbor_distances, core_distances, min_samples
    )

    neighbor_tree_row = get_tree_row_with_child(raw_tree, nearest_neighbor)
    potential_cluster = neighbor_tree_row["parent"]

    if neighbor_tree_row["lambda_val"] > lambda_:
        # Find appropriate cluster based on lambda of new point
        while (
            potential_cluster > tree_root
            and cluster_tree["lambda_val"][cluster_tree["child"] == potential_cluster]
            >= lambda_
        ):
            potential_cluster = cluster_tree["parent"][
                cluster_tree["child"] == potential_cluster
            ][0]

    if potential_cluster in cluster_map:
        cluster_label = cluster_map[potential_cluster]
    else:
        cluster_label = -1

    if cluster_label >= 0:
        max_lambda = max_lambdas[potential_cluster]

        if max_lambda > 0.0:
            lambda_ = min(max_lambda, lambda_)
            prob = lambda_ / max_lambda
        else:
            prob = 1.0
    else:
        prob = 0.0

    return cluster_label, prob


def approximate_predict(clusterer, points_to_predict):
    """
    Predict the cluster label of new points.

    The returned labels will be those of the original clustering found by
    ``clusterer``, and therefore are not (necessarily) the cluster labels that
    would be found by clustering the original data combined with
    ``points_to_predict``, hence the 'approximate' label.

    If you simply wish to assign new points to an existing clustering
    in the 'best' way possible, this is the function to use. If you
    want to predict how ``points_to_predict`` would cluster with
    the original data under HDBSCAN the most efficient existing approach
    is to simply recluster with the new point(s) added to the original dataset.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and
        either had ``prediction_data=True`` set, or called the
        ``generate_prediction_data`` method after the fact.

    points_to_predict : array, or array-like (n_samples, n_features)
        The new data points to predict cluster labels for. They should
        have the same dimensionality as the original dataset over which
        clusterer was fit.

    Returns
    -------
    labels : array (n_samples,)
        The predicted labels of the ``points_to_predict``.

    probabilities : array (n_samples,)
        The soft cluster scores for each of the ``points_to_predict``.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.membership_vector : Predict soft cluster
        membership.
    sklearn.cluster.hdbscan.prediction.all_points_membership_vectors : Predict
        soft cluster membership vectors for all points in the original dataset
        the clusterer was trained on.
    """
    if clusterer.prediction_data_ is None:
        raise ValueError(
            "Clusterer does not have prediction data!"
            " Try fitting with prediction_data=True set,"
            " or run generate_prediction_data on the clusterer"
        )

    points_to_predict = np.asarray(points_to_predict)

    if points_to_predict.shape[1] != clusterer.prediction_data_.raw_data.shape[1]:
        raise ValueError("New points dimension does not match fit data!")

    if clusterer.prediction_data_.cluster_tree.shape[0] == 0:
        warn(
            "Clusterer does not have any defined clusters, new data"
            " will be automatically predicted as noise."
        )
        labels = -1 * np.ones(points_to_predict.shape[0], dtype=np.int32)
        probabilities = np.zeros(points_to_predict.shape[0], dtype=np.float32)
        return labels, probabilities

    labels = np.empty(points_to_predict.shape[0], dtype=np.int32)
    probabilities = np.empty(points_to_predict.shape[0], dtype=np.float64)

    min_samples = clusterer.min_samples or clusterer.min_cluster_size
    neighbor_distances, neighbor_indices = clusterer.prediction_data_.tree.query(
        points_to_predict, k=2 * min_samples
    )

    for i in range(points_to_predict.shape[0]):
        label, prob = _find_cluster_and_probability(
            clusterer.condensed_tree_,
            clusterer.prediction_data_.cluster_tree,
            neighbor_indices[i],
            neighbor_distances[i],
            clusterer.prediction_data_.core_distances,
            clusterer.prediction_data_.cluster_map,
            clusterer.prediction_data_.max_lambdas,
            min_samples,
        )
        labels[i] = label
        probabilities[i] = prob

    return labels, probabilities


def approximate_predict_scores(clusterer, points_to_predict):
    """
    Predict the outlier score of new points.

    The returned scores will be based on the original clustering found by
    ``clusterer``, and therefore are not (necessarily) the outlier scores that
    would be found by clustering the original data combined with
    ``points_to_predict``, hence the 'approximate' label.

    If you simply wish to calculate the outlier scores for new points
    in the 'best' way possible, this is the function to use. If you
    want to predict the outlier score of ``points_to_predict`` with
    the original data under HDBSCAN the most efficient existing approach
    is to simply recluster with the new point(s) added to the original dataset.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and
        either had ``prediction_data=True`` set, or called the
        ``generate_prediction_data`` method after the fact.

    points_to_predict : array, or array-like (n_samples, n_features)
        The new data points to predict cluster labels for. They should
        have the same dimensionality as the original dataset over which
        clusterer was fit.

    Returns
    -------
    scores : array (n_samples,)
        The predicted scores of the ``points_to_predict``.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.membership_vector : Predict soft cluster
        membership.
    sklearn.cluster.hdbscan.prediction.all_points_membership_vectors : Predict
        soft cluster membership vectors for all points in the original dataset
        the clusterer was trained on.
    """
    try:
        clusterer.prediction_data_
    except AttributeError:
        raise ValueError(
            "Clusterer does not have prediction data!"
            " Try fitting with prediction_data=True set,"
            " or run generate_prediction_data on the clusterer"
        )

    points_to_predict = np.asarray(points_to_predict)

    if points_to_predict.shape[1] != clusterer.prediction_data_.raw_data.shape[1]:
        raise ValueError("New points dimension does not match fit data!")

    if clusterer.prediction_data_.cluster_tree.shape[0] == 0:
        warn(
            "Clusterer does not have any defined clusters, new data"
            " will be automatically predicted as outliers."
        )
        scores = np.ones(points_to_predict.shape[0], dtype=np.int32)
        return scores

    scores = np.empty(points_to_predict.shape[0], dtype=np.float64)

    min_samples = clusterer.min_samples or clusterer.min_cluster_size
    neighbor_distances, neighbor_indices = clusterer.prediction_data_.tree.query(
        points_to_predict, k=2 * min_samples
    )

    tree = clusterer.condensed_tree_._raw_tree

    parent_array = tree["parent"]

    tree_root = parent_array.min()
    max_lambdas = {}
    for parent in np.unique(tree["parent"]):
        max_lambdas[parent] = tree[tree["parent"] == parent]["lambda_val"].max()

    for n in np.argsort(parent_array):
        cluster = tree["child"][n]
        if cluster < tree_root:
            break

        parent = parent_array[n]
        if max_lambdas[cluster] > max_lambdas[parent]:
            max_lambdas[parent] = max_lambdas[cluster]

    for i in range(points_to_predict.shape[0]):
        neigh, lambda_ = _find_neighbor_and_lambda(
            neighbor_indices[i],
            neighbor_distances[i],
            clusterer.prediction_data_.core_distances,
            min_samples,
        )

        neighbor_tree_row = get_tree_row_with_child(tree, neigh)
        potential_cluster = neighbor_tree_row["parent"]

        if neighbor_distances[i].min() == 0:
            # the point is in the dataset, fix lambda for rounding errors
            lambda_ = neighbor_tree_row["lambda_val"]

        max_lambda = max_lambdas[potential_cluster]

        if max_lambda > 0.0:
            scores[i] = (max_lambda - lambda_) / max_lambda
        else:
            scores[i] = 0.0

    return scores


def membership_vector(clusterer, points_to_predict):
    """
    Predict soft cluster membership.

    Predicts sofr cluster membership, producing a vector for each point in
    ``points_to_predict`` that gives a probability that the given point is a
    member of a cluster for each of the selected clusters of the ``clusterer``.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and
        either had ``prediction_data=True`` set, or called the
        ``generate_prediction_data`` method after the fact.

    points_to_predict : array, or array-like (n_samples, n_features)
        The new data points to predict cluster labels for. They should
        have the same dimensionality as the original dataset over which
        clusterer was fit.

    Returns
    -------
    membership_vectors : array (n_samples, n_clusters)
        The probability that point ``i`` is a member of cluster ``j`` is
        in ``membership_vectors[i, j]``.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.approximate_predict : Predict the
        cluster label of new points.
    sklearn.cluster.hdbscan.prediction.all_points_membership_vectors : Predict
        soft cluster membership vectors for all points in the original dataset
        the clusterer was trained on.
    """

    points_to_predict = points_to_predict.astype(np.float64)
    clusters = np.array(
        sorted(list(clusterer.condensed_tree_._select_clusters()))
    ).astype(np.intp)

    result = np.empty((points_to_predict.shape[0], clusters.shape[0]), dtype=np.float64)

    min_samples = clusterer.min_samples or clusterer.min_cluster_size
    neighbor_distances, neighbor_indices = clusterer.prediction_data_.tree.query(
        points_to_predict, k=2 * min_samples
    )

    for i in range(points_to_predict.shape[0]):

        # We need to find where in the tree the new point would go
        # for the purposes of outlier membership approximation
        nearest_neighbor, lambda_ = _find_neighbor_and_lambda(
            neighbor_indices[i],
            neighbor_distances[i],
            clusterer.prediction_data_.core_distances,
            min_samples,
        )

        neighbor_tree_row = get_tree_row_with_child(
            clusterer.condensed_tree_._raw_tree, nearest_neighbor
        )

        if neighbor_tree_row["lambda_val"] <= lambda_:
            lambda_ = neighbor_tree_row["lambda_val"]

        distance_vec = dist_membership_vector(
            points_to_predict[i],
            clusterer.prediction_data_.exemplars,
            clusterer.prediction_data_.dist_metric,
        )
        outlier_vec = outlier_membership_vector(
            nearest_neighbor,
            lambda_,
            clusters,
            clusterer.condensed_tree_._raw_tree,
            clusterer.prediction_data_.leaf_max_lambdas,
            clusterer.prediction_data_.cluster_tree,
        )

        result[i] = distance_vec**0.5 * outlier_vec**2.0
        result[i] /= result[i].sum()

        result[i] *= prob_in_some_cluster(
            nearest_neighbor,
            lambda_,
            clusters,
            clusterer.condensed_tree_._raw_tree,
            clusterer.prediction_data_.leaf_max_lambdas,
            clusterer.prediction_data_.cluster_tree,
        )

    return result


def all_points_membership_vectors(clusterer):
    """
    Predict soft cluster membership for all points in the original dataset.

    This function is more efficient by making use of the fact that all points
    are already in the condensed tree, and processing in bulk.

    Parameters
    ----------
    clusterer : HDBSCAN
         A clustering object that has been fit to the data and
        either had ``prediction_data=True`` set, or called the
        ``generate_prediction_data`` method after the fact.
        This method does not work if the clusterer was trained
        with ``metric='precomputed'``.

    Returns
    -------
    membership_vectors : array (n_samples, n_clusters)
        The probability that point ``i`` of the original dataset is a member of
        cluster ``j`` is in ``membership_vectors[i, j]``.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.approximate_predict : Predict the
        cluster label of new points.
    sklearn.cluster.hdbscan.prediction.membership_vectors : Predict soft cluster
        membership.
    """
    clusters = np.array(
        sorted(list(clusterer.condensed_tree_._select_clusters()))
    ).astype(np.intp)
    all_points = clusterer.prediction_data_.raw_data

    # When no clusters found, return array of 0's
    if clusters.size == 0:
        return np.zeros(all_points.shape[0])

    distance_vecs = all_points_dist_membership_vector(
        all_points,
        clusterer.prediction_data_.exemplars,
        clusterer.prediction_data_.dist_metric,
    )
    outlier_vecs = all_points_outlier_membership_vector(
        clusters,
        clusterer.condensed_tree_._raw_tree,
        clusterer.prediction_data_.leaf_max_lambdas,
        clusterer.prediction_data_.cluster_tree,
    )
    in_cluster_probs = all_points_prob_in_some_cluster(
        clusters,
        clusterer.condensed_tree_._raw_tree,
        clusterer.prediction_data_.leaf_max_lambdas,
        clusterer.prediction_data_.cluster_tree,
    )

    result = distance_vecs * outlier_vecs
    row_sums = result.sum(axis=1)
    result = result / row_sums[:, np.newaxis]
    result *= in_cluster_probs[:, np.newaxis]

    return result
