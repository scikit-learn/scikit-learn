"""flat.py

Provides alternative functions to hdbscan.HDBSCAN and others to
1. Allow prediction on a flat clustering by specifying 'n_clusters'.
    This is done by choosing the best cluster_selection_epsilon that produces
    the required number of clusters without adding unnecessary outliers.
2. Makes approximate_predict, membership_vector, and
    all_points_membership_vectors consistent with cluster_selection_epsilon

Provides the following functions:
==================================
HDBSCAN_flat: trained HDBSCAN instance with 'n_clusters' clusters
    The attributes (labels, probabilities, prediction_data) are tuned to
    produce 'n_clusters' clusters.

approximate_predict_flat: labels and probabilities for novel points
    Allows selecting n_clusters for novel points, or using the
    original clustering (potentially specified using cluster_selection_epsilon)

membership_vector_flat: Soft-clustering probabilities for novel points
    Similar to approximate_predict_flat, but for soft-clustering.
    **Use with caution**

all_points_membership_vectors_flat: Soft-clustering probabilities
    Similar to membership_vector_flat, but for points in training set
    **Use with caution**
"""

import copy
from warnings import warn

import numpy as np
from ._hdbscan_tree import compute_stability, get_cluster_tree_leaves
from .hdbscan_ import HDBSCAN, _tree_to_labels
from ._trees import _bfs_from_cluster_tree
from .prediction import (
    PredictionData,
    _find_cluster_and_probability,
    _find_neighbor_and_lambda,
)
from ._prediction_utils import (
    get_tree_row_with_child,
    dist_membership_vector,
    outlier_membership_vector,
    prob_in_some_cluster,
    all_points_dist_membership_vector,
    all_points_outlier_membership_vector,
    all_points_prob_in_some_cluster,
)


def HDBSCAN_flat(
    X,
    n_clusters=None,
    cluster_selection_epsilon=0.0,
    clusterer=None,
    inplace=False,
    **kwargs,
):
    """
    Train a HDBSCAN clusterer by specifying `n_clusters`.

    Or, modify a trained clusterer to return specific `n_clusters`.

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)
        A feature array, or array of distances between samples if
        `metric='precomputed'`.

    n_clusters : int, default=None
        Number of clusters to produce. If `None`, revert to default `HDBSCAN`.

    cluster_selection_epsilon : float, default=0
        Core-distance below which to stop splitting clusters. This can
        indirectly impose `n_clusters`. This argument is ignored if
        `n_clusters` is supplied.

    clusterer : HDBSCAN, default=None
        If supplied, modify this clusterer to produce `n_clusters` clusters.

    inplace : bool, default=False
        If 'clusterer' parameter is supplied, and `inplace=True`, modify
        `clusterer` inplace. If `inplace=False`, return a modified copy of
        `clusterer`.

    **kwargs : keyword arguments
        All keyword arguments to pass to `HDBSCAN`.

    Returns
    -------
    new_clusterer : HDBSCAN
        New `HDBSCAN` instance; returned irrespective of `inplace`.

    Examples
    --------
    >>> from sklearn.cluster import HDBSCAN, HDBSCAN_flat
    >>> from sklearn.datasets import make_blobs
    >>> from sklearn.utils import shuffle
    >>> from sklearn.preprocessing import StandardScaler
    >>>
    >>> X, y = make_blobs(n_samples=200, random_state=10)
    >>> X, y = shuffle(X, y, random_state=7)
    >>> X = StandardScaler().fit_transform(X)
    >>>
    >>> # Extract flat clustering from HDBSCAN's hierarchy for 7 clusters
    >>> clusterer = HDBSCAN_flat(X, n_clusters=7,
    ...                             min_cluster_size=12, min_samples=8)
    >>> labels = clusterer.labels_
    >>> proba = clusterer.probabilities_
    >>>
    >>> # Use a previously initialized/trained HDBSCAN
    >>> old_clusterer = HDBSCAN(min_cluster_size=12, min_samples=8)
    >>> clusterer = HDBSCAN_flat(X, n_clusters=7,
    ...                             clusterer=old_clusterer, inplace=True)
    >>> labels = clusterer.labels_
    >>> proba = clusterer.probabilities_

    See Also
    ---------
    sklearn.cluster.hdbscan.HDBSCAN: Perform HDBSCAN clustering from vector
        array or distance matrix.
    sklearn.cluster.hdbscan.flat.re_init: Modify PredictionData of HDBSCAN to
        account for epsilon.
    """
    # Handle the trivial case first.
    if (n_clusters is None) and (cluster_selection_epsilon == 0.0):
        if (not isinstance(clusterer, HDBSCAN)) or (not inplace):
            # Always generate prediction_data to avoid later woes
            kwargs["prediction_data"] = True
            new_clusterer = HDBSCAN(**kwargs)
        else:
            new_clusterer = clusterer
            new_clusterer.prediction_data = True

        new_clusterer.fit(X)
        return new_clusterer

    if (n_clusters is not None) and (cluster_selection_epsilon != 0.0):
        warn(
            f"'cluster_selection_epsilon' (={cluster_selection_epsilon})"
            " is ignored when 'n_clusters' is supplied."
        )
        cluster_selection_epsilon = 0.0
        # This will later be chosen according to n_clusters

    if not isinstance(clusterer, HDBSCAN):
        # Initialize and train clusterer if one was not previously supplied.
        # Always generate prediction data
        kwargs["prediction_data"] = True
        new_clusterer = HDBSCAN(**kwargs)
        # We do not pass cluster_selection_epsilon here.
        # While this adds unnecessary computation, it makes the code
        #   easier to read and debug.
        new_clusterer.fit(X)
    else:
        if inplace:
            new_clusterer = clusterer
        else:
            new_clusterer = copy.deepcopy(clusterer)

        new_clusterer.prediction_data = True

        # Train on 'X'. Do this even if the supplied clusterer was trained,
        #   because we want to make sure it fits 'X'.
        new_clusterer.prediction_data = True
        new_clusterer.fit(X)

    if new_clusterer.cluster_selection_method == "eom":
        max_eom_clusters = len(new_clusterer.condensed_tree_._select_clusters())

    # Pick an epsilon value right after a split produces n_clusters,
    #   and the don't split further for smaller epsilon (larger lambda)
    if n_clusters is not None:
        if (new_clusterer.cluster_selection_method == "eom") and (
            n_clusters > max_eom_clusters
        ):
            warn(
                f"Cannot predict more than {max_eom_clusters} with cluster "
                "selection method 'eom'. Changing to method 'leaf'..."
            )
            new_clusterer.cluster_selection_method = "leaf"
        epsilon = select_epsilon(new_clusterer.condensed_tree_, n_clusters)
    else:
        # Or use the specified cluster_selection_epsilon
        epsilon = cluster_selection_epsilon

    new_clusterer.cluster_selection_epsilon = float(epsilon)

    # Extract tree related stuff, in order to re-assign labels
    single_linkage_tree = new_clusterer.single_linkage_tree_
    single_linkage_tree = single_linkage_tree.to_numpy()
    min_cluster_size = new_clusterer.min_cluster_size
    cluster_selection_method = new_clusterer.cluster_selection_method
    allow_single_cluster = new_clusterer.allow_single_cluster
    match_reference_implementation = False

    # Get labels according to the required cluster_selection_epsilon
    output = _tree_to_labels(
        None,
        single_linkage_tree,
        min_cluster_size,
        cluster_selection_method,
        allow_single_cluster,
        match_reference_implementation,
        cluster_selection_epsilon=epsilon,
    )

    # Reflect the related changes in HDBSCAN.
    (
        new_clusterer.labels_,
        new_clusterer.probabilities_,
        new_clusterer.cluster_persistence_,
        new_clusterer._condensed_tree,
        new_clusterer._single_linkage_tree,
    ) = output

    # PredictionData attached to HDBSCAN should also change.
    # A function re_init is defined in this module to handle this.
    re_init(
        new_clusterer.prediction_data_,
        new_clusterer.condensed_tree_,
        cluster_selection_epsilon=epsilon,
    )
    return new_clusterer


def approximate_predict_flat(
    clusterer,
    points_to_predict,
    n_clusters=None,
    cluster_selection_epsilon=None,
    prediction_data=None,
    return_prediction_data=False,
):
    """
    Predict the cluster label of new points at a particular flat clustering.

    The clustering produced is specified by `n_clusters`. This is a modified
    version of `hdbscan.approximate_predict` to allow selection of
    `n_clusters`.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and either had
        `prediction_data=True` set, or called the `generate_prediction_data`
        method after the fact.

    points_to_predict : array, or array-like (n_samples, n_features)
        The new data points to predict cluster labels for. They should
        have the same dimensionality as the original dataset over which
        `clusterer` was fit.

    n_clusters : int, default=None
        The number of clusters to have in the flat clustering
            (over the training data, not points_to_predict)
        Ignored when prediction_data is supplied.

    cluster_selection_epsilon : float, default=None
        Core-distance below which to stop splitting clusters. This can
        indirectly impose `n_clusters`. This argument is ignored if
        `n_clusters` is supplied.

    prediction_data : PredictionData, default=None
        If supplied, use this to predict clusters for points_to_predict.
        This allows predicting on multiple datasets without corrupting
        prediction data associated with `clusterer`.

        If neither `n_clusters`, nor `prediction_data` are supplied,
            then the `prediction_data` associated with `clusterer` is used.

    return_prediction_data : bool, default=False
        If True, return `prediction_data` along with labels and proba.

    Returns
    -------
    labels : array (n_samples,)
        The predicted labels of the ``points_to_predict``.

    probabilities : array (n_samples,)
        The soft cluster scores for each of the ``points_to_predict``.

    prediction_data : PredictionData, optional
        The `prediction_data` used to predict. Returned if
        `return_prediciton_data=True`.

    Examples
    --------
    >>> from sklearn.cluster import HDBSCAN, approximate_predict_flat
    >>> from sklearn.datasets import make_blobs
    >>> from sklearn.utils import shuffle
    >>> from sklearn.preprocessing import StandardScaler
    >>>
    >>> X, y = make_blobs(n_samples=200, random_state=10)
    >>> X, y = shuffle(X, y, random_state=7)
    >>> X = StandardScaler().fit_transform(X)
    >>>
    >>> hdb = HDBSCAN(prediction_data=True)
    >>> hdb.fit(X)
    HDBSCAN(prediction_data=True)
    >>> # From a fitted HDBSCAN model, predict for n_clusters=5
    >>> labels, proba = approximate_predict_flat(
    ...                     hdb, X, n_clusters=5)
    >>>
    >>> # Store prediciton data for later use.
    >>> labels, proba, pred_data = approximate_predict_flat(
    ...                             hdb, X, n_clusters=5,
    ...                             return_prediction_data=True)
    >>>
    >>> # Use this prediction data to predict on new points
    >>> labels1, proba1 = approximate_predict_flat(
    ...                             hdb, X,
    ...                             prediction_data=pred_data)
    See Also
    ---------
    sklearn.cluster.hdbscan.prediction.approximate_predict : Predict the
        cluster label of new points.
    """
    # Get number of fitted clusters for later use.
    n_clusters_fit = np.sum(np.unique(clusterer.labels_) >= 0)
    if n_clusters is not None:
        n_clusters = int(n_clusters)  # Ensure n_clusters is int

    # We'll need the condensed tree later...
    condensed_tree = clusterer.condensed_tree_

    # If none of the three arguments: prediction_data, n_clusters,
    #   and cluster_selection_epsilon are supplied,
    # then use clusterer's prediciton data directly
    if (
        (prediction_data is None)
        and ((n_clusters is None) or (n_clusters == n_clusters_fit))
        and (cluster_selection_epsilon is None)
    ):
        prediction_data = clusterer.prediction_data_

    # If either of n_clusters or cluster_selection_epsilon were supplied,
    #   then build prediction data from these by modifying clusterer's
    if not isinstance(prediction_data, PredictionData):
        if clusterer.prediction_data_ is None:
            raise ValueError(
                "Clusterer does not have prediction data!"
                " Try fitting with prediction_data=True set,"
                " or run generate_prediction_data on the clusterer"
            )
        # Get prediction data from clusterer
        prediction_data = clusterer.prediction_data_
        # Modify prediction_data to reflect new n_clusters
        # First, make a copy of prediction data to avoid modifying source
        prediction_data = copy.deepcopy(prediction_data)
        # Cluster selection method is hold by condensed_tree.
        # Change from 'eom' to 'leaf' if n_clusters is too large.
        if (condensed_tree.cluster_selection_method == "eom") and (
            (n_clusters is not None) and (n_clusters > n_clusters_fit)
        ):
            warn(
                f"Cannot predict more than {n_clusters_fit} with cluster "
                "selection method 'eom'. Changing to method 'leaf'..."
            )
            condensed_tree.cluster_selection_method = "leaf"
        # This change does not affect the tree associated with 'clusterer'
        # Re-initialize prediction_data for the specified n_clusters or epsilon
        re_init(
            prediction_data,
            condensed_tree,
            n_clusters=n_clusters,
            cluster_selection_epsilon=cluster_selection_epsilon,
        )

    # ============================================================
    # Now we're ready to use prediction_data
    # The rest of the code is copied from HDBSCAN's approximate_predict,
    #   but modified to use prediction_data instead of clusterer's attribute
    points_to_predict = np.asarray(points_to_predict)

    if points_to_predict.shape[1] != prediction_data.raw_data.shape[1]:
        raise ValueError("New points dimension does not match fit data!")

    if prediction_data.cluster_tree.shape[0] == 0:
        warn(
            "Prediction data does not have any defined clusters, new data"
            " will be automatically predicted as noise."
        )
        labels = -1 * np.ones(points_to_predict.shape[0], dtype=np.int32)
        probabilities = np.zeros(points_to_predict.shape[0], dtype=np.float32)
        if return_prediction_data:
            return labels, probabilities, prediction_data
        else:
            return labels, probabilities

    labels = np.empty(points_to_predict.shape[0], dtype=np.int32)
    probabilities = np.empty(points_to_predict.shape[0], dtype=np.float64)

    min_samples = clusterer.min_samples or clusterer.min_cluster_size
    neighbor_distances, neighbor_indices = prediction_data.tree.query(
        points_to_predict, k=2 * min_samples
    )

    for i in range(points_to_predict.shape[0]):
        label, prob = _find_cluster_and_probability(
            condensed_tree,
            prediction_data.cluster_tree,
            neighbor_indices[i],
            neighbor_distances[i],
            prediction_data.core_distances,
            prediction_data.cluster_map,
            prediction_data.max_lambdas,
            min_samples,
        )
        labels[i] = label
        probabilities[i] = prob

    if return_prediction_data:
        return labels, probabilities, prediction_data
    else:
        return labels, probabilities


def membership_vector_flat(
    clusterer,
    points_to_predict,
    prediction_data=None,
    n_clusters=None,
    cluster_selection_epsilon=0.0,
):
    """
    Predict soft cluster membership probabilities.

    Produces a vector for each point in ``points_to_predict`` that gives a
    probability that the given point is a member of a cluster for each of the
    selected clusters of the ``clusterer``.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and either had
        `prediction_data=True` set, or called the `generate_prediction_data`
        method after the fact.

    points_to_predict : array, or array-like (n_samples, n_features)
        The new data points to predict cluster labels for. They should
        have the same dimensionality as the original dataset over which
        clusterer was fit.

    prediction_data : PredictionData, default=None
        Prediction data associated with HDBSCAN for some flat clustering.

    n_clusters : int, default=None
        Number of clusters over which to compute membership probabilities.
        These clusters are obtained as a flat clustering at some
        `cluster_selection_epsilon`.

    cluster_selection_epsilon : float, default=0
        Core-distance below which to stop splitting clusters. This can
        indirectly impose `n_clusters`. This argument is ignored if
        `n_clusters` is supplied.

    Returns
    -------
    membership_vectors : array (n_samples, n_clusters)
        The probability that point ``i`` is a member of cluster ``j`` is
        in ``membership_vectors[i, j]``.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.membership_vectors : Predict soft cluster
        membership.
    sklearn.cluster.hdbscan.prediction.all_points_membership_vectors : Predict
        soft cluster membership vectors for all points in the original dataset
        the clusterer was trained on.

    Notes
    -----
    This function is an adaptation of hdbscan's membership_vector for
    `n_clusters`, `epsilon`. If neither `n_clusters` nor
    `cluster_selection_epsilon` are supplied, the `clusterer`'s original
    clustering is used.
    """
    points_to_predict = points_to_predict.astype(np.float64)
    # Extract condensed tree for later use
    condensed_tree = clusterer.condensed_tree_

    # Choose flat clustering based on cluster_selection_epsilon or n_clusters.
    # If neither is specified, use clusterer's cluster_selection_epsilon
    if (
        (n_clusters is None)
        and (cluster_selection_epsilon == 0.0)
        and (prediction_data is None)
    ):
        epsilon = clusterer.cluster_selection_epsilon
        # Use the same prediction_data as clusterer's
        prediction_data = clusterer.prediction_data_
    elif prediction_data is None:
        if n_clusters is not None:
            # Compute cluster_selection_epsilon so that a flat clustering
            #   produces a specified number of n_clusters
            # With method 'eom', we may fail to get 'n_clusters' clusters. So,
            try:
                epsilon = select_epsilon(condensed_tree, n_clusters)
            except AssertionError:
                warn(
                    f"Failed to predict {n_clusters} clusters with "
                    "cluster selection method 'eom'. Switching to 'leaf'..."
                )
                condensed_tree.cluster_selection_method = "leaf"
                epsilon = select_epsilon(condensed_tree, n_clusters)
        else:
            epsilon = cluster_selection_epsilon
        # Create another instance of prediction_data that is consistent
        #   with the selected value of epsilon.
        prediction_data = copy.deepcopy(clusterer.prediction_data_)
        re_init(prediction_data, condensed_tree, cluster_selection_epsilon=epsilon)

    # Flat clustering from prediction data
    clusters = clusters_from_prediction_data(prediction_data)

    # Initialize probabilities
    result = np.empty((points_to_predict.shape[0], clusters.shape[0]), dtype=np.float64)

    # k-NN for prediciton points to training set
    min_samples = clusterer.min_samples or clusterer.min_cluster_size
    neighbor_distances, neighbor_indices = prediction_data.tree.query(
        points_to_predict, k=2 * min_samples
    )

    # Loop over prediction points to compute probabilities
    for i in range(points_to_predict.shape[0]):
        # We need to find where in the tree the new point would go
        # for the purposes of outlier membership approximation
        nearest_neighbor, lambda_ = _find_neighbor_and_lambda(
            neighbor_indices[i],
            neighbor_distances[i],
            prediction_data.core_distances,
            min_samples,
        )

        # Find row in tree where nearest neighbor drops out,
        #   so we can get a lambda value for the nearest neighbor
        neighbor_tree_row = get_tree_row_with_child(
            condensed_tree._raw_tree, nearest_neighbor
        )

        # Assign lambda as min(lambda-to-neighbor, neighbor's-lambda-to-tree)
        # Equivalently, this assigns core distance for prediction point as
        #   max(dist-to-neighbor, neighbor's-dist-to-tree)
        if neighbor_tree_row["lambda_val"] <= lambda_:
            lambda_ = neighbor_tree_row["lambda_val"]

        # Probabilities based on distance to closest exemplar in each cluster:
        # Use new prediction_data that points to exemplars that are specific
        #   to the choice of n_clusters
        distance_vec = dist_membership_vector(
            points_to_predict[i], prediction_data.exemplars, prediction_data.dist_metric
        )
        # Probabilities based on how long the nearest exemplar persists in
        #   each cluster (with respect to most persistent exemplar)
        # Use new clusters that are defined by the choice of n_clusters.
        outlier_vec = outlier_membership_vector(
            nearest_neighbor,
            lambda_,
            clusters,
            condensed_tree._raw_tree,
            prediction_data.leaf_max_lambdas,
            prediction_data.cluster_tree,
        )

        # Merge the two probabilities to produce a single set of probabilities
        result[i] = distance_vec**0.5 * outlier_vec**2.0
        result[i] /= result[i].sum()

        # Include probability that the nearest neighbor belongs to a cluster
        result[i] *= prob_in_some_cluster(
            nearest_neighbor,
            lambda_,
            clusters,
            condensed_tree._raw_tree,
            prediction_data.leaf_max_lambdas,
            prediction_data.cluster_tree,
        )

    # Rename variable so it's easy to understand what's being returned
    membership_vectors = result
    return membership_vectors


def all_points_membership_vectors_flat(
    clusterer, prediction_data=None, n_clusters=None, cluster_selection_epsilon=None
):
    """
    Predict soft cluster membership vectors for all points in the dataset.

    This function predicts soft cluster membership vectors for all the points
    in the dataset that the clusterer was trained on. This function is more
    efficient by making use of the fact that all points are already in the
    condensed tree, and processing in bulk.

    Parameters
    ----------
    clusterer : HDBSCAN
        A clustering object that has been fit to the data and either had
        `prediction_data=True` set, or called the `generate_prediction_data`
        method after the fact. This method does not work if the clusterer was
        trained with `metric='precomputed'`.

    prediction_data : PredictionData, default=None
        Prediction data associated with HDBSCAN for some flat clustering.

    n_clusters : int, default=None
        Number of clusters over which to compute membership probabilities.
        These clusters are obtained as a flat clustering at some
        `cluster_selection_epsilon`.

    cluster_selection_epsilon : float, default=0
        Core-distance below which to stop splitting clusters. This can
        indirectly impose `n_clusters`. This argument is ignored if
        `n_clusters` is supplied.

    Returns
    -------
    membership_vectors : array (n_samples, n_clusters)
        The probability that point `i` of the original dataset is a member of
        cluster `j` is in `membership_vectors[i, j]`.

    See Also
    --------
    sklearn.cluster.hdbscan.prediction.all_points_membership_vectors : Predict
        soft cluster membership vectors for all points in the original dataset
        the clusterer was trained on.
    sklearn.cluster.hdbscan.prediction.membership_vectors : Predict soft cluster
        membership.

    Notes
    -----
    This function is an adaptation of hdbscan's `all_points_membership_vector`
    for `n_clusters`, `epsilon`. If neither `n_clusters` nor
    `cluster_selection_epsilon` are supplied, the `clusterer`'s original
    clustering is used.
    """
    # Extract condensed tree for later use
    condensed_tree = clusterer.condensed_tree_

    # Choose flat clustering based on cluster_selection_epsilon or n_clusters.
    # If neither is specified, use clusterer's cluster_selection_epsilon
    if (n_clusters is None) and (cluster_selection_epsilon is None):
        epsilon = clusterer.cluster_selection_epsilon
        # Use the same prediction_data as clusterer's
        prediction_data = clusterer.prediction_data_
    elif prediction_data is None:
        if n_clusters is not None:
            # Compute cluster_selection_epsilon so that a flat clustering
            #   produces a specified number of n_clusters
            # With method 'eom', we may fail to get 'n_clusters' clusters. So,
            try:
                epsilon = select_epsilon(condensed_tree, n_clusters)
            except AssertionError:
                warn(
                    f"Failed to predict {n_clusters} clusters with "
                    "cluster selection method 'eom'. Switching to 'leaf'..."
                )
                condensed_tree.cluster_selection_method = "leaf"
                epsilon = select_epsilon(condensed_tree, n_clusters)
        else:
            epsilon = cluster_selection_epsilon
        # Create another instance of prediction_data that is consistent
        #   with the selected value of epsilon.
        prediction_data = copy.deepcopy(clusterer.prediction_data_)
        re_init(prediction_data, condensed_tree, cluster_selection_epsilon=epsilon)

    # Flat clustering at the chosen epsilon from prediction_data
    clusters = clusters_from_prediction_data(prediction_data)

    all_points = prediction_data.raw_data

    # When no clusters found, return array of 0's
    if clusters.size == 0:
        return np.zeros(all_points.shape[0])

    # Probabilities based on distance to closest exemplar in each cluster:
    # Use new prediction_data that points to exemplars that are specific
    #   to the choice of n_clusters
    distance_vecs = all_points_dist_membership_vector(
        all_points, prediction_data.exemplars, prediction_data.dist_metric
    )

    # Probabilities based on how long the point persists in
    #   each cluster (with respect to most persistent exemplar)
    # Use new clusters that are defined by the choice of n_clusters.
    outlier_vecs = all_points_outlier_membership_vector(
        clusters,
        condensed_tree._raw_tree,
        prediction_data.leaf_max_lambdas,
        prediction_data.cluster_tree,
    )

    # Include probability that the point belongs to a cluster
    in_cluster_probs = all_points_prob_in_some_cluster(
        clusters,
        condensed_tree._raw_tree,
        prediction_data.leaf_max_lambdas,
        prediction_data.cluster_tree,
    )

    # Aggregate the three probabilities to produce membership vectors
    result = distance_vecs * outlier_vecs
    row_sums = result.sum(axis=1)
    result = result / row_sums[:, np.newaxis]
    result *= in_cluster_probs[:, np.newaxis]

    # Re-name variable to clarify what's being returned.
    membership_vectors = result
    return membership_vectors


def select_epsilon(condensed_tree, n_clusters):
    """
    Pick optimal epsilon from condensed tree based on n_clusters,
        calls functions specific to 'eom' or 'leaf' selection methods
    """
    cluster_selection_method = condensed_tree.cluster_selection_method
    if cluster_selection_method == "eom":
        return select_epsilon_eom(condensed_tree, n_clusters)
    if cluster_selection_method == "leaf":
        return select_epsilon_leaf(condensed_tree, n_clusters)
    raise ValueError(
        'Invalid Cluster Selection Method: %s\nShould be one of: "eom", "leaf"\n'
    )


def select_epsilon_eom(condensed_tree, n_clusters):
    """
    Select epsilon so that persistence-based clustering,
        after truncating the tree at the above epsilon,
        has exactly 'n_clusters' clusters
    """
    # With method 'eom', max clusters are produced for epsilon=0,
    #   as computed by
    eom_base_clusters = condensed_tree._select_clusters()
    max_clusters = len(eom_base_clusters)
    # Increasing epsilon can only reduce the number of ouput clusters.

    assert n_clusters <= max_clusters, (
        f"Cannot produce more than {max_clusters} with method 'eom'. "
        + "Use method 'leaf' instead to extract flat clustering."
    )

    tree = condensed_tree._raw_tree
    # To select epsilon, consider all values where clusters are split
    cluster_lambdas = tree["lambda_val"][tree["child_size"] > 1]
    candidate_epsilons = 1.0 / np.unique(cluster_lambdas) - 1.0e-12
    # Subtract the extra e-12 to avoid numerical errors in comparison
    # Then, we avoid splitting for all epsilon below this.
    candidate_epsilons = np.sort(candidate_epsilons)[::-1]

    for epsilon in candidate_epsilons:
        sel_clusters = _new_select_clusters(condensed_tree, epsilon)
        if len(sel_clusters) == n_clusters:
            break
    else:
        raise RuntimeError("Could not find epsilon")

    return epsilon


def select_epsilon_leaf(condensed_tree, n_clusters):
    """
    Select epsilon so that the leaves of condensed tree,
        after truncating at the above epsilon,
        has exactly 'n_clusters' clusters
    """
    # Use an epsilon value that produces the right number of clusters.
    # The condensed tree of HDBSCAN has this information.
    # Extract the lambda levels (=1/distance) from the condensed tree
    lambdas = condensed_tree._raw_tree["lambda_val"]
    # We don't want values that produce a large cluster and
    #   just one or two individual points.
    child_sizes = condensed_tree._raw_tree["child_size"]
    child_sizes = child_sizes.astype(int)
    # Keep only those lambda values corresponding to cluster separation;
    #   i.e., with child_sizes > 1
    lambdas = lambdas[child_sizes > 1]
    # Get the unique values, because when two clusters fall out of one,
    #   the entry with lambda is repeated.
    lambdas = np.unique(lambdas.astype(float))
    if n_clusters > len(lambdas) + 1:
        warn(
            f"HDBSCAN can only compute {len(lambdas)+1} clusters. "
            f"Setting n_clusters to {len(lambdas)+1}..."
        )
        n_clusters = len(lambdas) + 1

    # lambda values are sorted by np.unique.
    # Now, get epsilon (distance threshold) as 1/lambda
    epsilon = 1.0 / lambdas[n_clusters - 2]
    # At this epsilon, n_clusters have been split.
    # Stop splits at epsilons smaller than this.
    # To allow for numerical errors,
    return epsilon - 1.0e-12


def re_init(predData, condensed_tree, n_clusters=None, cluster_selection_epsilon=0.0):
    """
    Modify PredictionData of HDBSCAN to account for epsilon.
    epsilon is the cluster_selection_epsilon that controls granularity
        of clusters; Large epsilon => More clusters

    Parameters
    ----------
    predData: PredictionData
        Contains data to use for predicting novel points.
        Defined in the HDBSCAN module

    condensed_tree: CondensedTree
        Tree structure that contains hierarchical clustering.
        Defined in the HDBSCAN module

    n_clusters: int, optional, default=None
        If specified, use this to obtain cluster_selection_epsilon
            from CondensedTree; Overrides cluster_selection_epsilon parameter

    cluster_selection_epsilon: float, default=0.
        In cluster tree, nodes are not split further beyond (>=) this value.
        epsilon is the inverse of core distance.

    Returns
    -------
    None
    """
    # predData must be a pre-trained PredictionData instance from hdbscan
    # If n_clusters is specified, compute cluster_selection_epsilon;
    if n_clusters is not None:
        cluster_selection_epsilon = select_epsilon(condensed_tree, n_clusters)

    # This is the key modification:
    # Select clusters according to selection method and epsilon.
    selected_clusters = _new_select_clusters(condensed_tree, cluster_selection_epsilon)
    # _new_select_clusters is a modification of get_clusters
    #   from hdbscan._hdbscan_tree

    # raw tree, used later to get exemplars and lambda values
    raw_condensed_tree = condensed_tree._raw_tree

    # Re-do the cluster map: Map cluster numbers in tree (N, N+1, ..)
    #    to the cluster labels produced as output
    predData.cluster_map = {
        int(c): n for n, c in enumerate(sorted(list(selected_clusters)))
    }
    predData.reverse_cluster_map = {n: c for c, n in predData.cluster_map.items()}

    # Re-compute lambdas and exemplars for selected clusters;
    predData.max_lambdas = {}
    predData.exemplars = []

    for cluster in selected_clusters:
        # max_lambda <=> smallest distance <=> most persistent point(s)
        predData.max_lambdas[cluster] = raw_condensed_tree["lambda_val"][
            raw_condensed_tree["parent"] == cluster
        ].max()

        # Map all sub-clusters of selected cluster to the selected cluster's
        #       label in output.
        # Map lambdas too...
        for sub_cluster in predData._clusters_below(cluster):
            predData.cluster_map[sub_cluster] = predData.cluster_map[cluster]
            predData.max_lambdas[sub_cluster] = predData.max_lambdas[cluster]

        # Create set of exemplar points for later use.
        # Novel points are assigned based on cluster of closest exemplar.
        cluster_exemplars = np.array([], dtype=np.int64)
        # For each selected cluster, get all of its leaves,
        #       and leaves of leaves, and so on...
        for leaf in predData._recurse_leaf_dfs(cluster):
            # Largest lambda => Most persistent points
            leaf_max_lambda = raw_condensed_tree["lambda_val"][
                raw_condensed_tree["parent"] == leaf
            ].max()
            # Get the most persistent points
            points = raw_condensed_tree["child"][
                (raw_condensed_tree["parent"] == leaf)
                & (raw_condensed_tree["lambda_val"] == leaf_max_lambda)
            ]
            # Add most persistent points as exemplars
            cluster_exemplars = np.hstack([cluster_exemplars, points])

        # Add exemplars for each leaf of each selected cluster.
        predData.exemplars.append(predData.raw_data[cluster_exemplars])
    return


def _new_select_clusters(
    condensed_tree,
    cluster_selection_epsilon,
    allow_single_cluster=False,
    match_reference_implementation=False,
):
    """
    Adaptation of get_clusters from hdbscan._hdbscan_tree.
    Avoids the label and proba computation at the end,
        and returns only the selected clusters instead.
    """
    tree = condensed_tree._raw_tree
    cluster_selection_method = condensed_tree.cluster_selection_method
    stability = compute_stability(tree)

    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        node_list = sorted(stability.keys(), reverse=True)[:-1]
        # (exclude root)

    cluster_tree = tree[tree["child_size"] > 1]
    is_cluster = {cluster: True for cluster in node_list}

    if cluster_selection_method == "eom":
        for node in node_list:
            child_selection = cluster_tree["parent"] == node
            subtree_stability = np.sum(
                [stability[child] for child in cluster_tree["child"][child_selection]]
            )
            if subtree_stability > stability[node]:
                is_cluster[node] = False
                stability[node] = subtree_stability
            else:
                for sub_node in _bfs_from_cluster_tree(cluster_tree, node):
                    if sub_node != node:
                        is_cluster[sub_node] = False

        if cluster_selection_epsilon != 0.0:
            eom_clusters = set([c for c in is_cluster if is_cluster[c]])
            selected_clusters = epsilon_search(
                eom_clusters,
                cluster_tree,
                cluster_selection_epsilon,
                allow_single_cluster,
            )
            for c in is_cluster:
                if c in selected_clusters:
                    is_cluster[c] = True
                else:
                    is_cluster[c] = False

    elif cluster_selection_method == "leaf":
        leaves = set(get_cluster_tree_leaves(cluster_tree))
        if len(leaves) == 0:
            for c in is_cluster:
                is_cluster[c] = False
            is_cluster[tree["parent"].min()] = True

        if cluster_selection_epsilon != 0.0:
            selected_clusters = epsilon_search(
                leaves, cluster_tree, cluster_selection_epsilon, allow_single_cluster
            )
        else:
            selected_clusters = leaves

        for c in is_cluster:
            if c in selected_clusters:
                is_cluster[c] = True
            else:
                is_cluster[c] = False
    else:
        raise ValueError(
            'Invalid Cluster Selection Method: %s\nShould be one of: "eom", "leaf"\n'
        )

    clusters = set([int(c) for c in is_cluster if is_cluster[c]])
    return clusters


def epsilon_search(
    leaves, cluster_tree, cluster_selection_epsilon, allow_single_cluster
):
    selected_clusters = []
    processed = []

    for leaf in leaves:
        eps = 1 / cluster_tree["lambda_val"][cluster_tree["child"] == leaf][0]
        if eps < cluster_selection_epsilon:
            if leaf not in processed:
                epsilon_child = traverse_upwards(
                    cluster_tree, cluster_selection_epsilon, leaf, allow_single_cluster
                )
                if hasattr(epsilon_child, "__len__"):
                    epsilon_child = epsilon_child[0]

                selected_clusters.append(epsilon_child)

                for sub_node in _bfs_from_cluster_tree(cluster_tree, epsilon_child):
                    if sub_node != epsilon_child:
                        processed.append(sub_node)
        else:
            selected_clusters.append(leaf)

    return set(selected_clusters)


def traverse_upwards(
    cluster_tree, cluster_selection_epsilon, leaf, allow_single_cluster
):
    root = cluster_tree["parent"].min()
    parent = cluster_tree[cluster_tree["child"] == leaf]["parent"]
    if parent == root:
        if allow_single_cluster:
            return parent
        else:
            return leaf  # return node closest to root

    parent_eps = 1 / cluster_tree[cluster_tree["child"] == parent]["lambda_val"]
    if parent_eps > cluster_selection_epsilon:
        return parent
    else:
        return traverse_upwards(
            cluster_tree, cluster_selection_epsilon, parent, allow_single_cluster
        )


def clusters_from_prediction_data(prediction_data):
    """
    Extract selected clusters from PredictionData instance.
    """
    return np.array(sorted(list(prediction_data.reverse_cluster_map.values()))).astype(
        np.intp
    )
