import numpy as np
import scipy.sparse.csgraph
import scipy.spatial.distance
import scipy.stats


def compute_pair_to_pair_dists(X, metric):
    """
    Computes the pairwise distance matrix between samples in the input array
      `X` using the specified distance metric.

    Parameters:
    - X (numpy.ndarray): Sample embeddings with shape (N, D).
    - metric (str): Distance metric to compute dissimilarity between observations.

    Returns:
    - dists (numpy.ndarray): Pairwise distance matrix with shape (N, N).
    """
    dists = scipy.spatial.distance.cdist(X, X, metric=metric)
    np.maximum(dists, 1e-12, out=dists)
    np.fill_diagonal(dists, val=np.inf)
    return dists


def get_subarray(arr, inds_a=None, inds_b=None):
    """
    Retrieves a subarray from the input array based on specified indices.

    Parameters:
    - arr (numpy.ndarray): Input array.
    - inds_a (numpy.ndarray, optional): Indices for the first dimension.
    Defaults to None.
    - inds_b (numpy.ndarray, optional): Indices for the second dimension.
    If None, defaults to `inds_a`.

    Returns:
    - subarray (numpy.ndarray): Subarray based on specified indices.
    """
    if inds_a is None:
        return arr
    if inds_b is None:
        inds_b = inds_a
    inds_a_mesh, inds_b_mesh = np.meshgrid(inds_a, inds_b)
    return arr[inds_a_mesh, inds_b_mesh]


def get_internal_objects(mutual_reach_dists):
    """
    Identifies internal nodes and corresponding edge weights in a Minimum
    Spanning Tree (MST) represented by `mutual_reach_dists`.

    Parameters:
    - mutual_reach_dists (numpy.ndarray): Matrix representing mutual
    reachability distances.

    Returns:
    - internal_node_inds (numpy.ndarray): Indices of internal nodes.
    - internal_edge_weights (numpy.ndarray): Edge weights corresponding
    to internal nodes.
    """
    mst = scipy.sparse.csgraph.minimum_spanning_tree(mutual_reach_dists)
    mst = mst.toarray()

    is_mst_edges = mst > 0.0

    internal_node_inds = (is_mst_edges + is_mst_edges.T).sum(axis=0) > 1
    internal_node_inds = np.flatnonzero(internal_node_inds)

    internal_edge_weights = get_subarray(mst, inds_a=internal_node_inds)

    return internal_node_inds, internal_edge_weights


def compute_cluster_core_distance(dists, d):
    """
    Computes the core distances for each sample in a given distance matrix.

    Parameters:
    - dists (numpy.ndarray): Pairwise distance matrix.
    - d (int): Exponent value for the core distance computation.

    Returns:
    - core_dists (numpy.ndarray): Core distances for each sample.
    """
    n, m = dists.shape

    if n == m and n > 800:
        from ...neighbors import NearestNeighbors

        nn = NearestNeighbors(n_neighbors=801, metric="precomputed")
        dists, _ = nn.fit(np.nan_to_num(dists, posinf=0.0)).kneighbors(
            return_distance=True
        )
        n = dists.shape[1]

    core_dists = np.power(dists, -d).sum(axis=-1, keepdims=True) / (n - 1 + 1e-12)

    np.clip(core_dists, a_min=1e-12, a_max=1e12, out=core_dists)

    np.power(core_dists, -1.0 / d, out=core_dists)

    return core_dists


def compute_mutual_reach_dists(
    dists, d, is_symmetric, cls_inds_a=None, cls_inds_b=None
):
    """
    Computes mutual reachability distances based on the given distance matrix
      and clustering indices.

    Parameters:
    - dists (numpy.ndarray): Pairwise distance matrix.
    - d (float): Exponent value for core distance computation.
    - is_symmetric (bool): Indicates whether the computation is for symmetric
    mutual reachability distances.
    - cls_inds_a (numpy.ndarray, optional): Indices for the first cluster.
    Defaults to None.
    - cls_inds_b (numpy.ndarray, optional): Indices for the second cluster.
    Defaults to None.

    Returns:
    - mutual_reach_dists (numpy.ndarray): Matrix of mutual reachability distances.
    """
    cls_dists = get_subarray(dists, inds_a=cls_inds_a, inds_b=cls_inds_b)

    if is_symmetric:
        core_dists_a = core_dists_b = compute_cluster_core_distance(
            d=d, dists=cls_dists
        )

    else:
        core_dists_a = compute_cluster_core_distance(d=d, dists=cls_dists)
        core_dists_b = compute_cluster_core_distance(d=d, dists=cls_dists.T).T

    mutual_reach_dists = cls_dists.copy()
    np.maximum(mutual_reach_dists, core_dists_a, out=mutual_reach_dists)
    np.maximum(mutual_reach_dists, core_dists_b, out=mutual_reach_dists)

    return mutual_reach_dists


def fn_density_sparseness(cls_inds, dists, d):
    """
    Computes the density sparseness of a cluster based on its indices and the
      pairwise distance matrix.

    Parameters:
    - cls_inds (numpy.ndarray): Indices of samples in the cluster.
    - dists (numpy.ndarray): Pairwise distance matrix.
    - d (int): Exponent value for core distance computation.

    Returns:
    - dsc (float): Density sparseness of the cluster.
    - internal_node_inds (numpy.ndarray): Indices of internal nodes in the cluster.
    """
    if cls_inds.size <= 3:
        return 0.0, np.empty(0, dtype=int)

    mutual_reach_dists = compute_mutual_reach_dists(dists=dists, d=d, is_symmetric=True)
    internal_node_inds, internal_edge_weights = get_internal_objects(mutual_reach_dists)

    dsc = float(internal_edge_weights.max())
    internal_node_inds = cls_inds[internal_node_inds]

    return dsc, internal_node_inds


def fn_density_separation(cls_i, cls_j, dists, d):
    """
    Computes the density separation between two clusters based on their
    indices and the pairwise distance matrix.

    Parameters:
    - cls_i (int): Cluster ID of the first cluster.
    - cls_j (int): Cluster ID of the second cluster.
    - dists (numpy.ndarray): Pairwise distance matrix.
    - d (int): Exponent value for core distance computation.

    Returns:
    - cls_i (int): Cluster ID of the first cluster.
    - cls_j (int): Cluster ID of the second cluster.
    - dspc_ij (float): Density separation between the two clusters.
    """
    mutual_reach_dists = compute_mutual_reach_dists(
        dists=dists, d=d, is_symmetric=False
    )
    dspc_ij = float(mutual_reach_dists.min()) if mutual_reach_dists.size else np.inf
    return cls_i, cls_j, dspc_ij


def _check_duplicated_samples(X, threshold=1e-9):
    """
    Checks for duplicated samples in the input array `X` based
    on a specified threshold.

    Parameters:
    - X (numpy.ndarray): Input array containing samples.
    - threshold (float, optional): Threshold for considering samples as duplicated.
      Defaults to 1e-9.

    Raises:
    - ValueError: If duplicated samples are found in `X`.
    """
    if X.shape[0] <= 1:
        return

    from ...neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(X)
    dists, _ = nn.kneighbors(return_distance=True)

    if np.any(dists < threshold):
        raise ValueError("Duplicated samples have been found in X.")
