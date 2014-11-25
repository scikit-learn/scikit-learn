# Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause
from __future__ import division

import warnings
import numpy as np
from scipy import sparse
from math import sqrt

from ..metrics.pairwise import euclidean_distances
from ..base import TransformerMixin, ClusterMixin, BaseEstimator
from ..externals.six.moves import xrange
from ..utils import check_array
from ..utils.extmath import row_norms, safe_sparse_dot
from .hierarchical import AgglomerativeClustering


def _iterate_sparse_X(X):
    """
    This little hack returns a densified row when iterating over a sparse
    matrix, insted of constructing a sparse matrix for every row that is
    expensive.
    """
    n_samples = X.shape[0]
    X_indices = X.indices
    X_data = X.data
    X_indptr = X.indptr

    for i in xrange(n_samples):
        row = np.zeros(X.shape[1])
        startptr, endptr = X_indptr[i], X_indptr[i + 1]
        nonzero_indices = X_indices[startptr: endptr]
        row[nonzero_indices] = X_data[startptr: endptr]
        yield row


class _CFNode(object):
    """Each node in a CFTree is called a CFNode.

    The CFNode can have a maximum of branching_factor
    number of CFSubclusters.

    Parameters
    ----------
    threshold : float
        Threshold needed for a new subcluster to enter a CFSubcluster.

    branching_factor : int
       Maximun number of CF subclusters in each node.

    is_leaf : bool
        We need to know if the CFNode is a leaf or not, in order to
        retrieve the final subclusters.

    Attributes
    ----------
    subclusters_ : array-like
        list of subclusters for a particular CFNode.

    prev_leaf_ : _CFNode
        prev_leaf. Useful only if is_leaf is True.

    next_leaf_ : _CFNode
        next_leaf. Useful only if is_leaf is True.
        the final subclusters.

    init_centroids_ : ndarray, shape (branching_factor + 1, n_features)
        manipulate init_centroids throughout rather than centroids_ since
        the centroids are just a view of the init_centroids_ .

    init_sq_norm_ : ndarray, shape (branching_factor + 1,)
        manipulate squared_norms throughout. similar to init_centroids_.

    centroids_ : ndarray
        view of init_centroids_.

    squared_norm_ : ndarray
        view of init_sq_norm_.

    n_ : int
        number of subclusters in the node.
    """
    def __init__(self, threshold, branching_factor, is_leaf, n_features):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.is_leaf = is_leaf
        self.n_features = n_features

        # The list of subclusters, centroids and squared norms
        # to manipulate throughout.
        self.subclusters_ = []
        self.init_centroids_ = np.zeros((branching_factor + 1, n_features))
        self.init_sq_norm_ = np.zeros((branching_factor + 1))
        self.n_ = 0
        self.squared_norm_ = []
        self.prev_leaf_ = None
        self.next_leaf_ = None

    def update(self, subcluster):
        self.subclusters_.append(subcluster)
        self.init_centroids_[self.n_] = subcluster.centroid_
        self.init_sq_norm_[self.n_] = subcluster.sq_norm_
        self.n_ += 1

        # Keep centroids and squared norm as views. In this way
        # if we change init_centroids and init_sq_norm_, it is
        # sufficient,
        self.centroids_ = self.init_centroids_[:self.n_, :]
        self.squared_norm_ = self.init_sq_norm_[:self.n_]

    def split_node(self, child_node, parent_subcluster):
        r"""
        Split the node if there is no place for a new subcluster.
        """
        dist = euclidean_distances(
            child_node.centroids_, Y_norm_squared=child_node.squared_norm_,
            squared=True)
        n_clusters = dist.shape[0]

        farthest_idx = np.unravel_index(
            dist.argmax(), (n_clusters, n_clusters))
        dist_idx = dist[[farthest_idx]]
        newsubcluster1 = _CFSubcluster()
        newsubcluster2 = _CFSubcluster()

        threshold = self.threshold
        branching_factor = self.branching_factor
        new_node1 = _CFNode(threshold, branching_factor,
                            is_leaf=child_node.is_leaf,
                            n_features=child_node.n_features)
        new_node2 = _CFNode(threshold, branching_factor,
                            is_leaf=child_node.is_leaf,
                            n_features=child_node.n_features)
        newsubcluster1.child_ = new_node1
        newsubcluster2.child_ = new_node2

        if child_node.is_leaf:
            if child_node.prev_leaf_ is not None:
                child_node.prev_leaf_.next_leaf_ = new_node1
            new_node1.prev_leaf_ = child_node.prev_leaf_
            new_node1.next_leaf_ = new_node2
            new_node2.prev_leaf_ = new_node1
            new_node2.next_leaf_ = child_node.next_leaf_
            if child_node.next_leaf_ is not None:
                child_node.next_leaf_.prev_leaf_ = new_node2

        for idx, subcluster in enumerate(child_node.subclusters_):
            if dist_idx[0][idx] > dist_idx[1][idx]:
                new_node1.update(subcluster)
                newsubcluster1.update(subcluster)
            else:
                new_node2.update(subcluster)
                newsubcluster2.update(subcluster)

        ind = self.subclusters_.index(parent_subcluster)
        self.subclusters_.remove(parent_subcluster)
        self.init_centroids_[ind: self.n_ - 1, :] = \
            self.init_centroids_[ind + 1: self.n_, :]
        self.init_sq_norm_[ind: self.n_ - 1] = \
            self.init_sq_norm_[ind + 1: self.n_]
        self.n_ -= 1
        self.update(newsubcluster1)
        self.update(newsubcluster2)

    def insert_cf_subcluster(self, subcluster):
        """
        Insert a new subcluster into the node
        """
        if not self.subclusters_:
            self.update(subcluster)
            return False

        # We need to find the closest subcluster among all the
        # subclusters so that we can insert our new subcluster.
        dist_matrix = safe_sparse_dot(
            self.centroids_, subcluster.centroid_)
        dist_matrix *= -2.
        dist_matrix += self.squared_norm_
        closest_index = np.argmin(dist_matrix)
        closest_subcluster = self.subclusters_[closest_index]

        # If the subcluster has a child, we need a recursive strategy.
        if closest_subcluster.child_ is not None:
            split_child = closest_subcluster.child_.insert_cf_subcluster(
                subcluster)

            if not split_child:
                # If it is determined that the child need not be split, we
                # can just update the closest_subcluster
                closest_subcluster.update(subcluster)
                self.init_centroids_[closest_index] = \
                    self.subclusters_[closest_index].centroid_
                self.init_sq_norm_[closest_index] = \
                    self.subclusters_[closest_index].sq_norm_
                return False

            # things not too good. we need to redistribute the subclusters in
            # our child node, and add a new subcluster in the parent
            # subcluster to accomodate the new child.
            else:
                self.split_node(closest_subcluster.child_, closest_subcluster)

                if len(self.subclusters_) > self.branching_factor:
                    return True
                return False

        # good to go!
        else:
            should_merge = closest_subcluster.merge_subcluster(
                subcluster, self.threshold)
            if should_merge:
                self.init_centroids_[closest_index] = \
                    closest_subcluster.centroid_
                self.init_sq_norm_[closest_index] = \
                    closest_subcluster.sq_norm_
                return False

            # not close to any other subclusters, and we still
            # have space, so add.
            elif len(self.subclusters_) < self.branching_factor:
                self.update(subcluster)
                return False

            # We do not have enough space nor is it closer to an
            # other subcluster. We need to split.
            else:
                self.update(subcluster)
                return True


class _CFSubcluster(object):
    """Each subcluster in a CFNode is called a CFSubcluster.

    A CFSubcluster can have a CFNode has its child.

    Parameters
    ----------
    X : ndarray, optional
        Sample. This is kept optional to allow initialization of empty
        subclusters.

    index : int, optional
        Index of the array in the original data. This enables to
        retrieve the final subclusters.

    Attributes
    ----------
    n_ : int
        Number of samples that belong to each subcluster.

    ls_ : ndarray
        Linear sum of all the samples in a subcluster. Prevents holding
        all sample data in memory.

    ss_ : float
        Sum of the squared l2 norms of all samples belonging to a subcluster.

    centroid_ : ndarray
        Centroid of the subcluster. Prevent recomputing of centroids when
        self.centroids_ is called.

    child_ : _CFNode
        Child Node of the subcluster. Once a given _CFNode is set as the child
        of the _CFNode, it is set to self.child_.

    sq_norm_ : ndarray
        Squared norm of the subcluster. Used to prevent recomputing when
        pairwise minimum distances are computed.
    """
    def __init__(self, X=None):
        self.X = X
        if X is None:
            self.n_ = 0
            self.ls_ = self.ss_ = 0.0
        else:
            self.n_ = 1
            self.ls_ = self.centroid_ = X
            self.ss_ = np.dot(self.ls_, self.ls_)
            self.sq_norm_ = np.dot(self.ls_, self.ls_)
        self.child_ = None

    def update(self, subcluster):
        self.n_ += subcluster.n_
        self.ls_ = self.ls_ + subcluster.ls_
        self.ss_ = self.ss_ + subcluster.ss_
        self.centroid_ = self.ls_ / self.n_
        self.sq_norm_ = np.dot(self.centroid_, self.centroid_)

    def merge_subcluster(self, nominee_cluster, threshold):
        """Check if a cluster is worthy enough to be merged. If
           yes than merge."""
        new_ss = self.ss_ + nominee_cluster.ss_
        new_ls = self.ls_ + nominee_cluster.ls_
        new_n = self.n_ + nominee_cluster.n_
        new_centroid = new_ls / new_n
        dot_product = -2 * np.dot(new_ls, new_centroid)
        new_norm = np.dot(new_centroid, new_centroid)
        new_radius = (new_ss + dot_product) / new_n + new_norm
        if new_radius <= threshold ** 2:
            self.n_, self.ls_, self.ss_, self.centroid_, self.sq_norm_ =  (
                new_n, new_ls, new_ss, new_centroid, new_norm)
            return True
        return False

    def radius(self):
        """Return radius of the subcluster"""
        dot_product = -2 * np.dot(self.ls_, self.centroid_)
        return sqrt(((self.ss_ + dot_product) / self.n_) + self.sq_norm_)


class Birch(BaseEstimator, TransformerMixin, ClusterMixin):
    """Implements the Birch clustering algorithm.

    Every new sample is inserted into the root of the Clustering Feature
    Tree. It is then clubbed together with the subcluster that has the
    centroid closest to the new sample. This is done recursively till it
    ends up at the subcluster of the leaf of the tree has the closest centroid.

    Parameters
    ----------
    threshold : float, default 0.5
        The radius of the subcluster obtained by merging a new sample and the
        closest subcluster should be greater than the square of the threshold.
        If the radius is greater than the square of the threshold, than a new
        subcluster is started.

    branching_factor : int, default 50
        Maximum number of CF subclusters in each node. If a new samples enters
        such that the number of subclusters exceed the branching_factor then
        the node has to be split. The corresponding parent also has to be
        split and if the number of subclusters in the parent is greater than
        the branching factor, then it has to be split recursively.

    global_clusters : int, instance of sklearn.cluster model, default None
        Number of clusters after the final clustering step, which treats the
        subclusters from the leaves as new samples. By default the global
        clustering step is AgglomerativeClustering with global_clusters set to 3.
        By default, this final clustering step is not performed and the
        subclusters are returned as they are.

    compute_labels : bool, default True
        Whether or not to compute labels for each fit.

    Attributes
    ----------
    root_ : _CFNode
        Root of the CFTree.

    dummy_leaf_ : _CFNode
        Start pointer to all the leaves.

    subcluster_centers_ : ndarray,
        Centroids of all subclusters read directly from the leaves.

    subcluster_labels_ : ndarray,
        Labels assigned to the centroids of the subclusters after
        they are clustered globally.

    labels_ : ndarray, shape (n_samples,)
        Array of labels assigned to the input data.
        if partial_fit is used instead of fit, they are assigned to the
        last batch of data.

    Examples
    --------
    >>> from sklearn.cluster import Birch
    >>> X = [[0, 1], [0.3, 1], [-0.3, 1], [0, -1], [0.3, -1], [-0.3, -1]]
    >>> brc = Birch(branching_factor=50, global_clusters=None, threshold=0.5,
    ... compute_labels=True)
    >>> brc.fit(X)
    Birch(branching_factor=50, compute_labels=True, global_clusters=None,
       threshold=0.5)
    >>> brc.predict(X)
    array([0, 0, 0, 1, 1, 1])

    References
    ----------
    * Tian Zhang, Raghu Ramakrishnan, Maron Livny
      BIRCH: An efficient data clustering method for large databases.
      http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf

    * Roberto Perdisci
      JBirch - Java implementation of BIRCH clustering algorithm
      https://code.google.com/p/jbirch/
    """

    def __init__(self, threshold=0.5, branching_factor=50, global_clusters=3,
                 compute_labels=True):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.global_clusters = global_clusters
        self.compute_labels = compute_labels
        self.partial_fit_ = False
        self.root_ = None

    def fit(self, X, y=None):
        """
        Build a CF Tree for the input data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data.
        """
        X = check_array(X, accept_sparse='csr')
        threshold = self.threshold
        branching_factor = self.branching_factor

        if branching_factor <= 1:
            raise ValueError("Branching_factor should be greater than one.")
        n_samples, n_features = X.shape
        # If partial_fit is called for the first time or if fit is called,
        # we need to initialize the root.
        if not self.partial_fit_ or (self.root_ is None and self.partial_fit_):
            # The first root is the leaf. Manipulate this object throughout.
            self.root_ = _CFNode(threshold, branching_factor, is_leaf=True,
                                 n_features=n_features)

            # To enable getting back subclusters.
            self.dummy_leaf_ = _CFNode(threshold, branching_factor,
                                       is_leaf=True, n_features=n_features)
            self.dummy_leaf_.next_leaf_ = self.root_

        # Cannot vectorize. Enough to convince to use cython.
        if not sparse.issparse(X):
            iter_func = iter
        else:
            iter_func = _iterate_sparse_X

        for sample in iter_func(X):
            subcluster = _CFSubcluster(sample)
            split = self.root_.insert_cf_subcluster(subcluster)

            if split:
                subclusters_pairwise = euclidean_distances(
                    self.root_.centroids_,
                    Y_norm_squared=self.root_.squared_norm_,
                    squared=True)

                # Separating the two farthest clusters.
                n_clusters = subclusters_pairwise.shape[0]
                farthest_idx = np.unravel_index(
                    subclusters_pairwise.argmax(),
                    (n_clusters, n_clusters))
                dist_idx = subclusters_pairwise[[farthest_idx]]

                new_node1 = _CFNode(
                    threshold, branching_factor, is_leaf=self.root_.is_leaf,
                    n_features=n_features)
                new_node2 = _CFNode(
                    threshold, branching_factor, is_leaf=self.root_.is_leaf,
                    n_features=n_features)

                new_subcluster1 = _CFSubcluster()
                new_subcluster2 = _CFSubcluster()
                new_subcluster1.child_ = new_node1
                new_subcluster2.child_ = new_node2

                for idx, subcluster in enumerate(self.root_.subclusters_):
                    if dist_idx[0][idx] > dist_idx[1][idx]:
                        new_node1.update(subcluster)
                        new_subcluster1.update(subcluster)
                    else:
                        new_node2.update(subcluster)
                        new_subcluster2.update(subcluster)

                if self.root_.is_leaf:
                    self.dummy_leaf_.next_leaf_ = new_node1
                    new_node1.prev_leaf_ = self.dummy_leaf_
                    new_node1.next_leaf_ = new_node2
                    new_node2.prev_leaf_ = new_node1

                del self.root_
                self.root_ = _CFNode(threshold, branching_factor,
                                     is_leaf=False,
                                     n_features=n_features)
                self.root_.update(new_subcluster1)
                self.root_.update(new_subcluster2)

        centroids = np.concatenate([
            leaf.centroids_ for leaf in self.get_leaves()])
        self.subcluster_centers_ = centroids

        self._global_clustering(X)
        return self

    def get_leaves(self):
        """
        Retrieve the leaves of the CF Node.

        Returns
        -------
        leaves: array-like
            List of the leaf nodes.
        """
        leaf_ptr = self.dummy_leaf_.next_leaf_
        leaves = []
        while leaf_ptr:
            leaves.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf_
        return leaves

    def _check_fit(self, X):
        if not hasattr(self, 'subcluster_centers_'):
            raise ValueError("Fit training data before predicting")
        if X.shape[1] != self.subcluster_centers_.shape[1]:
            raise ValueError(
                "Training data and predicted data do "
                "not have same no. of features.")

    def predict(self, X):
        """
        Predict data using the centroids_ of subclusters.

        Avoid computation of the row norms of X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data.

        Returns
        -------
        labels: ndarray, shape(n_samples)
            Labelled data.
        """
        X = check_array(X, accept_sparse='csr')
        self._check_fit(X)
        reduced_distance = safe_sparse_dot(X, self.subcluster_centers_.T)
        reduced_distance *= -2
        reduced_distance += self._subcluster_norms
        return self.subcluster_labels_[np.argmin(reduced_distance, axis=1)]

    def partial_fit(self, X=None, y=None):
        """
        Online learning. Prevents rebuilding of CFTree from scratch.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features), None
            Input data. If X is not provided, only the global clustering
            step is done.
        """
        self.partial_fit_ = True
        if X is None:
            # Perform just the final global clustering step.
            self._global_clustering()
            return self
        else:
            return self.fit(X)

    def transform(self, X, y=None):
        """
        Transform X into subcluster centroids dimension.

        Each dimension represents the distance between each cluster centroid.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data.

        Returns
        -------
        X_trans: {array-like, sparse matrix}, shape (n_samples, n_clusters)
            Transformed data.
        """
        if not hasattr(self, 'subcluster_centers_'):
            raise ValueError("Fit training data before predicting")
        return euclidean_distances(X, self.subcluster_centers_)

    def _global_clustering(self, X=None):
        """
        Global clustering for the subclusters obtained after fitting
        """
        clusters = self.global_clusters
        centroids = self.subcluster_centers_
        compute_labels = (X is not None) and self.compute_labels

        # Preprocessing for the global clustering.
        not_enough_centroids = False
        if hasattr(clusters, 'fit_predict'):
            global_cluster = clusters
        elif isinstance(clusters, int):
            global_cluster = AgglomerativeClustering(
                n_clusters=clusters)
            # There is no need to perform the global clustering step.
            if len(centroids) < clusters:
                not_enough_centroids = True
        elif clusters is not None:
            raise ValueError("n_clusters should be an instance of "
                             "ClusterMixin or an int")

        # To use in predict to avoid recalculation.
        if compute_labels:
            self._subcluster_norms = row_norms(
                self.subcluster_centers_, squared=True)

        if self.global_clusters is None or not_enough_centroids:
            self.subcluster_labels_ = np.arange(len(centroids))
            if not_enough_centroids:
                warnings.warn(
                    "Number of subclusters found (%d) by Birch is less "
                    "than (%d). Decrease the threshold."
                    % (len(centroids), self.global_clusters))
        else:
            # The global clustering step that clusters the subclusters of
            # the leaves. It assumes the centroids of the subclusters as
            # samples and finds the final centroids.
            self.subcluster_labels_ = global_cluster.fit_predict(
                self.subcluster_centers_)

        if compute_labels:
            self.labels_ = self.predict(X)
