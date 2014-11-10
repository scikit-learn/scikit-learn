# Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause
from __future__ import division

import warnings
import numpy as np
from scipy import sparse

from ..metrics.pairwise import euclidean_distances
from ..base import TransformerMixin, ClusterMixin
from ..externals.six.moves import xrange
from ..utils import check_array
from ..utils.extmath import safe_sparse_dot
from .hierarchical import AgglomerativeClustering

def iterate_X(X):
    n_samples = X.shape[0]
    X_issparse = sparse.issparse(X)
    if X_issparse:
        X_indices = X.indices
        X_data = X.data
        X_indptr = X.indptr

    for i in xrange(n_samples):
        if X_issparse:
            row = np.zeros(X.shape[1])
            startptr, endptr = X_indptr[i], X_indptr[i + 1]
            nonzero_indices = X_indices[startptr: endptr]
            row[nonzero_indices] = X_data[startptr: endptr]
            yield row
        else:
            yield X[i]


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

    centroids_ : array-like
        list of centroids for a particular CFNode.

    squared_norm_ : array-like
        list of squared norms for a particular CFNode.
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
                self.subclusters_[closest_index].update(subcluster)
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
            closest_threshold = subcluster.sq_norm_  + dist_matrix[closest_index]
            if self.threshold ** 2 >= closest_threshold:
                self.subclusters_[closest_index].update(subcluster)
                self.init_centroids_[closest_index] = \
                    self.subclusters_[closest_index].centroid_
                self.init_sq_norm_[closest_index] = \
                    self.subclusters_[closest_index].sq_norm_
                return False

            # not close to any other subclusters, and we still have space, so add.
            elif len(self.subclusters_) < self.branching_factor:
                self.update(subcluster)
                return False

            # We do not have enough space nor is it closer to an other subcluster.
            # We need to split.
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
            self.ls_ = 0.0
        else:
            self.n_ = 1
            self.ls_ = self.centroid_ = X
            self.sq_norm_ = np.dot(self.ls_, self.ls_)
        self.child_ = None

    def update(self, subcluster):
        self.n_ += subcluster.n_
        self.ls_ = self.ls_ + subcluster.ls_
        self.centroid_ = self.ls_ / self.n_
        self.sq_norm_ = np.dot(self.centroid_, self.centroid_)


class Birch(TransformerMixin, ClusterMixin):
    """Implements the Birch clustering algorithm.

    Every new sample is inserted into the root of the Clustering Feature
    Tree. It is then clubbed together with the subcluster that has the
    centroid closest to the new sample. This is done recursively till it
    ends up at the subcluster of the leaf of the tree has the closest centroid.

    Parameters
    ----------
    threshold : float, default 1.0
        The minimum distance between a new sample and the closest subcluster
        centroid for it to be part of the subcluster. If the distance
        is greater than this threshold, than a new subcluster is started.

    branching_factor : int, default 8
        Maximun number of CF subclusters in each node. If a new samples enters
        such that the number of subclusters exceed the branching_factor then
        the node has to be split. The corresponding parent also has to be
        split and if the number of subclusters in the parent is greater than
        the branching factor, then it has to be split recursively.

    n_clusters : int, instance of ClusterMixin or None, default 3
        Number of clusters after the final clustring step, which treats the
        subclusters from the leaves as new samples. By default the global
        clustering step is AgglomerativeClustering with n_clusters set to 3.
        If set to None, this final clustering step is not performed and the
        subclusters are returned as they are.

    Attributes
    ----------
    root_ : _CFNode
        Root of the CFTree.

    dummy_leaf_ : CFNode
        Start pointer to all the leaves.

    centroids_ : ndarray
        Centroids of all subclusters read directly from the leaves.

    labels_ : ndarray
        Array of labels assigned to the input data.

    References
    ----------
    * Tian Zhang, Raghu Ramakrishnan, Maron Livny
      BIRCH: An efficient data clustering method for large databases.
      http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf

    * Roberto Perdisci
      JBirch - Java implementation of BIRCH clustering algorithm
      https://code.google.com/p/jbirch/
    """

    def __init__(self, threshold=1.0, branching_factor=8, n_clusters=3):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.n_clusters = n_clusters
        self.partial_fit_ = False
        self.root_ = None

    def fit(self, X):
        """
        Build a CF Tree for the input data.

        Parameters
        ----------
        X : ndarray
            Input data.
        """
        X = check_array(X, accept_sparse='csr')
        threshold = self.threshold
        branching_factor = self.branching_factor

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
        for sample in iterate_X(X):
            subcluster = _CFSubcluster(sample)
            split = self.root_.insert_cf_subcluster(subcluster)

            if split:
                subclusters_pairwise = euclidean_distances(
                    self.root_.centroids_,
                    Y_norm_squared=self.root_.squared_norm_,
                    squared=True)

                # Separating the two farthest clusters.
                n_clusters = subclusters_pairwise.shape[0]
                farthest_idx = np.unravel_index(subclusters_pairwise.argmax(),
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

        # To get labels_ and centroids_
        leaves = self.get_leaves()
        centroids = list()
        weights = list()
        for leaf in leaves:
            centroids.append(leaf.centroids_)
            for sf in leaf.subclusters_:
                weights.append(sf.n_)

        # Weights assigned to each centroid for the global clustering step.
        centroids = np.concatenate(centroids)
        weights = np.asarray(weights)

        # Preprocessing for the global clustering.
        not_enough_centroids = False
        if isinstance(self.n_clusters, ClusterMixin):
            clf = self.n_clusters
        elif isinstance(self.n_clusters, int):
            clf = AgglomerativeClustering(n_clusters=self.n_clusters)
            # There is no need to perform the global clustering step.
            if len(centroids) < self.n_clusters:
                not_enough_centroids = True
        elif self.n_clusters is not None:
            raise ValueError("n_clusters should be an instance of "
                             "ClusterMixin or an int")            

        if self.n_clusters is None or not_enough_centroids:
            self.centroids_ = centroids
            self.labels_ = self.predict(X)
            if not_enough_centroids:
                warnings.warn(
                    "Number of clusters %s found by Birch is lesser than "
                    "%s. Decrease the threshold."
                    % (len(centroids), self.n_clusters))
            return self

        # The global clustering step that clusters the subclusters of
        # the leaves. It assumes the centroids of the subclusters as
        # samples and finds the final centroids.
        labels = clf.fit_predict(centroids)
        n_clusters = len(np.unique(labels))
        new_centroids = np.empty((n_clusters, X.shape[1]))
        for i in xrange(n_clusters):
            mask = labels == i
            new_centroids[i] = np.average(centroids[mask], axis=0,
                                          weights=weights[mask])
        self.centroids_ = new_centroids
        self.labels_ = self.predict(X)
        return self

    def get_leaves(self):
        """
        Retrieve the leaves.
        """
        leaf_ptr = self.dummy_leaf_.next_leaf_
        leaves = []
        while leaf_ptr:
            leaves.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf_
        return leaves

    def predict(self, X):
        """
        Predict data using the centroids_ of subclusters.
 
        Parameters
        ----------
        X : ndarray
            Input data.
        """
        X_dot_centroid = safe_sparse_dot(X, self.centroids_.T)
        X_dot_centroid *= -2
        X_dot_centroid += np.sum(self.centroids_ ** 2, axis=1)
        return np.argmin(X_dot_centroid, axis=1)

    def partial_fit(self, X):
        """
        Online learning. Prevents rebuilding of CFTree from scratch.

        Parameters
        ----------
        X : ndarray
            Input data.
        """
        self.partial_fit_ = True
        return self.fit(X)
