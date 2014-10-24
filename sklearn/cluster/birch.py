# Authors: Manoj Kumar <manojkumarsivaraj334@gmail.com>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import numpy as np

from ..metrics.pairwise import (pairwise_distances,
                                pairwise_distances_argmin_min)
from ..base import TransformerMixin
from ..utils import check_array


class CFNode(object):
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

    prev_leaf_ : CFNode
        prev_leaf. Useful only if is_leaf is True.

    next_leaf_ : CFNode
        next_leaf. Useful only if is_leaf is True.
        the final subclusters.
    """
    def __init__(self, threshold, branching_factor, is_leaf):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.is_leaf = is_leaf

        # The list of subclusters to manipulate throughout
        self.subclusters_ = []
        self.prev_leaf_ = None
        self.next_leaf_ = None

    def get_centroids(self):
        return [sc.ls_ / sc.n_ for sc in self.subclusters_]

    def split_node(self, child_node, parent_subcluster):
        r"""
        Split the node if there is no place for a new subcluster.
        """
        dist = pairwise_distances(child_node.get_centroids())
        max_ = dist.argmax()
        farthest_idx1 = max_ / dist.shape[0]
        farthest_idx2 = max_ % dist.shape[0]

        dist_idx = dist[[farthest_idx1, farthest_idx2]]

        newsubcluster1 = CFSubcluster()
        newsubcluster2 = CFSubcluster()

        threshold = self.threshold
        branching_factor = self.branching_factor
        new_node1 = CFNode(threshold, branching_factor, child_node.is_leaf)
        new_node2 = CFNode(threshold, branching_factor, child_node.is_leaf)
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
                new_node1.subclusters_.append(subcluster)
                newsubcluster1.update(subcluster)
            else:
                new_node2.subclusters_.append(subcluster)
                newsubcluster2.update(subcluster)

        self.subclusters_.remove(parent_subcluster)
        self.subclusters_.append(newsubcluster1)
        self.subclusters_.append(newsubcluster2)

    def insert_cf_subcluster(self, subcluster):
        """
        Insert a new subcluster into the nide
        """
        if not self.subclusters_:
            self.subclusters_.append(subcluster)
            return False

        # We need to find the closest subcluster among all the
        # subclusters so that we can insert our new subcluster.

        subcluster_centroids = self.get_centroids()
        closest_index, closest_threshold = \
            pairwise_distances_argmin_min(subcluster.ls_, subcluster_centroids)

        # Index returned is a numpy array.
        closest_index = closest_index[0]
        closest_subcluster = self.subclusters_[closest_index]

        # If the subcluster has a child, we need a recursive strategy.
        if closest_subcluster.child_ is not None:
            split_child = closest_subcluster.child_.insert_cf_subcluster(
                subcluster)

            if not split_child:
                # If it is determined that the child need not be split, we
                # can just update the closest_subcluster
                self.subclusters_[closest_index].update(subcluster)
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
        elif self.threshold >= closest_threshold:
            closest_subcluster.update(subcluster)
            return False

        # not close to any other subclusters, and we still have space, so add.
        elif len(self.subclusters_) < self.branching_factor:
            self.subclusters_.append(subcluster)
            return False

        # We do not have enough space nor is it closer to an other subcluster.
        # We need to split.
        else:
            self.subclusters_.append(subcluster)
            return True


class CFSubcluster(object):
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
    sample_indices_ : array-like
        Array of sample indices that belong to each subcluster.

    n_ : int
        Number of samples that belong to each subcluster.

    ls_ : ndarray
        Linear sum of all the samples in a subcluster. Prevents holding
        all sample data in memory.

    ss_ : ndarray
        Squared sum of all the samples in a subcluster. Prevents holding
        all sample data in memory.

    child_ : CFNode
        Child Node of the subcluster.
    """
    def __init__(self, X=None, index=None):
        self.X = X
        self.sample_indices_ = []
        if index is not None:
            self.sample_indices_ = [index]
        if X is None:
            self.n_ = 0
            self.ls_ = 0.0
            self.ss_ = 0.0
        else:
            self.n_ = 1
            self.ls_ = X
            self.ss_ = X ** 2
        self.child_ = None

    def update(self, subcluster):
        self.sample_indices_ += subcluster.sample_indices_
        self.n_ += subcluster.n_
        self.ls_ += subcluster.ls_
        self.ss_ += subcluster.ss_


class Birch(TransformerMixin):
    """Implements the Birch algorithm.

    Insert a new sample is inserted into the CF Tree, the sample
    ultimately ends up at the leaf which has the closest centroid,
    and the parent sublusters are updated.

    TODO: Implement Merging Refinement
    TODO: Implement recomputing threshold and rebuilding the tree if
    the Memory Limit is breached.

    Parameters
    ----------
    threshold : float, default 1.0
        The minimum distance between a new sample and the closest subcluster
        centroid for it to be part of the subcluster. If the distance
        is greater than this threshold, than a new subcluster is started.

    branching_factor : int, default 8
        Maximun number of CF subclusters in each node. If a new samples enters
        such that the number of subclusters exceed the branching_factor then
        the node and the corresponding parent has to be updated.

    copy : bool, default True
        If set to False, X will be written inplace.

    Attributes
    ----------
    root_ : CFNode
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

    def __init__(self, threshold=1.0, branching_factor=8, copy=True):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.copy = copy

    def fit(self, X):
        """
        Build a CF Tree for the input data.

        Parameters
        ----------
        X : ndarray
            Input data.
        """
        X = check_array(X, copy=self.copy)
        threshold = self.threshold
        branching_factor = self.branching_factor

        # The first root is the leaf. Manipulate this object throughout.
        self.root_ = CFNode(threshold, branching_factor, True)

        # To enable getting back subclusters.
        self.dummy_leaf_ = CFNode(threshold, branching_factor, True)
        self.dummy_leaf_.next_leaf_ = self.root_

        # Cannot vectorize. Enough to convince to use cython.
        for ind, sample in enumerate(X):
            subcluster = CFSubcluster(sample, ind)

            split = self.root_.insert_cf_subcluster(subcluster)

            if split:
                subclusters_pairwise = pairwise_distances(
                    self.root_.get_centroids())
                max_ = subclusters_pairwise.argmax()
                farthest_idx1 = max_ / subclusters_pairwise.shape[0]
                farthest_idx2 = max_ % subclusters_pairwise.shape[0]
                dist_idx = subclusters_pairwise[[
                    farthest_idx1, farthest_idx2]]

                new_node1 = CFNode(
                    threshold, branching_factor, self.root_.is_leaf)
                new_node2 = CFNode(
                    threshold, branching_factor, self.root_.is_leaf)

                new_subcluster1 = CFSubcluster()
                new_subcluster2 = CFSubcluster()
                new_subcluster1.child_ = new_node1
                new_subcluster2.child_ = new_node2

                for idx, subcluster in enumerate(self.root_.subclusters_):
                    if dist_idx[0][idx] > dist_idx[1][idx]:
                        new_node1.subclusters_.append(subcluster)
                        new_subcluster1.update(subcluster)
                    else:
                        new_node2.subclusters_.append(subcluster)
                        new_subcluster2.update(subcluster)

                if self.root_.is_leaf:
                    self.dummy_leaf_.next_leaf_ = new_node1
                    new_node1.prev_leaf_ = self.dummy_leaf_
                    new_node1.next_leaf_ = new_node2
                    new_node2.prev_leaf_ = new_node1

                del self.root_
                self.root_ = CFNode(threshold, branching_factor, False)
                self.root_.subclusters_.append(new_subcluster1)
                self.root_.subclusters_.append(new_subcluster2)

        # Impedance matching with sklearn cluster API
        leaves = self.get_leaves()
        clusters = [sub.sample_indices_ for leaf in leaves
                    for sub in leaf.subclusters_]
        centroids = [leaf.get_centroids() for leaf in leaves]
        centroids = np.concatenate(centroids)
        labels = np.empty(X.shape[0], dtype=np.int)
        for k, cluster in enumerate(clusters):
            labels[cluster] = k
        self.labels_ = labels
        self.centroids_ = centroids
        return self

    def get_leaves(self):
        r"""
        Retrieve the leaves.
        """
        leaf_ptr = self.dummy_leaf_.next_leaf_
        leaves = []
        while leaf_ptr:
            leaves.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf_
        return leaves
