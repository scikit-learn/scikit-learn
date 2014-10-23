"""
A Proof of Concept for Birch according to http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf
The optional requirements are not handled, i.e
1. Rebuilding of tree if its higher than the required memory.
2. Merging refinement.

"""
import numpy as np
from math import sqrt

from sklearn.metrics.pairwise import pairwise_distances, pairwise_distances_argmin_min
from sklearn.base import TransformerMixin


class CFNode(object):
    r"""
    Each node in a CFTree is called a CFNode.

    The CFNode has a branching_factor number of CFSubclusters.

    Parameters
    ==========
    threshold : float
        Threshold needed for a new subcluster to enter a CFSubcluster.

    branching_factor : int
        Number of CF subclusters in each node.

    is_leaf : bool
        We need to know if the CFNode is a leaf or not, in order to
        retrieve the final subclusters.
    """

    def __init__(self, threshold, branching_factor, is_leaf):
        self.threshold = threshold
        self.branching_factor = branching_factor
        self.is_leaf = is_leaf

        # The list of subclusters to manipulate throughout
        self.subclusters = []
        self.prev_leaf = None
        self.next_leaf = None

    def set_prev_leaf(self, node):
        self.prev_leaf = node

    def set_next_leaf(self, node):
        self.next_leaf = node

    def splitNode(self, child_node, parent_subcluster):
        r"""
        Split the node if there is no place for a new subcluster.
        """
        subclusters_pairwise = pairwise_distances([subcluster.ls / subcluster.n
            for subcluster in child_node.subclusters])
        max_ = subclusters_pairwise.argmax()
        farthest_idx1 = max_ / subclusters_pairwise.shape[0]
        farthest_idx2 = max_ % subclusters_pairwise.shape[0]

        dist_idx = subclusters_pairwise[np.array([farthest_idx1, farthest_idx2])]

        newsubcluster1 = CFSubcluster()
        newsubcluster2 = CFSubcluster()

        newNode1 = CFNode(self.threshold,
            self.branching_factor, child_node.is_leaf)
        newNode2 = CFNode(self.threshold,
            self.branching_factor, child_node.is_leaf)
        newsubcluster1.child = newNode1
        newsubcluster2.child = newNode2

        if child_node.is_leaf:
            if child_node.prev_leaf is not None:
                child_node.prev_leaf.set_next_leaf(newNode1)
            newNode1.set_prev_leaf(child_node.prev_leaf)
            newNode1.set_next_leaf(newNode2)
            newNode2.set_prev_leaf(newNode1)
            newNode2.set_next_leaf(child_node.next_leaf)
            if child_node.next_leaf is not None:
                child_node.next_leaf.set_prev_leaf(newNode2)

        for idx, subcluster in enumerate(child_node.subclusters):
            if dist_idx[0][idx] > dist_idx[1][idx]:
                newNode1.subclusters.append(subcluster)
                newsubcluster1.updatesubcluster(subcluster)
            else:
                newNode2.subclusters.append(subcluster)
                newsubcluster2.updatesubcluster(subcluster)

        self.subclusters.remove(parent_subcluster)
        self.subclusters.append(newsubcluster1)
        self.subclusters.append(newsubcluster2)


    def insert_cf_subcluster(self, subcluster):
        if not self.subclusters:
            self.subclusters.append(subcluster)
            return False

        # We need to find the closest subcluster among all the
        # subclusters so that we can insert our new subcluster.

        subcluster_centroids = [(old_subcluster.ls / old_subcluster.n) for old_subcluster in self.subclusters]
        closest_index, closest_threshold = pairwise_distances_argmin_min(
            subcluster.ls, subcluster_centroids)
        closest_subcluster = self.subclusters[closest_index]

        # If the subcluster has a child, we need a recursive strategy.
        if closest_subcluster.child is not None:
            split_child = closest_subcluster.child.insert_cf_subcluster(subcluster)

            if not split_child:
                # If it is determined that the child need not be split, we
                # can just update the closest_subcluster
                self.subclusters[closest_index].updatesubcluster(subcluster)
                return False

            # things not too good. we need to redistribute the subclusters in our
            # child node, and add a new subcluster in the parent subcluster to
            # accomodate the new child. 
            else:
                self.splitNode(closest_subcluster.child, closest_subcluster)

                if len(self.subclusters) > self.branching_factor:
                    return True
                return False

        # good to go!
        elif self.threshold >= closest_threshold:
            closest_subcluster.updatesubcluster(subcluster)
            return False

        # not close to any other subclusters, and we still have space, so add.
        elif len(self.subclusters) < self.branching_factor:
            self.subclusters.append(subcluster)
            return False

        # We do not have enough space nor is it closer to an other subcluster.
        # We need to split.
        else:
            self.subclusters.append(subcluster)
            return True


class CFSubcluster(object):
    r"""
    Each subcluster(subcluster) in a CFNode is called a CFSubcluster.

    A CFSubcluster has a CFNode has its child.

    Parameters
    ==========
    X : ndarray, optional
        Sample.

    index : int, optional
        Integer of the array in the original data.

    """
    def __init__(self, X=None, index=None):
        self.X = X
        self.sample_indices = []
        if index is not None:
            self.sample_indices = [index]
        if X is None:
            self.n = 0
            self.ls = 0.0
            self.ss = 0.0
        else:
            self.n = 1
            self.ls = X
            self.ss = X**2
        self.child = None

    def setChild(self, node):
        self.child = node

    def updatesubcluster(self, subcluster):
        self.sample_indices.extend(subcluster.sample_indices)
        self.n += subcluster.n
        self.ls += subcluster.ls
        self.ss += subcluster.ss


class Birch(TransformerMixin):
    r"""
    Tries to implement the Birch algorithm.

    Parameters
    ==========
    threshold : float
        Threshold needed for a new subcluster to enter a CFSubcluster.

    branching_factor : int
        Number of CF subclusters in each node.
    """
    def __init__(self, threshold=1.0, branching_factor=8):
        self.threshold = threshold
        self.branching_factor = branching_factor

    def fit(self, X):

        # The first root is the leaf. Manipulate this object throughout.
        self.root_ = CFNode(self.threshold, self.branching_factor, True)

        # To enable getting back subclusters.
        self.dummy_leaf = CFNode(self.threshold, self.branching_factor, True)
        self.dummy_leaf.set_next_leaf(self.root_)
        # Cannot vectorize. Enough to convince to use cython.
        for ind, sample in enumerate(X):
            subcluster = CFSubcluster(sample, ind)

            split = self.root_.insert_cf_subcluster(subcluster)

            if split:
                subclusters_pairwise = pairwise_distances([
                    subcluster.ls / subcluster.n for subcluster in self.root_.subclusters]
                    )
                max_ = subclusters_pairwise.argmax()
                farthest_idx1 = max_ / subclusters_pairwise.shape[0]
                farthest_idx2 = max_ % subclusters_pairwise.shape[0]
                dist_idx = subclusters_pairwise[np.array([
                    farthest_idx1, farthest_idx2])]

                newsubcluster1 = CFSubcluster()
                newsubcluster2 = CFSubcluster()

                # The New nodes are leaves.
                newNode1 = CFNode(self.threshold, self.branching_factor, self.root_.is_leaf)
                newNode2 = CFNode(self.threshold, self.branching_factor, self.root_.is_leaf)
                newsubcluster1.child = newNode1
                newsubcluster2.child = newNode2

                for idx, subcluster in enumerate(self.root_.subclusters):
                    if dist_idx[0][idx] > dist_idx[1][idx]:
                        newNode1.subclusters.append(subcluster)
                        newsubcluster1.updatesubcluster(subcluster)
                    else:
                        newNode2.subclusters.append(subcluster)
                        newsubcluster2.updatesubcluster(subcluster)

                if self.root_.is_leaf:
                    self.dummy_leaf.set_next_leaf(newNode1)
                    newNode1.set_prev_leaf(self.dummy_leaf)
                    newNode1.set_next_leaf(newNode2)
                    newNode2.set_prev_leaf(newNode1)

                del self.root_
                self.root_ = CFNode(self.threshold, self.branching_factor, False)
                self.root_.subclusters.append(newsubcluster1)
                self.root_.subclusters.append(newsubcluster2)


    def getSubclusters(self):
        r"""
        Retrieve the subclusters.
        """
        leaf_ptr = self.dummy_leaf.next_leaf
        subclusters = []
        while leaf_ptr:
            subclusters.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf
        return subclusters
