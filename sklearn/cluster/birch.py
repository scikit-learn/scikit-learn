"""
A Proof of Concept for Birch according to http://www.cs.sfu.ca/CourseCentral/459/han/papers/zhang96.pdf
The optional requirements are not handled, i.e
1. Rebuilding of tree if its higher than the required memory.
2. Merging refinement.

"""

# XXX : add license + author names

import numpy as np

from sklearn.metrics.pairwise import pairwise_distances, pairwise_distances_argmin_min
from sklearn.base import TransformerMixin


class CFNode(object):
    r"""Each node in a CFTree is called a CFNode.

    The CFNode has a branching_factor number of CFSubclusters.

    Parameters
    ----------
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

    def get_centroids(self):
        return [sc.ls / sc.n for sc in self.subclusters]

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
        newsubcluster1.child = new_node1
        newsubcluster2.child = new_node2

        if child_node.is_leaf:
            if child_node.prev_leaf is not None:
                child_node.prev_leaf.next_leaf = new_node1
            new_node1.prev_leaf = child_node.prev_leaf
            new_node1.next_leaf = new_node2
            new_node2.prev_leaf = new_node1
            new_node2.next_leaf = child_node.next_leaf
            if child_node.next_leaf is not None:
                child_node.next_leaf.prev_leaf = new_node2

        for idx, subcluster in enumerate(child_node.subclusters):
            if dist_idx[0][idx] > dist_idx[1][idx]:
                new_node1.subclusters.append(subcluster)
                newsubcluster1.update(subcluster)
            else:
                new_node2.subclusters.append(subcluster)
                newsubcluster2.update(subcluster)

        self.subclusters.remove(parent_subcluster)
        self.subclusters.append(newsubcluster1)
        self.subclusters.append(newsubcluster2)

    def insert_cf_subcluster(self, subcluster):
        """XXX
        """
        if not self.subclusters:
            self.subclusters.append(subcluster)
            return False

        # We need to find the closest subcluster among all the
        # subclusters so that we can insert our new subcluster.

        subcluster_centroids = self.get_centroids()
        closest_index, closest_threshold = \
            pairwise_distances_argmin_min(subcluster.ls, subcluster_centroids)
        closest_subcluster = self.subclusters[closest_index]

        # If the subcluster has a child, we need a recursive strategy.
        if closest_subcluster.child is not None:
            split_child = closest_subcluster.child.insert_cf_subcluster(subcluster)

            if not split_child:
                # If it is determined that the child need not be split, we
                # can just update the closest_subcluster
                self.subclusters[closest_index].update(subcluster)
                return False

            # things not too good. we need to redistribute the subclusters in our
            # child node, and add a new subcluster in the parent subcluster to
            # accomodate the new child.
            else:
                self.split_node(closest_subcluster.child, closest_subcluster)

                if len(self.subclusters) > self.branching_factor:
                    return True
                return False

        # good to go!
        elif self.threshold >= closest_threshold:
            closest_subcluster.update(subcluster)
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
    r"""Each subcluster(subcluster) in a CFNode is called a CFSubcluster.

    A CFSubcluster has a CFNode has its child.

    Parameters
    ----------
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
            self.ss = X ** 2
        self.child = None

    def update(self, subcluster):
        self.sample_indices += subcluster.sample_indices
        self.n += subcluster.n
        self.ls += subcluster.ls
        self.ss += subcluster.ss


class Birch(TransformerMixin):
    r"""Implements the Birch algorithm.

    Parameters
    ----------
    threshold : float
        Threshold needed for a new subcluster to enter a CFSubcluster.

    branching_factor : int
        Number of CF subclusters in each node.

    Attributes
    ----------
    XXX
    """
    def __init__(self, threshold=1.0, branching_factor=8):
        self.threshold = threshold
        self.branching_factor = branching_factor

    def fit(self, X):
        threshold = self.threshold
        branching_factor = self.branching_factor

        # The first root is the leaf. Manipulate this object throughout.
        self.root_ = CFNode(threshold, branching_factor, True)

        # To enable getting back subclusters.
        self.dummy_leaf = CFNode(threshold, branching_factor, True)
        self.dummy_leaf.next_leaf = self.root_
        # Cannot vectorize. Enough to convince to use cython.
        for ind, sample in enumerate(X):
            subcluster = CFSubcluster(sample, ind)

            split = self.root_.insert_cf_subcluster(subcluster)

            if split:
                subclusters_pairwise = pairwise_distances(self.root_.get_centroids())
                max_ = subclusters_pairwise.argmax()
                farthest_idx1 = max_ / subclusters_pairwise.shape[0]
                farthest_idx2 = max_ % subclusters_pairwise.shape[0]
                dist_idx = subclusters_pairwise[[farthest_idx1, farthest_idx2]]

                new_node1 = CFNode(threshold, branching_factor, self.root_.is_leaf)
                new_node2 = CFNode(threshold, branching_factor, self.root_.is_leaf)

                new_subcluster1 = CFSubcluster()
                new_subcluster2 = CFSubcluster()
                new_subcluster1.child = new_node1
                new_subcluster2.child = new_node2

                for idx, subcluster in enumerate(self.root_.subclusters):
                    if dist_idx[0][idx] > dist_idx[1][idx]:
                        new_node1.subclusters.append(subcluster)
                        new_subcluster1.update(subcluster)
                    else:
                        new_node2.subclusters.append(subcluster)
                        new_subcluster2.update(subcluster)

                if self.root_.is_leaf:
                    self.dummy_leaf.next_leaf = new_node1
                    new_node1.prev_leaf = self.dummy_leaf
                    new_node1.next_leaf = new_node2
                    new_node2.prev_leaf = new_node1

                del self.root_
                self.root_ = CFNode(threshold, branching_factor, False)
                self.root_.subclusters.append(new_subcluster1)
                self.root_.subclusters.append(new_subcluster2)


            # Impedance matching with sklearn cluster API
            clusters = [sub.sample_indices for leaf in self.get_subclusters()
                        for sub in leaf.subclusters]
            centroids = [leaf.get_centroids() for leaf in self.get_subclusters()]
            centroids = np.concatenate(centroids)
            labels = np.empty(X.shape[0], dtype=np.int)
            for k, cluster in enumerate(clusters):
                labels[cluster] = k
            self.labels_ = labels
            self.centroids_ = centroids

    def get_subclusters(self):
        r"""
        Retrieve the subclusters.
        """
        leaf_ptr = self.dummy_leaf.next_leaf
        subclusters = []
        while leaf_ptr:
            subclusters.append(leaf_ptr)
            leaf_ptr = leaf_ptr.next_leaf
        return subclusters


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    from sklearn.datasets.samples_generator import make_blobs

    ##############################################################################
    # Generate sample data
    np.random.seed(0)

    batch_size = 45
    centers = [[1, 1], [-1, -1], [1, -1]]
    n_clusters = len(centers)
    X, labels_true = make_blobs(n_samples=3000, centers=centers, cluster_std=0.7)

    ##############################################################################
    # Compute clustering with Birch
    birch = Birch(threshold=1.0, branching_factor=8)
    birch.fit(X.copy())

    ##############################################################################
    # Plot result
    from itertools import cycle
    colors = cycle(['#4EACC5', '#FF9C34', '#4E9A06'])

    labels = birch.labels_
    centroids = birch.centroids_
    n_clusters = np.unique(labels).size
    print "n_clusters : %d" % n_clusters

    plt.figure()
    for this_centroid, k, col in zip(centroids, range(n_clusters), colors):
        mask = labels == k
        plt.plot(X[mask, 0], X[mask, 1], 'w',
                 markerfacecolor=col, marker='.')
        plt.plot(this_centroid[0], this_centroid[1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=10)
    plt.title('Birch')
    plt.show()
