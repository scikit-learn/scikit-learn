"""Markov clustering algorithm."""

# Author: Uri Goren, uri@goren4u.com

from scipy.sparse import eye, csr_matrix
from ..preprocessing import normalize
from ..metrics.pairwise import pairwise_distances
from collections import defaultdict


def hash_list_of_ints(ints):
    return ",".join(map(str, sorted(ints)))


class MarkovClustering:
    def __init__(self, metric="cosine", bias=1, inflation_power=2,
                 inflation_threshold=1e-10, self_loops_weight=0.01,
                 expansion_power=2, iteration_limit=100):
        """
        Initializing similarity matrix
        Either by setting it (metric = None)
        Or by deriving it from a distance function (similarity = bias-distance)
        """
        self.labels_ = None
        self.metric = metric
        self.bias = bias
        self.inflation_power = inflation_power
        self.inflation_threshold = inflation_threshold
        self.self_loops_weight = self_loops_weight
        self.expansion_power = expansion_power
        self.iteration_limit = iteration_limit

    def normalize(self, T):
        return normalize(T, norm='l1', axis=1)

    def self_loops(self, T):
        return eye(T.shape[0]) * self.self_loops_weight + T

    def expansion(self, T):
        ret = T
        for _ in range(1, self.expansion_power):
            ret = ret * T
        return ret

    def inflation(self, T):
        for i in range(len(T.data)):
            if T.data[i] < self.inflation_threshold:
                T.data[i] = 0
            else:
                T.data[i] = T.data[i] ** self.inflation_power
        return T

    def fit(self, X, verbose=False):
        if self.metric is None:
            T = X
        else:
            T = csr_matrix(self.bias-pairwise_distances(X, metric=self.metric))
        iterations = 0
        prev_T = csr_matrix(T.shape)
        T = self.self_loops(T)
        while (iterations < self.iteration_limit) and ((prev_T - T).nnz != 0):
            prev_T = T.copy()
            iterations += 1
            T = self.normalize(T)
            T = self.expansion(T)
            T = self.inflation(T)
            if verbose:
                print("========Iteration #{i}=======".format(i=iterations))
                print(T.toarray())
        self.labels_ = self.extract_labels(T)
        return self

    def extract_labels(self, T):
        M = T.tocoo()
        rows = defaultdict(set)
        for i, d in enumerate(M.data):
            if d == 0:
                continue
            rows[M.row[i]].add(M.col[i])
        row_hashes = [hash_list_of_ints(rows[i]) for i in range(M.shape[0])]
        d = dict([(l, i) for i, l in enumerate(set(row_hashes))])
        labels = [d[row_hashes[i]] for i in range(M.shape[0])]
        return labels

    def clusters(self, labels=None):
        ret = defaultdict(set)
        for i, c in enumerate(self.labels_):
            if labels is None:
                ret[c].add(i)
            else:
                ret[c].add(labels[i])
        return ret

if __name__ == "__main__":
    A = [[1, 2, 0], [2, 1, 0], [0, 0, 4], [0, 0, 9]]
    clusters = MarkovClustering().fit(A, verbose=True).clusters()
    print(clusters)
