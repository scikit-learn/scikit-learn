import numpy as np

from sklearn.base import BaseEstimator
from sklearn.datasets import load_iris
from sklearn.feature_selection import SelectorMixin
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import check_random_state


class ReliefF(SelectorMixin, BaseEstimator):
    def __init__(
        self,
        k,
        prior="empirical",
        m="all",
        sigma=np.inf,
        neighbors_metric="manhattan",
        neighbors_algorithm="auto",
        metric_params=None,
        random_state=None,
    ):
        self.k = k
        self.prior = prior
        self.m = m
        self.sigma = sigma
        self.random_state = random_state
        self.neighbors_metric = neighbors_metric
        self.neighbors_algorithm = neighbors_algorithm
        self.metric_params = metric_params

    def fit(self, X, y):
        X, y = self._validate_data(
            X, y, accept_sparse=["csr", "csc"], multi_output=True
        )
        n_samples, n_features = X.shape[0], X.shape[1]
        random_state = check_random_state(self.random_state)
        self.weights_ = np.zeros(n_features)
        random_samples = np.arange(n_samples)
        random_state.shuffle(random_samples)
        if self.m != "all":
            random_samples = random_samples[:self.m]
        else:
            self.m = n_samples
        feature_scale = X.max(axis=0) - X.min(axis=0)
        values, counts = np.unique(y, return_counts=True)
        probabilities = counts / n_samples
        label_probs = dict(zip(values, probabilities))
        X = (X - np.mean(X, axis=0)) / feature_scale

        for i in random_samples:
            X_sample = X[i].reshape(1, -1)
            y_sample = y[i]
            sample_class_prob = label_probs[y_sample]
            matching_indices = np.where((y == y_sample) & (np.arange(n_samples) != i))
            X_hits = X[matching_indices]
            self._update_weights(X_hits, X_sample, feature_scale)
            for other_label, other_class_prob in label_probs.items():
                if other_label == y_sample:
                    continue

                class_indices = np.where(y == other_label)
                X_class = X[class_indices]
                prior_ratio = (other_class_prob / (1 - sample_class_prob))
                self._update_weights(X_class, X_sample, feature_scale, prior_ratio=prior_ratio, hits=False)

    def _update_weights(self, X_class, X_sample, feature_scale, prior_ratio=1, hits=True):
        neighbor_indices = self._get_neighbors(X_class, X_sample)
        num_neighbors = len(neighbor_indices)
        distances = self._get_distances(num_neighbors)
        for i, neighbor_idx in enumerate(neighbor_indices):
            X_neighbor = X_class[neighbor_idx]
            diff = self._get_differences(X_sample.reshape(-1), X_neighbor, feature_scale)
            update_values = prior_ratio * (diff / self.m) * distances[i]
            if hits:
                self.weights_ -= update_values
            else:
                self.weights_ += update_values

    def _get_neighbors(self, X_class, X_sample):
        nbrs = NearestNeighbors(
            n_neighbors=self.k,
            algorithm=self.neighbors_algorithm,
            metric=self.neighbors_metric,
            metric_params=self.metric_params,
        )
        nbrs.fit(X_class)
        neighbor_indices = nbrs.kneighbors(X_sample, n_neighbors=self.k, return_distance=False)
        return neighbor_indices.reshape(-1)

    def _get_differences(self, X_sample, X_neighbor, feature_scale):
        return np.abs(X_sample - X_neighbor)

    def _get_distances(self, neighbor_count):
        dist = np.exp(-((np.arange(1, neighbor_count + 1) / self.sigma) ** 2))
        return dist / dist.sum()

    def _get_support_mask(self):
        pass
