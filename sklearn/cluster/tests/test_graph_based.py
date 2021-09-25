"""
Tests for Graph-Based clustering algorithms
"""

import unittest
from typing import Dict, Tuple

import numpy as np
from sklearn import datasets
from parameterized import parameterized
from sklearn.metrics import rand_score
from sklearn.preprocessing import StandardScaler

from sklearn.cluster import (
    ConnectedComponentsClustering,
    SpanTreeConnectedComponentsClustering,
)


np.random.seed(0)


def prepare_sklearn_clustering_datasets() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:

    # ============
    # https://scikit-learn.org/stable/auto_examples/cluster/plot_cluster_comparison.html
    # Generate datasets. We choose the size big enough to see the scalability
    # of the algorithms, but not too big to avoid too long running times
    # ============

    np.random.seed(0)

    n_samples = 1500
    noisy_circles = datasets.make_circles(
        n_samples=n_samples,
        factor=0.5,
        noise=0.05,
    )
    noisy_moons = datasets.make_moons(
        n_samples=n_samples,
        noise=0.05,
    )
    blobs = datasets.make_blobs(
        n_samples=n_samples,
        random_state=8,
    )
    no_structure = (
        np.random.rand(n_samples, 2),
        np.zeros(n_samples, dtype=int),
    )

    # Anisotropicly distributed data
    random_state = 170
    X, y = datasets.make_blobs(
        n_samples=n_samples,
        random_state=random_state,
    )
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    aniso = (X_aniso, y)

    # blobs with varied variances
    varied = datasets.make_blobs(
        n_samples=n_samples,
        cluster_std=[1.0, 2.5, 0.5],
        random_state=random_state,
    )

    sklearn_clustering_datasets = {}

    sklearn_clustering_datasets["noisy_circles"] = noisy_circles
    sklearn_clustering_datasets["noisy_moons"] = noisy_moons
    sklearn_clustering_datasets["blobs"] = blobs
    sklearn_clustering_datasets["no_structure"] = no_structure
    sklearn_clustering_datasets["aniso"] = aniso
    sklearn_clustering_datasets["varied"] = varied

    return sklearn_clustering_datasets


class TestConnectedComponentsClustering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sklearn_clustering_datasets = prepare_sklearn_clustering_datasets()
        cls.sklearn_clustering_datasets = sklearn_clustering_datasets

    @parameterized.expand(
        [
            (1.5, 1, [0, 0, 0]),
            (1.25, 1, [0, 0, 0]),
            (1.0, 3, [0, 1, 2]),
            (0.75, 3, [0, 1, 2]),
        ]
    )
    def test_mini_dataset(self, threshold, n_clusters, labels):

        X = np.array([[0, 1], [1, 0], [1, 1]])

        clustering = ConnectedComponentsClustering(
            threshold=threshold,
            metric="euclidean",
            n_jobs=-1,
        )

        clustering.fit(X)

        labels_pred = clustering.labels_
        n_clusters_pred = len(np.unique(labels_pred))

        self.assertEqual(n_clusters_pred, n_clusters)
        self.assertTrue(np.allclose(labels_pred, labels))

        labels_pred_2 = clustering.fit_predict(X)

        self.assertTrue(np.allclose(labels_pred_2, labels))

    def test_sklearn_clustering_datasets(self):

        for dataset_name, dataset in self.sklearn_clustering_datasets.items():

            X, y = dataset

            # normalize dataset for easier parameter selection
            X = StandardScaler().fit_transform(X)

            clustering = ConnectedComponentsClustering(
                threshold=0.275,
                metric="euclidean",
                n_jobs=-1,
            )

            labels_pred = clustering.fit_predict(X)

            score = rand_score(labels_true=y, labels_pred=labels_pred)

            if dataset_name in ["aniso", "varied"]:
                self.assertNotEqual(score, 1.0)
            else:
                self.assertEqual(score, 1.0)


class TestSpanTreeConnectedComponentsClustering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        # sklearn_clustering_datasets
        sklearn_clustering_datasets = prepare_sklearn_clustering_datasets()
        cls.sklearn_clustering_datasets = sklearn_clustering_datasets

        # n_cluster for each dataset
        datasets_n_cluster = {}

        datasets_n_cluster["noisy_circles"] = 2
        datasets_n_cluster["noisy_moons"] = 2
        datasets_n_cluster["blobs"] = 3
        datasets_n_cluster["no_structure"] = 1
        datasets_n_cluster["aniso"] = 3
        datasets_n_cluster["varied"] = 3

        cls.datasets_n_cluster = datasets_n_cluster

    @parameterized.expand([(1,), (2,), (3,)])
    def test_mini_dataset(self, n_clusters):

        X = np.array([[0, 1], [1, 0], [1, 1]])

        clustering = SpanTreeConnectedComponentsClustering(
            n_clusters=n_clusters,
            metric="euclidean",
            n_jobs=-1,
        )

        clustering.fit(X)

        labels_pred = clustering.labels_
        n_clusters_pred = len(np.unique(labels_pred))

        self.assertEqual(n_clusters_pred, n_clusters)

        labels_pred_2 = clustering.fit_predict(X)
        n_clusters_pred = len(np.unique(labels_pred_2))

        self.assertEqual(n_clusters_pred, n_clusters)

    def test_sklearn_clustering_datasets(self):

        for dataset_name, dataset in self.sklearn_clustering_datasets.items():

            X, y = dataset

            # normalize dataset for easier parameter selection
            X = StandardScaler().fit_transform(X)

            clustering = SpanTreeConnectedComponentsClustering(
                n_clusters=self.datasets_n_cluster[dataset_name],
                metric="euclidean",
                n_jobs=-1,
            )

            labels_pred = clustering.fit_predict(X)

            score = rand_score(labels_true=y, labels_pred=labels_pred)

            if dataset_name in ["aniso", "varied"]:
                self.assertNotEqual(score, 1.0)
            else:
                self.assertEqual(score, 1.0)
