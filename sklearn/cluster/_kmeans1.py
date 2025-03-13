import numpy as np
from sklearn.cluster import KMeans as KMeansBase
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted

class KMeans(KMeansBase):
    """
    KMeans clustering with optional automatic cluster selection using the elbow method.

    Parameters
    ----------
    n_clusters : int or str, default=8
        The number of clusters to form as well as the number of centroids to generate.
        If "auto", the number of clusters is determined using the elbow method.
    max_auto_clusters : int, default=10
        Maximum number of clusters to test when n_clusters="auto".
        Ignored if n_clusters is an integer.
    init : {'k-means++', 'random'}, callable or array-like of shape (n_clusters, n_features), default='k-means++'
        Method for initialization.
    n_init : int, default=10
        Number of times the k-means algorithm will be run with different centroid seeds.
    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm for a single run.
    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference in the cluster centers.
    verbose : int, default=0
        Verbosity mode.
    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization.
    copy_x : bool, default=True
        When pre-computing distances, copy the data to avoid modifying it.
    algorithm : {'lloyd', 'elkan'}, default='lloyd'
        K-means algorithm to use. 'lloyd' is the standard algorithm, 'elkan' is more efficient for dense data.

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers.
    labels_ : ndarray of shape (n_samples,)
        Labels of each point.
    inertia_ : float
        Sum of squared distances of samples to their closest cluster center.
    n_iter_ : int
        Number of iterations run.
    """

    def __init__(self, n_clusters=8, *, max_auto_clusters=10, init='k-means++', 
                 n_init=10, max_iter=300, tol=1e-4, verbose=0, random_state=None, 
                 copy_x=True, algorithm='lloyd'):
        super().__init__(
            n_clusters=n_clusters, init=init, n_init=n_init, max_iter=max_iter, 
            tol=tol, verbose=verbose, random_state=random_state, copy_x=copy_x, 
            algorithm=algorithm
        )
        self.max_auto_clusters = max_auto_clusters

    def _auto_cluster_selection(self, X):
        """
        Determine the optimal number of clusters using the elbow method.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training data.

        Returns
        -------
        optimal_k : int
            The optimal number of clusters based on the elbow method.
        """
        inertias = []
        cluster_range = range(1, self.max_auto_clusters + 1)
        
        # Fit KMeans for each k and collect inertia
        for k in cluster_range:
            model = KMeansBase(
                n_clusters=k, init=self.init, n_init=self.n_init, 
                max_iter=self.max_iter, tol=self.tol, verbose=self.verbose, 
                random_state=self.random_state, copy_x=self.copy_x, 
                algorithm=self.algorithm  # Use 'lloyd' or user-specified algorithm
            )
            model.fit(X)
            inertias.append(model.inertia_)

        # Calculate the elbow point using second derivative
        diffs = np.diff(inertias)  # First derivative
        diffs2 = np.diff(diffs)    # Second derivative
        if len(diffs2) == 0:  # Edge case: too few points to compute curvature
            return 1
        optimal_k = cluster_range[np.argmax(diffs2) + 1]  # Max curvature
        return optimal_k

    def fit(self, X, y=None, sample_weight=None):
        """
        Compute k-means clustering, with automatic cluster selection if n_clusters="auto".

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster.
        y : Ignored
            Not used, present for API consistency.
        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        # Validate input
        X = check_array(X, accept_sparse='csr', dtype=[np.float64, np.float32])

        # Handle n_clusters="auto"
        if isinstance(self.n_clusters, str) and self.n_clusters == "auto":
            if self.max_auto_clusters < 1:
                raise ValueError("max_auto_clusters must be positive")
            self.n_clusters = self._auto_cluster_selection(X)
            if self.verbose:
                print(f"Auto-selected {self.n_clusters} clusters using elbow method")

        # Proceed with standard KMeans fitting
        return super().fit(X, y, sample_weight)

# Example usage (for testing locally):
if __name__ == "__main__":
    from sklearn.datasets import make_blobs
    X, _ = make_blobs(n_samples=300, centers=4, cluster_std=1, random_state=42)
    model = KMeans(n_clusters="auto", max_auto_clusters=10, random_state=42)
    model.fit(X)
    print(f"Optimal clusters: {model.n_clusters}, Inertia: {model.inertia_}")