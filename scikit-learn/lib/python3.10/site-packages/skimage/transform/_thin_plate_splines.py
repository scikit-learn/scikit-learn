import numpy as np
from scipy.spatial import distance_matrix

from .._shared.utils import check_nD


class ThinPlateSplineTransform:
    """Thin-plate spline transformation.

    Given two matching sets of points, source and destination, this class
    estimates the thin-plate spline (TPS) transformation which transforms
    each point in source into its destination counterpart.

    Attributes
    ----------
    src : (N, 2) array_like
        Coordinates of control points in source image.

    References
    ----------
    .. [1] Bookstein, Fred L. "Principal warps: Thin-plate splines and the
           decomposition of deformations," IEEE Transactions on pattern analysis
           and machine intelligence 11.6 (1989): 567–585.
           DOI:`10.1109/34.24792`
           https://user.engineering.uiowa.edu/~aip/papers/bookstein-89.pdf

    Examples
    --------
    >>> import skimage as ski

    Define source and destination control points such that they simulate
    rotating by 90 degrees and generate a meshgrid from them:

    >>> src = np.array([[0, 0], [0, 5], [5, 5], [5, 0]])
    >>> dst = np.array([[5, 0], [0, 0], [0, 5], [5, 5]])

    Estimate the transformation:

    >>> tps = ski.transform.ThinPlateSplineTransform()
    >>> tps.estimate(src, dst)
    True

    Appyling the transformation to `src` approximates `dst`:

    >>> np.round(tps(src))
    array([[5., 0.],
           [0., 0.],
           [0., 5.],
           [5., 5.]])

    Create a meshgrid to apply the transformation to:

    >>> grid = np.meshgrid(np.arange(5), np.arange(5))
    >>> grid[1]
    array([[0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1],
           [2, 2, 2, 2, 2],
           [3, 3, 3, 3, 3],
           [4, 4, 4, 4, 4]])

    >>> coords = np.vstack([grid[0].ravel(), grid[1].ravel()]).T
    >>> transformed = tps(coords)
    >>> np.round(transformed[:, 1]).reshape(5, 5).astype(int)
    array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]])
    """

    def __init__(self):
        self._estimated = False
        self._spline_mappings = None
        self.src = None

    def __call__(self, coords):
        """Estimate the transformation from a set of corresponding points.

        Parameters
        ----------
        coords : (N, 2) array_like
            x, y coordinates to transform

        Returns
        -------
        transformed_coords: (N, D) array
            Destination coordinates
        """
        if self._spline_mappings is None:
            msg = (
                "Transformation is undefined, define it by calling `estimate` "
                "before applying it"
            )
            raise ValueError(msg)
        coords = np.array(coords)

        if coords.ndim != 2 or coords.shape[1] != 2:
            msg = "Input `coords` must have shape (N, 2)"
            raise ValueError(msg)

        radial_dist = self._radial_distance(coords)
        transformed_coords = self._spline_function(coords, radial_dist)

        return transformed_coords

    @property
    def inverse(self):
        raise NotImplementedError("Not supported")

    def estimate(self, src, dst):
        """Estimate optimal spline mappings between source and destination points.

        Parameters
        ----------
        src : (N, 2) array_like
            Control points at source coordinates.
        dst : (N, 2) array_like
            Control points at destination coordinates.

        Returns
        -------
        success: bool
            True indicates that the estimation was successful.

        Notes
        -----
        The number N of source and destination points must match.
        """
        check_nD(src, 2, arg_name="src")
        check_nD(dst, 2, arg_name="dst")

        if src.shape[0] < 3 or dst.shape[0] < 3:
            msg = "Need at least 3 points in in `src` and `dst`"
            raise ValueError(msg)
        if src.shape != dst.shape:
            msg = f"Shape of `src` and `dst` didn't match, {src.shape} != {dst.shape}"
            raise ValueError(msg)

        self.src = src
        n, d = src.shape

        dist = distance_matrix(src, src)
        K = self._radial_basis_kernel(dist)
        P = np.hstack([np.ones((n, 1)), src])
        n_plus_3 = n + 3
        L = np.zeros((n_plus_3, n_plus_3), dtype=np.float32)
        L[:n, :n] = K
        L[:n, -3:] = P
        L[-3:, :n] = P.T
        V = np.vstack([dst, np.zeros((d + 1, d))])
        try:
            self._spline_mappings = np.linalg.solve(L, V)
        except np.linalg.LinAlgError:
            return False
        return True

    def _radial_distance(self, coords):
        """Compute the radial distance between input points and source points."""
        dists = distance_matrix(coords, self.src)
        return self._radial_basis_kernel(dists)

    def _spline_function(self, coords, radial_dist):
        """Estimate the spline function in X and Y directions."""
        n = self.src.shape[0]
        w = self._spline_mappings[:n]
        a = self._spline_mappings[n:]
        transformed_coords = a[0] + np.dot(coords, a[1:]) + np.dot(radial_dist, w)
        return transformed_coords

    @staticmethod
    def _radial_basis_kernel(r):
        """Compute the radial basis function for thin-plate splines.

        Parameters
        ----------
        r : (4, N) ndarray
            Input array representing the Euclidean distance between each pair of
            two collections of control points.

        Returns
        -------
        U : (4, N) ndarray
            Calculated kernel function U.
        """
        _small = 1e-8  # Small value to avoid divide-by-zero
        r_sq = r**2
        U = np.where(r == 0.0, 0.0, r_sq * np.log(r_sq + _small))
        return U
