# -*- coding: utf-8 -*-
""" Multiview Spectral Clustering.

It provides a function which produces
a single clustering assigment, but considering all the input data from
different views.
"""

___author___ = 'Maria Araceli Burgueño Caballero'
___email___ = "mburgueno@uoc.edu"
___status___ = "Pre-Production"


import numpy as np
import warnings
from numbers import Number
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.multiview.cpcmv import MVCPC
from sklearn.base import BaseEstimator
from sklearn.utils import check_array
from sklearn.cluster import KMeans


def laplacian_ng(mat):
    """ Laplacian following Ng et al formulation.

    Parameters
    ----------

    mat : ndarray
        Input data (array or matrix).

    Returns
    -------

    output : ndarray
        Laplacian formulation (matrix form).
    """

    mat = check_array(mat)

    D = np.sum(mat, axis=1)
    sqriD = np.diag(1 / np.sqrt(D))
    return np.dot(np.dot(sqriD, mat), sqriD)


def suggested_sigma(distances, neighbour=None):
    """ Gives the suggested sigma to spectral cluster given a distance matrix.

    Parameters
    ----------

    distances: ndarray
        Input data.
    neighbour: numeric, default None
        If None, select avg distance according to Luxburg criterium (avg
        distance to the log(number of samples)-th neighbour). If it is an
        integer > 1, avg distance to the neighbour-th closest sample. If it
        is a real number 0<neighbour<=1: avg distance to the
        (neighbour-th*number of samples).

    Returns
    -------

    result : numeric
        Suggested sigma according to distance matrix.
    """

    distances = check_array(distances)

    n = None
    if (neighbour is None):
        n = np.ceil(np.log(distances.shape[0])) - 1
    else:
        n = neighbour - 1

    dist_ord = np.sort(distances, axis=0)
    # Compute the mean removing NAs and infinite values just in case
    # If it is 0 then return 1 (or we will get errors)
    result = np.mean(dist_ord[int(n), np.isfinite(dist_ord[int(n), :])])
    if result == 0.0:
        result = 1.0
    return result


def distance_gaussian_similarity(distances, sigma):
    """ Given a distance matrix, compute the gaussian similarity of its values.

    Parameters
    ----------
    distances : ndarray.
        Input data.
    sigma :     numeric.
        Sigma parameter for the gaussian function.

    Returns
    -------

    result : ndarray.
        Gaussian similarity matrix.
    """

    distances = check_array(distances)

    my_factor = -1 / (2 * sigma**2)
    result = np.exp(distances**2 * my_factor)
    # 0's are dangerous afterwards, they should be replaced by something safer
    result[result == 0] = 1e-16
    return result


def mvsc(x, is_distance, k, sigmas=None, neighbours=None, clustering=True):
    """ Multiview spectral clustering on a list of matrices or distance matrices.

    Computes the multiview spectral clustering of data on a list of matrices
    or distance matrices (or a mix of both), supposed to be different views of
    the same data.
    In the case of plain data matrices, euclidean distance will be used to
    generate distance matrices for that data view.

    Notes
    -----

    All input views must have the same number of samples (rows).

    Parameters
    ----------

    x : list
        A list of feature matrices or distance matrices (or a mix of both).
    is_distance : array-like.
        A list or array which indicates whether a matrix with the same index
        in x is a distance matrix (true value) or not (false value).
    k : int
        Number of desired clusters.
    neighbours: Either None, an integer value or a vector of int, default: None
        They correspond to the expected number of neighbours per point, used
        to estimate the sigma values of the Gaussian radial basis function.
        If it is NULL then the default sigma computation is used (average
        distance to the log(n)-th neighbour, with n = number of samples).
        If it is a single value then the same number of neighbours is used
        on all input views, else each value in the vector is applied to the
        corresponding input view. Does not have effect if sigma is different
        from None.
    clustering: boolean
        Tells mvsc if it has to perform the clustering on the projection or to
        skip the clustering step of the algorithm.

    Returns
    -------

    tuple
        A tuple with four elements:

        clustering is a vector of integers with the clustering assignment of
        each sample (not included if clustering = FALSE)

        evalues is a matrix with the eigenvalues of the common principal
        components (CPC) step

        evectors is a matrix with the eigenvectors of the CPC step

        sigmas is a vector with the sigmas used on the Gaussian radial basis
        function of each input view.

    Raises
    ------

        ValueError: Matrices are not square matrices, k value is negative or
        data samples and is_distance parameters do not have the same length.

    """

    nviews = len(x)  # Number of input matrices
    if sigmas is not None:
        if isinstance(sigmas, Number):
            sigmas = [sigmas] * nviews
        elif len(sigmas) == 1:
            sigmas = sigmas * nviews
    if neighbours is not None:
        if isinstance(neighbours, Number):
            neighbours = [neighbours] * nviews
        elif len(neighbours) == 1:
            neighbours = neighbours * nviews
    # Placeholder to store the actual sigmas used
    my_sigmas = [0] * nviews

    # Compute the joint diagonal matrix of the similarity matrices
    # First we have to create a p x p x n array with the laplacian matrices
    # p = number of samples, n = number of views
    num_points = x[0].shape[0]
    lap_matrix = np.zeros((nviews, num_points, num_points))
    for i in np.arange(nviews):  # len(x)
        # Compute distance if necessary, then the kernel
        if not is_distance[i]:
            views_dist = euclidean_distances(x[i])
        else:
            views_dist = x[i]

        if sigmas is not None:
            my_sigmas[i] = sigmas[i]
        elif neighbours is not None:
            my_sigmas[i] = suggested_sigma(views_dist, neighbours[i])
        else:
            my_sigmas[i] = suggested_sigma(views_dist)

        view_grbf = distance_gaussian_similarity(views_dist, my_sigmas[i])
        view_lsym = laplacian_ng(view_grbf)
        lap_matrix[i] = view_lsym
        # end of for loop

    # Now we have to compute CPC on the laplacians to get the eigen values
    # and vectors
    mvcpc = MVCPC(k=k)
    cpc_result_evalues, cpc_result_evectors = mvcpc.fit_transform(
        lap_matrix)
    cpc_result_evectors = (cpc_result_evectors.T /
                           np.linalg.norm(cpc_result_evectors, axis=1)).T

    # Run KMeans on the eigenvectors (only first K columns are computed) and
    # return everything
    kmeans = KMeans(n_clusters=k, max_iter=(10 + round(np.log(num_points))))
    kmeans_clust = kmeans.fit(cpc_result_evectors[:, :k]).labels_
    return (kmeans_clust, cpc_result_evalues, cpc_result_evectors, my_sigmas)


class MVSC(BaseEstimator):
    """
    Multiview spectral clustering on a list of matrices or distance matrices.

    Computes the multiview spectral clustering of data on a list of matrices
    or distance matrices (or a mix of both), supposed to be different views of
    the same data.
    In the case of plain data matrices, euclidean distance will be used to
    generate distance matrices for that data view.

    Notes
    -----

    All input views must have the same number of samples (rows).

    Parameters
    ----------

    k : int
        Number of desired clusters.
    sigmas: Either None, an integer value or a vector of int, default: None
        They correspond to the sigma parameter in the Gaussian radial basis
        function. If it is None then the default sigma computation is used
        (average distance to the log(n)-th neighbour, with n = number of
        samples), unless neighbours has a value different from None.
        If it is a single number then the same sigma is applied to all input
        views. If it is a vector each value in it is applied to the
        corresponding input view.
    neighbours: Either None, an integer value or a vector of int, default: None
        They correspond to the expected number of neighbours per point, used
        to estimate the sigma values of the Gaussian radial basis function.
        If it is NULL then the default sigma computation is used (average
        distance to the log(n)-th neighbour, with n = number of samples).
        If it is a single value then the same number of neighbours is used
        on all input views, else each value in the vector is applied to the
        corresponding input view. Does not have effect if sigma is different
        from None.
    clustering: boolean
        Tells mvsc if it has to perform the clustering on the projection or to
        skip the clustering step of the algorithm.

    Attributes
    ----------
    embedding_ : ndarray
        Clustering of the nviews input data.
    evalues_ : ndarray
        Eigenvalues computed during spectral clustering.
    evectors_ :ndarray
        Eigenvectors computed during spectral clustering.
    sigmas_ : ndarray
        Best sigmas used for calculating Gaussian similarity.

    References
    ----------

        Ng, Andrew Y, Michael I Jordan, and Yair Weiss. 2001. “On spectral
        clustering: Analysis and an algorithm.” *Nips* 14 (14). MIT Press:
        849–56. doi:10.1.1.19.8100.

        Planck, Max, and Ulrike Von Luxburg. 2006. “A Tutorial on Spectral
        Clustering A Tutorial on Spectral Clustering.” *Statistics and
        Computing* 17 (March). Springer US: 395–416.
        doi:10.1007/s11222-007-9033-z.

        Shi, Jianbo, and Jitendra Malik. 2005. “Normalized Cuts and Image
        Segmentation Normalized Cuts and Image Segmentation.” *Pattern
        Analysis and Machine Intelligence, IEEE Transactions* on 22 (March):
        888–905. doi:10.1109/CVPR.1997.609407.

        Trendafilov, Nickolay T. 2010. “Stepwise estimation of common principal
        components.” *Computational Statistics and Data Analysis* 54 (12):
        3446–57. doi:10.1016/j.csda.2010.03.010.
    """

    def __init__(self, k=2, sigmas=None, neighbours=None, clustering=True):
        self.k = k
        self.sigmas = sigmas
        self.neighbours = neighbours
        self.clustering = clustering

    def fit(self, x, is_distance):
        """
        Computes the multiview spectral clustering and return the clustering,
        eigenvalues, eienvectors and sigmas used in the computation.

        Notes
        -----

        All input views must have the same number of samples (rows).

        Parameters
        ----------

        x : list
            A list of feature matrices or distance matrices (or a mix of both).
        is_distance : array-like.
            A list or array which indicates whether a matrix with the same
            index in x is a distance matrix (true value) or not (false value).
        k : int
            Number of desired clusters.
        """
        self.fit_transform(x, is_distance)
        return self

    def fit_transform(self, x, is_distance):
        """
        Computes the multiview spectral clustering and return the clustering,
        eigenvalues, eienvectors and sigmas used in the computation.

        Notes
        -----

        All input views must have the same number of samples (rows).

        Parameters
        ----------

        x : list
            A list of feature matrices or distance matrices (or a mix of both).
        is_distance : array-like.
            A list or array which indicates whether a matrix with the same
            index in x is a distance matrix (true value) or not (false value).


        Returns
        -------

        tuple
            A tuple with four elements:

            clustering is a vector of integers with the clustering assignment
            of each sample (not included if clustering = FALSE)

            evalues is a matrix with the eigenvalues of the common principal
            components (CPC) step

            evectors is a matrix with the eigenvectors of the CPC step

            sigmas is a vector with the sigmas used on the Gaussian radial
            basis
            function of each input view.

        Raises
        ------

            ValueError: Matrices are not square matrices, k value is negative
            or data samples and is_distance parameters do not have the same
            length.

        Examples
        --------

        >>> import numpy as np
        >>> m = np.array([[1, 4, 7], [2, 5, 8], [3, 6, 9]])
        >>> q = np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])
        >>> r = np.array([[2, 1, 8], [4, 5, 6], [3, 7, 9]]).T
        >>> matrices = [m, q, r]
        >>> is_distance = [False, False, False]
        >>> mvsc = MVSC(k=3)
        >>> mvsc.fit_transform(matrices, is_distance)
            (array([1, 2, 0]), array([[ 0.99983709, 0.99983709,  0.99947615],
                                      [ 0.49076485, 0.49076485,  0.44022256],
                                      [ 0.10945481, 0.10945481,  0.15255827]]),
            array([[-0.56674541,  0.64092999, -0.51769527],
                   [-0.61728928,  0.08583247,  0.78204011],
                   [-0.54566802, -0.76278538, -0.34699406]]),
            [1.7320508075688774, 1.7320508075688774, 5.2779168675293677])
        """

        if len(x) != len(is_distance):
            raise ValueError("Data samples and is_distance lengths does not "
                             "match. Data sample length: %d, is_distance "
                             "length: %d" % (len(x), len(is_distance)))
        if self.k > x[0].shape[0]:
            self.k = x[0].shape[0]
            warnings.warn("k is greater than matrix dimension. k=%d is "
                          "computed instead." % x[0].shape[0])
        elif self.k < 0:
            raise ValueError("k value must be between 0 and number of samples"
                             " of data matrix.")

        for i in np.arange(len(x) - 1):
            for j in np.arange(i + 1, len(x)):
                if x[i].shape[0] != x[j].shape[0]:
                    raise ValueError("Input data matrices have no same number "
                                     "of samples (rows).")

        kmeans_clust, cpcresult_evalues, cpcresult_evectors, my_sigmas = mvsc(
            x, is_distance, self.k, self.sigmas, self.neighbours,
            self.clustering)
        self.embedding_ = kmeans_clust
        self.evalues_ = cpcresult_evalues
        self.evectors_ = cpcresult_evectors
        self.sigmas_ = my_sigmas
        return (kmeans_clust, cpcresult_evalues,
                cpcresult_evectors, my_sigmas)

    def get_params(self, deep=True):
        return {"sigmas": self.sigmas, "neighbours": self.neighbours,
                "clustering": self.clustering}

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            self.setattr(parameter, value)
        return self
