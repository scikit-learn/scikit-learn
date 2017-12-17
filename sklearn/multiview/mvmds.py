# -*- coding: utf-8 -*-
""" Multiview Multidimensional Scaling.

Multiview MDS (multidimensional scaling) provides a representation of the
pattern of proximities among a set of matrices. Likewise, these matrices are
composed of different attributes, which have a pattern of proximities as well.
Coding of two principal functions found in mvmds file. The first function,
preprocess mvmds, preprocess data for multiview MDs algorithm The second,
md_mds, computes the MDS algorithm itself for a set of matrices."""


import numpy as np
from sklearn.multiview import MVCPC
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.base import BaseEstimator
from sklearn.utils import check_array
import warnings


def preprocess_mvmds(values):
    """ Preprocess of a distance matrix for mvMDS (square + double center).

    Parameters
    ----------

    values: ndarray
        The distance matrix to be preprocessed.

    Returns
    -------

    p_matrix: nparray
        The preprocessed matrix (same dimensions)
    """

    values = check_array(values)

    p_matrix = values ** 2
    row_means = np.mean(p_matrix, axis=1)
    col_means = np.mean(p_matrix, axis=0)
    all_means = np.mean(p_matrix)

    p_matrix += all_means
    p_matrix.T[:, ] -= row_means[:]
    p_matrix -= col_means
    return p_matrix


def mvmds(x, is_distance, k=2):
    """
    Multiview MDS on a list of matrices.

    Multiview multidimensional scaling (mvmds) receives two or more
    feature matrices or distance matrices, according to is_distance parameters
    and produces a low-dimensional representation of the samples according to
    the combined information in all the input data.
    In the case of plain data matrices, euclidean distance will be used to
    generate distance matrices for that data view.
    Mv_MDS preprocess input matrix wih square and double center.

    Notes
    -----

    All input views must have the same number of samples (rows).

    Parameters
    ----------

    x : list
        A list of data matrices. Matrices can be raw data or distancematrices.
        In the case of plain data matrices, euclidean distance will be used to
        generate a distance matrix for that data view.
    is_distance : list
        Each boolean value indicates wheter the matrix in x with the same
        index is a distance matrix or not.
    k : int, default: 2
        Number of desired dimensions of the low-dimensional projection.

    Returns
    -------
    common : ndarray
        A n x k matrix with the k-dimensional projection, where n is the
        number of samples in the dataset.

    Raises
    ------

        ValueError: Matrices are not square matrices, k value is negative or
        data samples and is_distance parameters do not have the same length.
    """

    # Number of rows/observations
    num_obs = x[0].shape[0]
    my_mat = np.zeros((len(x), num_obs, num_obs))
    for i in np.arange(len(x)):
        if not is_distance[i]:
            # It is not a distance matrix
            my_view = euclidean_distances(x[i])
        else:
            my_view = x[i]
        my_views2 = preprocess_mvmds(my_view)
        my_mat[i] = -my_views2 / 2
    cpc = MVCPC(k=k)
    common = cpc.fit_transform(my_mat)
    return common[1]


class MVMDS(BaseEstimator):
    """Multiview Multidimensional Scaling.

    Multiview multidimensional scaling (mvmds) receives two or more
    feature matrices or distance matrices, according to is_distance parameters
    and produces a low-dimensional representation of the samples according to
    the combined information in all the input data.
    In the case of plain data matrices, euclidean distance will be used to
    generate distance matrices for that data view.
    Mv_MDS preprocess input matrix wih square and double center.

    Notes
    -----

    All input views must have the same number of samples (rows).

    Parameters
    ----------

    k : int, default: 2
        Number of desired dimensions of the low-dimensional projection.

    Attributes
    ----------

    components_ : ndarray
        Principal components of the dataset input.

    References
    ----------
        Kruskal, J B. 1964. “Multidimensional scaling by optimizing goodness
        of fit to a nonmetric hypothesis.” *Psychometrika* 29 (1): 1–27.
        doi:10.1007/BF02289565.

        Trendafilov, Nickolay T. 2010. “Stepwise estimation of common principal
        components.” *Computational Statistics and Data Analysis* 54 (12):
        3446–57. doi:10.1016/j.csda.2010.03.010.
    """

    def __init__(self, k=2):
        self.k = k

    def fit(self, x, is_distance):
        """
        Computes euclidean distances of X according to ``is_distance``,
        preprocess these data and fit them.

        Parameters
        ----------

        x : list
            A list of data matrices. Matrices can be raw data or distance
            matrices. In the case of plain data matrices, euclidean distance
            will be used to generate a distance matrix for that data view.
        is_distance : list
            Each boolean value indicates wheter the matrix in x with the same
            index is a distance matrix or not.

        """

        self.fit_transform(x, is_distance)
        return self

    def fit_transform(self, x, is_distance):
        """
        Computes euclidean distances of X according to ``is_distance``,
        preprocess these data, fit and return them.

        Parameters
        ----------

        x : list
            A list of data matrices. Matrices can be raw data or distance
            matrices. In the case of plain data matrices, euclidean distance
            will be used to generate a distance matrix for that data view.
        is_distance : list
            Each boolean value indicates wheter the matrix in x with the same
            index is a distance matrix or not.

        Returns
        -------
        common : ndarray
            A n x k matrix with the k-dimensional projection, where n is the
            number of samples in the dataset.

        Raises
        ------

            ValueError: Matrices are not square matrices, k value is negative
            or data samples and is_distance parameters do not have the same
            length.

        Examples
        --------

        >>> import numpy as np
        >>> m = np.array([[1, 4], [2, 5], [3, 6]])
        >>> r = np.array([[2, 4], [1, 5], [8, 6]])
        >>> matrices = [m, r]
        >>> is_distance = [False, False]
        >>> mv_mds = MVMDS(k=2)
        >>> mv_mds.fit_transform(matrices, is_distance)
        array([[-0.55030705 -0.60318224]
               [-0.24721761  0.77817101]
               [ 0.79752467 -0.17498877]])
        >>>
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

        common = mvmds(x, is_distance, self.k)
        self.components_ = common
        return self.components_


############################################################
#                           MAIN                           #
############################################################
# # values = np.array([[1, 2, 3, 4], [5, 6, 7, 8],
# #                    [9, 10, 11, 12], [13, 14, 15, 16]]).T
# # print(values)
# # p_matrix = preprocess_mvmds(values)
# # print(p_matrix)
# m = np.array([[1, 4], [2, 5], [3, 6]])
# r = np.array([[2, 4], [1, 5], [8, 6]])
# # q = np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])
# # # r = np.array([[2, 4, 3], [1, 5, 7], [8, 6, 9]])
# matrices = [m, r]
# distance_bool = [False, False]
# print(matrices)
# mv_mds = MVMDS()
# print(mv_mds.fit(matrices, distance_bool).components_)
# x = mv_mds.fit_transform(matrices, distance_bool)
# print(x)
