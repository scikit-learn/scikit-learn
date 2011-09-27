# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Vlad Niculae <vlad@vene.ro>

import numpy as np

from ..decomposition import PCA
from ..decomposition.dict_learning import BaseDictionaryLearning
from ..cluster import KMeans
from ..metrics.pairwise import euclidean_distances


class KMeansCoder(BaseDictionaryLearning):
    """K-means based dictionary learning

    The fit method receives an array of signals, whitens them using
    a PCA transform and run a KMeans algorithm to extract the patch
    cluster centers.

    The input is then correlated with each individual "patch-center"
    treated as a convolution kernel. The transformation can be done
    using various sparse coding methods, tresholding or the triangle
    k-means non-linearity.

    This estimator only implements the unsupervised feature extraction
    part of the referenced paper. Image classification can then be
    performed by training a linear SVM model on the output of this
    estimator.

    Parameters
    ----------
    n_atoms: int,
        number of centers extracted by the kmeans algorithm

    whiten: boolean, optional: default True
        perform a whitening PCA on the data at feature extraction time

    n_components: int, optional: default None
        number of components to keep after whitening individual samples

    max_iter: int, default 100
        maximum number of iterations to run the k-means algorithm

    n_init, int, default 1
        number of times to initialize the k-means algorithm in order to
        avoid convergence to local optima

    n_prefit: int, default 5
        dimension of reduced curriculum space in which to prefit the k-means
        algorithm for increased performance.
        This is used only when `whiten=True`.

    tol: float, default 1e-4
        tolerance for numerical errors

    local_contrast: boolean, optional: default True
        perform local contrast normalization on the extracted patch

    verbose: bool, default False
        whether to display verbose output

    transform_algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        Algorithm used to transform the data.
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection X.T * Y

    transform_n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    split_sign: bool, default False
        whether to split the transformed feature vectors into positive and
        negative components, such that the downstream classification algorithms
        can assign different weights depending on the sign

    n_jobs: int,
        number of parallel jobs to run

    Attributes
    ----------
    components_: array of shape n_atoms, n_features
        centers extracted by k-means from the patch space

    Reference
    ---------
    An Analysis of Single-Layer Networks in Unsupervised Feature Learning
    Adam Coates, Honglak Lee and Andrew Ng. In NIPS*2010 Workshop on
    Deep Learning and Unsupervised Feature Learning.
    http://robotics.stanford.edu/~ang/papers/nipsdlufl10-AnalysisSingleLayerUnsupervisedFeatureLearning.pdf

    """
    def __init__(self, n_atoms, whiten=True, n_components=None,
                 max_iter=100, n_init=1, n_prefit=5, tol=1e-4,
                 local_contrast=True, n_drop_components=0, verbose=False,
                 transform_algorithm='omp', transform_n_nonzero_coefs=None,
                 transform_alpha=None, split_sign=False, n_jobs=1):
        self.n_atoms = n_atoms
        self.whiten = whiten
        self.max_iter = max_iter
        self.n_init = n_init
        self.n_components = n_components
        self.local_contrast = local_contrast
        self.n_prefit = n_prefit
        self.verbose = verbose
        self.tol = tol
        self.transform_algorithm = transform_algorithm
        self.transform_n_nonzero_coefs = transform_n_nonzero_coefs
        self.transform_alpha = transform_alpha
        self.split_sign = split_sign
        self.n_jobs = n_jobs

    def local_contrast_normalization(self, X):
        """Normalize the patch-wise variance of the signal

        Parameters
        ----------
        X: array-like, shape n_samples, n_features
            Data to be normalized

        Returns
        -------
        X:
            Data after individual normalization of the samples
        """
        # XXX: this should probably be extracted somewhere more general
        # center all colour channels together
        X = X.reshape((X.shape[0], -1))
        X -= X.mean(axis=1)[:, None]

        X_std = X.std(axis=1)
        # Cap the divisor to avoid amplifying samples that are essentially
        # a flat surface into full-contrast salt-and-pepper garbage.
        # the actual value is a wild guess
        # This trick is credited to N. Pinto
        min_divisor = (2 * X_std.min() + X_std.mean()) / 3
        X /= np.maximum(min_divisor, X_std).reshape(
            (X.shape[0], 1))
        return X

    def fit(self, X, y=None, **kwargs):
        """Fit the encoder on a collection of data, e.g. image patches.

        Parameters
        ----------
        X: array-like, shape: n_samples, n_features
            the patch data to be fitted

        Returns
        -------
        self: object
            Returns the object itself
        """
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape
        # normalize each patch individually
        if self.local_contrast:
            if self.verbose:
                print "Local contrast normalization of the data"
            X = self.local_contrast_normalization(X)

        # kmeans model to find the filters
        if self.verbose:
            print "About to extract atoms from %d samples" % n_samples
        kmeans = KMeans(k=self.n_atoms, init='k-means++',
                        max_iter=self.max_iter, n_init=self.n_init,
                        tol=self.tol, verbose=self.verbose)

        if self.whiten:
            if self.verbose:
                print "Whitening PCA of the samples"
            self.pca = pca = PCA(whiten=True, n_components=self.n_components)
            pca.fit(X)
            X = pca.transform(X)

            # compute the KMeans centers
            if 0 < self.n_prefit < n_features:
                if self.verbose:
                    print "First KMeans in simplified curriculum space"
                # starting the kmeans on a the projection to the first singular
                # components: curriculum learning trick by Andrej Karpathy
                kmeans.fit(X[:, :self.n_prefit])

                # warm restart by padding previous centroids with zeros
                # with full dimensionality this time
                kmeans.init = np.zeros((self.n_atoms, n_features),
                                       dtype=kmeans.cluster_centers_.dtype)
                kmeans.init[:, :self.n_prefit] = kmeans.cluster_centers_
                if self.verbose:
                    print "Second KMeans in full whitened sample space"
                kmeans.set_params(n_init=1).fit(X)
            else:
                if self.verbose:
                    print "KMeans in full original sample space"
                # regular kmeans fit (without the curriculum trick)
                kmeans.fit(X)

            # project back the centers in original, non-whitened space (useful
            # for qualitative inspection of the filters)
            self.components_ = self.pca.inverse_transform(
                                                       kmeans.cluster_centers_)
        else:
            # find the kernel in the raw original dimensional space
            # TODO: experiment with component wise scaling too
            self.pca = None
            kmeans.fit(X)
            self.components_ = kmeans.cluster_centers_

        self.kmeans = kmeans
        self.inertia_ = kmeans.inertia_
        return self

    def transform(self, X, y=None):
        """Map a collection of samples into the feature space

        """ + BaseDictionaryLearning.transform.__doc__
        if self.local_contrast:
            # TODO: make it inplace by default explictly
            X = self.local_contrast_normalization(X)

        return BaseDictionaryLearning.transform(self, X, y)
