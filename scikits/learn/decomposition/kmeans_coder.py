# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Vlad Niculae <vlad@vene.ro>

import numpy as np

from scikits.learn.decomposition import PCA
from scikits.learn.cluster import KMeans

class KMeansCoder():
    def __init__(self, n_centers=400, whiten=True, n_components=None,
                 n_pools=2, max_iter=100, n_init=1, n_prefit=5, tol=1e-4,
                 local_contrast=True, n_drop_components=0, verbose=False):
        self.n_centers = n_centers
        self.whiten = whiten
        self.n_pools = n_pools
        self.max_iter = max_iter
        self.n_init = n_init
        self.n_components = n_components
        self.local_contrast = local_contrast
        self.n_prefit = n_prefit
        self.verbose = verbose
        self.tol = tol
        self.n_drop_components = n_drop_components

    def local_contrast_normalization(self, patches):
        """Normalize the patch-wise variance of the signal"""
        # center all colour channels together
        patches = patches.reshape((patches.shape[0], -1))
        patches -= patches.mean(axis=1)[:, None]

        patches_std = patches.std(axis=1)
        # Cap the divisor to avoid amplifying patches that are essentially
        # a flat surface into full-contrast salt-and-pepper garbage.
        # the actual value is a wild guess
        # This trick is credited to N. Pinto
        min_divisor = (2 * patches_std.min() + patches_std.mean()) / 3
        patches /= np.maximum(min_divisor, patches_std).reshape(
            (patches.shape[0], 1))
        return patches

    def fit(self, X, y=None, **kwargs):
        """Fit the encoder on a collection of data, e.g. image patches.

        Parameters
        ----------
        X: array-like, shape: n_samples, *patch_size
            the patch data to be fitted
        """
        patches = np.atleast_2d(X)
        n_patches = len(patches)
        patches = patches.reshape((n_patches, -1))
        # normalize each patch individually
        if self.local_contrast:
            if self.verbose:
                print "Local contrast normalization of the patches"
            patches = self.local_contrast_normalization(patches)

        # kmeans model to find the filters
        if self.verbose:
            print "About to extract filters from %d patches" % n_patches
        kmeans = KMeans(k=self.n_centers, init='k-means++',
                        max_iter=self.max_iter, n_init=self.n_init,
                        tol=self.tol, verbose=self.verbose)

        if self.whiten:
            # whiten the patch space
            if self.verbose:
                print "Whitening PCA of the patches"
            self.pca = pca = PCA(whiten=True, n_components=self.n_components)
            pca.fit(patches)

            # implement a band-pass filter by dropping the first eigen
            # values which are generally low frequency components
            if self.n_drop_components:
                pca.components_ = pca.components_[self.n_drop_components:, :]
            patches = pca.transform(patches)

            # compute the KMeans centers
            if 0 < self.n_prefit < patches.shape[1]:
                if self.verbose:
                    print "First KMeans in simplified curriculum space"
                # starting the kmeans on a the projection to the first singular
                # components: curriculum learning trick by Andrej Karpathy
                kmeans.fit(patches[:, :self.n_prefit])

                # warm restart by padding previous centroids with zeros
                # with full dimensionality this time
                kmeans.init = np.zeros((self.n_centers, patches.shape[1]),
                                       dtype=kmeans.cluster_centers_.dtype)
                kmeans.init[:, :self.n_prefit] = kmeans.cluster_centers_
                if self.verbose:
                    print "Second KMeans in full whitened patch space"
                kmeans.fit(patches, n_init=1)
            else:
                if self.verbose:
                    print "KMeans in full original patch space"
                # regular kmeans fit (without the curriculum trick)
                kmeans.fit(patches)

            # project back the centers in original, non-whitened space (useful
            # for qualitative inspection of the filters)
            self.filters_ = self.pca.inverse_transform(kmeans.cluster_centers_)
        else:
            # find the kernel in the raw original dimensional space
            # TODO: experiment with component wise scaling too
            self.pca = None
            kmeans.fit(patches)
            self.filters_ = kmeans.cluster_centers_

        self.kmeans = kmeans
        self.inertia_ = kmeans.inertia_
        return self

    def transform(self, X):
        """Map a collection of 2D images into the feature space"""
        raise NotImplemented("Ignore the stuff below for now")

        X = self._check_images(X)
        n_samples, n_rows, n_cols, n_channels = X.shape
        n_filters = self.filters_.shape[0]
        ps = (8, 8) #self.patch_size

        pooled_features = np.zeros((X.shape[0], self.n_pools, self.n_pools,
                                    n_filters), dtype=X.dtype)

        n_rows_adjusted = n_rows - self.patch_size + 1
        n_cols_adjusted = n_cols - self.patch_size + 1

        # TODO: the following is really slow: profile me, optimize me and then
        # parallelize me by using joblib.Parallel.map on the toplevel loop:

        for r in xrange(n_rows_adjusted):
            if self.verbose:
                print "Extracting features for row #%d/%d" % (
                    r + 1, n_rows_adjusted)

            for c in xrange(n_cols_adjusted):

                patches = X[:, r:r + ps, c:c + ps, :].reshape((n_samples, -1))

                if self.local_contrast:
                    # TODO: make it inplace by default explictly
                    patches = self.local_contrast_normalization(patches)

                if self.whiten:
                    # TODO: make it possible to pass pre-allocated array
                    patches = self.pca.transform(patches)

                # extract distance from each patch to each cluster center
                # TODO: make it possible to reuse pre-allocated distance array
                filters = self.kmeans.cluster_centers_
                distances = euclidean_distances(patches, filters)

                # triangle features
                distance_means = distances.mean(axis=1)[:, None]
                features = np.maximum(0, distance_means - distances)

                # features are pooled over image regions
                out_r = r * self.n_pools / n_rows_adjusted
                out_c = c * self.n_pools / n_cols_adjusted
                pooled_features[:, out_r, out_c, :] += features

        # downstream classifiers expect a 2 dim shape
        return pooled_features.reshape(pooled_features.shape[0], -1)
