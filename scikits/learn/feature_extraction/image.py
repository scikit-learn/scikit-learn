"""
Utilities to extract features from images.
"""

# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np
import math
from scipy import sparse
from ..utils.fixes import in1d
from ..base import BaseEstimator
from ..pca import PCA
from ..cluster import KMeans

################################################################################
# From an image to a graph

def _make_edges_3d(n_x, n_y, n_z=1):
    """Returns a list of edges for a 3D image.

    Parameters
    ===========
    n_x: integer
        The size of the grid in the x direction.
    n_y: integer
        The size of the grid in the y direction.
    n_z: integer, optional
        The size of the grid in the z direction, defaults to 1
    """
    vertices = np.arange(n_x*n_y*n_z).reshape((n_x, n_y, n_z))
    edges_deep = np.vstack((vertices[:, :, :-1].ravel(),
                            vertices[:, :, 1:].ravel()))
    edges_right = np.vstack((vertices[:, :-1].ravel(), vertices[:, 1:].ravel()))
    edges_down = np.vstack((vertices[:-1].ravel(), vertices[1:].ravel()))
    edges = np.hstack((edges_deep, edges_right, edges_down))
    return edges


def _compute_gradient_3d(edges, img):
    n_x, n_y, n_z = img.shape
    gradient = np.abs(img[edges[0]/(n_y*n_z), \
                                (edges[0] % (n_y*n_z))/n_z, \
                                (edges[0] % (n_y*n_z))%n_z] - \
                           img[edges[1]/(n_y*n_z), \
                                (edges[1] % (n_y*n_z))/n_z, \
                                (edges[1] % (n_y*n_z)) % n_z])
    return gradient


# XXX: Why mask the image after computing the weights?

def _mask_edges_weights(mask, edges, weights):
    """Apply a mask to weighted edges"""
    inds = np.arange(mask.size)
    inds = inds[mask.ravel()]
    ind_mask = np.logical_and(in1d(edges[0], inds),
                              in1d(edges[1], inds))
    edges, weights = edges[:, ind_mask], weights[ind_mask]
    maxval = edges.max()
    order = np.searchsorted(np.unique(edges.ravel()), np.arange(maxval+1))
    edges = order[edges]
    return edges, weights


def img_to_graph(img, mask=None, return_as=sparse.coo_matrix, dtype=None):
    """Graph of the pixel-to-pixel gradient connections

    Edges are weighted with the gradient values.

    Parameters
    ==========
    img: ndarray, 2D or 3D
        2D or 3D image
    mask : ndarray of booleans, optional
        An optional mask of the image, to consider only part of the
        pixels.
    return_as: np.ndarray or a sparse matrix class, optional
        The class to use to build the returned adjacency matrix.
    dtype: None or dtype, optional
        The data of the returned sparse matrix. By default it is the
        dtype of img
    """
    img = np.atleast_3d(img)
    if dtype is None:
        dtype = img.dtype
    n_x, n_y, n_z = img.shape
    edges   = _make_edges_3d(n_x, n_y, n_z)
    weights = _compute_gradient_3d(edges, img)
    if mask is not None:
        edges, weights = _mask_edges_weights(mask, edges, weights)
        img = img.squeeze()[mask]
    else:
        img = img.ravel()
    n_voxels = img.size
    diag_idx = np.arange(n_voxels)
    i_idx = np.hstack((edges[0], edges[1]))
    j_idx = np.hstack((edges[1], edges[0]))
    graph = sparse.coo_matrix((np.hstack((weights, weights, img)),
                              (np.hstack((i_idx, diag_idx)),
                               np.hstack((j_idx, diag_idx)))),
                              (n_voxels, n_voxels),
                              dtype=dtype)
    if return_as is np.ndarray:
        return graph.todense()
    return return_as(graph)


################################################################################
# From an image to a set of small image patches

def extract_patches2d(images, image_size, patch_size, offsets=(0, 0)):
    """Reshape a collection of 2D images into a collection of patches

    The extracted patches are not overlapping to avoid having to copy any
    memory.

    TODO: right now images are graylevel only: add channels to the right
    position according to the natural memory layout from PIL or scikits.image

    Parameters
    ----------
    images: array with shape (n_images, i_h, i_w) or (n_images, i_h * i_w)
        the original image data

    image_size: tuple of ints (i_h, i_w)
        the dimensions of the images

    patch_size: tuple of ints (p_h, p_w)
        the dimensions of one patch

    offsets: tuple of ints (o_h, o_w), optional (0, 0) by default
        location of the first extracted patch

    Returns
    -------
    patches: array with shape (n_patches, *patch_size)
    """
    i_h, i_w = image_size
    p_h, p_w = patch_size

    images = np.atleast_2d(images)
    n_images = images.shape[0]
    images = images.reshape((n_images, i_h, i_w))

    # handle offsets and compute remainder to find total number of patches
    o_h, o_w = offsets
    n_h, r_h = divmod(i_h - o_h,  p_h)
    n_w, r_w = divmod(i_w - o_w,  p_w)
    n_patches = n_images * n_h * n_w

    # extract the image areas that can be sliced into whole patches
    max_h = -r_h or None
    max_w = -r_w or None
    images = images[:, o_h:max_h, o_w:max_w]

    # slice the images into patches
    patches = images.reshape((n_images, n_h, p_h, n_w, p_w))

    # reorganize the patches into the expected shape
    patches = patches.transpose((2, 4, 0, 1, 3)).reshape((p_h, p_w, n_patches))

    # one more transpose to put the n_patches as the first dom
    return patches.transpose((2, 0, 1))


class ConvolutionalKMeansEncoder(BaseEstimator):
    """Unsupervised sparse feature extractor for 2D images

    The fit method extracts patches from the images, whiten them using
    a PCA transform and run a KMeans algorithms to extract centers of
    the same (small) shape.

    The input is then correlated with each individual "patch-center"
    treated as a convolution kernel. The activations are them sparse coded
    the "triangle" variant of a KMeans such that of center for each c[k]
    we define an transformation function f_k over the input patches:

        f_k(x) = max(0, np.mean(z) - z[k])

    where z[k] = linalg.norm(x - c[k]) and x is an input patch.

    Activations are then sum-pooled over the 4 quadrants of the original
    image space.

    The transform operation is performed by applying the convolutional,
    sum-pooling and SVM prediction steps.

    This estimator only implements the unsupervised feature extraction
    part of the referenced paper. Image classification can then be
    performed by training a linear SVM model on the output of this
    estimator.

    Parameters
    ----------
    n_centers: int, optional: default 1000
        number of centers extracted by the kmeans algorithm

    patch_size: tuple of int, optional: default 6
        the size of the square patches / convolution kernels learned by
        kmeans

    step_size: int, optional: 1
        number of pixels to shift between two consecutive patches (a.k.a.
        stride)

    whiten: boolean, optional: default True
        perform a whitening PCA on the patches at feature extraction time

    pools: int, optional: default 2
        number equal size areas to perform the sum-pooling of features
        over: pools=2 means 4 quadrants, pools=3 means 6 areas and so on

    Reference
    ---------
    An Analysis of Single-Layer Networks in Unsupervised Feature Learning
    Adam Coates, Honglak Lee and Andrew Ng. In NIPS*2010 Workshop on
    Deep Learning and Unsupervised Feature Learning.
    http://robotics.stanford.edu/~ang/papers/nipsdlufl10-AnalysisSingleLayerUnsupervisedFeatureLearning.pdf
    """

    def __init__(self, n_centers=1000, image_size=None, patch_size=6,
                 step_size=1, whiten=True, pools=2, max_iter=1, n_init=1):
        self.n_centers = n_centers
        self.patch_size = patch_size
        self.step_size = step_size
        self.whiten = True
        self.pools = pools
        self.image_size = image_size
        self.max_iter = max_iter
        self.n_init = n_init

    def _check_images(self, X):
        """Check that X can seen as a consistent collection of images"""
        X = np.atleast_2d(X)
        n_samples = X.shape[0]

        if self.image_size is None:
            if len(X.shape) == 3:
                self.image_size = X.shape[1:]
            elif len(X.shape) > 3:
                raise ValueError("%r is not a valid images shape" % (X.shape,))
            else:
                # assume square images
                _, n_features = X.shape
                size = math.sqrt(n_features)
                if size ** 2 != n_features:
                    raise ValueError("images with shape %r are not squares: "
                                     "the image size must be made explicit" %
                                     (X.shape,))
                self.image_size = (size, size)

        return X.reshape((n_samples, -1))

    def fit(self, X):
        """Fit the feature extractor on a collection of 2D images"""
        fit_transform(X)
        return self

    def fit_transform(self, X):
        """Fit the model while returning the transformed input"""
        X = self._check_images(X)

        # step 1: extract the patches
        offsets = [(o, o) for o in range(0, self.patch_size - 1, self.step_size)]
        patch_size = (self.patch_size, self.patch_size)

        # this list of patches does not copy the memory allocated for raw image
        # data
        patches_by_offset = [extract_patches2d(
            X, self.image_size, patch_size, offsets=o) for o in offsets]

        # TODO: compute pca and kmeans taking other offsets into account to
        # step 2: whiten the patch space
        patches = patches_by_offset[0]
        patches = patches.reshape((patches.shape[0], -1))
        pca = PCA(whiten=self.whiten).fit(patches)
        patches_pca = pca.transform(patches)

        # step 3: compute the KMeans centers
        kmeans = KMeans(k=self.n_centers, max_iter=self.max_iter,
                        n_init=self.n_init)
        kmeans.fit(patches_pca)
        self.inertia_ = kmeans.inertia_

        # step 4: project back the centers in original, non-whitened space
        self.kernels_ = (np.dot(kmeans.cluster_centers_, pca.components_)
                         + pca.mean_)

    def transform(self, X):
        """Map a collection of 2D images into the feature space"""
        X = self._check_images(X)
        raise NotImplementedError("implement me!")

