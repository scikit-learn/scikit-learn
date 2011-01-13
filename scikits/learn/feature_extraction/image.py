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
from ..cluster.k_means_ import all_pairs_l2_distance_squared

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

    Examples
    --------

    >>> image = np.arange(16).reshape((1, 4, 4))
    >>> image
    array([[[ 0,  1,  2,  3],
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11],
            [12, 13, 14, 15]]])

    >>> patches = extract_patches2d(image, (4, 4), (2, 2))
    >>> patches.shape
    (4, 2, 2)

    >>> patches[0]
    array([[0, 1],
           [4, 5]])

    >>> patches[1]
    array([[2, 3],
           [6, 7]])

    >>> patches[2]
    array([[ 8,  9],
           [12, 13]])

    >>> patches[3]
    array([[10, 11],
           [14, 15]])

    >>> patches = extract_patches2d(image, (4, 4), (2, 2), (0, 1))
    >>> patches.shape
    (2, 2, 2)

    >>> patches[0]
    array([[1, 2],
           [5, 6]])

    >>> patches[1]
    array([[ 9, 10],
           [13, 14]])
    """
    i_h, i_w = image_size[:2]
    p_h, p_w = patch_size

    images = np.atleast_2d(images)
    n_images = images.shape[0]
    images = images.reshape((n_images, i_h, i_w, -1))
    n_colors = images.shape[-1]

    # handle offsets and compute remainder to find total number of patches
    o_h, o_w = offsets
    n_h, r_h = divmod(i_h - o_h,  p_h)
    n_w, r_w = divmod(i_w - o_w,  p_w)
    n_patches = n_images * n_h * n_w

    # extract the image areas that can be sliced into whole patches
    max_h = -r_h or None
    max_w = -r_w or None
    images = images[:, o_h:max_h, o_w:max_w, :]

    # put the color dim before the sliceable dims
    images = images.transpose((0, 3, 1, 2))

    # slice the last two dims of the images into patches
    patches = images.reshape((n_images, n_colors, n_h, p_h, n_w, p_w))

    # reorganize the dims to put n_image, n_h, and n_w at the end so that
    # reshape will combine them all in n_patches
    patches = patches.transpose((3, 5, 1, 0, 2, 4))
    patches = patches.reshape((p_h, p_w, n_colors, n_patches))

    # one more transpose to put the n_patches as the first dom
    patches = patches.transpose((3, 0, 1, 2))

    # remove the color dimension if useless
    if patches.shape[-1] == 1:
        return patches.reshape((n_patches, p_h, p_w))
    else:
        return patches

def most_square_shape(N):
    """return integers (height, width) with area N that is closest to square, 
    given that width >= height.

    This returns nice tiling dimensions for plotting filters.
    """
    for i in xrange(int(np.sqrt(N)),0, -1):
        if 0 == N % i:
            return (i, N/i)
def tile_images(imgs, flip=False, return_img=True, scale_each=False, eps=1e-8):
    """Return an image or RGB ndarray in which the `imgs` are arranged (tiled) for viewing.
    """
    # returns an ndarray that is 
    w = imgs.transpose(0,3,1,2) # old code compatibility...
    if w.shape[1] != 3:
        raise NotImplementedError('not rgb', w.shape)
    if w.shape[2] != w.shape[3]:
        # i haven't debugged this well enough to know that i've got it working for nonsquare
        raise NotImplementedError('not square', w.shape)
    def move_to_01(X):
        return (X - X.min()) / (np.max(X.max() - X.min(), eps))

    if not scale_each:
        w = np.asarray(255 * move_to_01(w), dtype='uint8')
    trows, tcols= most_square_shape(w.shape[0])
    outrows = trows * w.shape[2] + trows-1
    outcols = tcols * w.shape[3] + tcols-1
    out = np.zeros((outrows, outcols,3), dtype='uint8')

    tr_stride= 1+w.shape[1]
    for tr in range(trows):
        for tc in range(tcols):
            # this is supposed to flip the filters back into the image coordinates
            # as well as put the channels in the right place, but I don't know if it
            # really does that
            w_rc = w[tr*tcols+tc].transpose(1,2,0)[::-1 if flip else 1,::-1 if flip else 1]
            if scale_each:
                w_rc = np.asarray(255 * move_to_01(w_rc), dtype='uint8')
            out[tr*(1+w.shape[2]):tr*(1+w.shape[2])+w.shape[2],
                    tc*(1+w.shape[3]):tc*(1+w.shape[3])+w.shape[3]] = w_rc
    if return_img:
        from PIL import Image
        return Image.fromarray(out, 'RGB')
    else:
        return out

class ConvolutionalKMeansEncoder(BaseEstimator):
    """Unsupervised sparse feature extractor for 2D images

    The fit method extracts patches from the images, whiten them using
    a PCA transform and run a KMeans algorithm to extract the patch
    cluster centers.

    The input is then correlated with each individual "patch-center"
    treated as a convolution kernel. The activations are then
    sparse-encoded using the "triangle" variant of a KMeans such that
    for each center c[k] we define a transformation function f_k over
    the input patches:

        f_k(x) = max(0, np.mean(z) - z[k])

    where z[k] = linalg.norm(x - c[k]) and x is an input patch.

    Activations are then sum-pooled over the 4 quadrants of the original
    image space.

    The transform operation is performed by applying the triangle kmeans
    sparse-encoding and sum-pooling.

    This estimator only implements the unsupervised feature extraction
    part of the referenced paper. Image classification can then be
    performed by training a linear SVM model on the output of this
    estimator.

    Parameters
    ----------
    n_centers: int, optional: default 400
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

    def __init__(self, n_centers=400, image_size=None, patch_size=6,
                 step_size=1, whiten=True, n_components=None,
                 pools=2, max_iter=1, n_init=1, kmeans_init_algo='k-means++',
                 verbose=1, n_drop_components=0,
                 center_mode='all',
                 lcn_cap_divisor=True,
                 ):
        self.n_centers = n_centers
        self.patch_size = patch_size
        self.step_size = step_size
        self.whiten = whiten
        self.pools = pools
        self.image_size = image_size
        self.max_iter = max_iter
        self.n_init = n_init
        self.n_components = n_components
        self.n_drop_components = n_drop_components
        self.kmeans_init_algo=kmeans_init_algo
        self.verbose = verbose
        self.center_mode = center_mode
        self.lcn_cap_divisor=lcn_cap_divisor
        if center_mode not in ('all', 'none', 'channel'):
            raise ValueError('invalid center mode', center_mode)

    def _check_images(self, X):
        """Check that X can seen as a consistent collection of images"""
        X = np.atleast_2d(X)
        n_samples = X.shape[0]

        if self.image_size is None:
            if len(X.shape) > 3:
                self.image_size = X.shape[1:3]
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

    def local_contrast_normalization(self, patches):
        if self.center_mode == 'all': # center all colour channels together
            patches = patches.reshape((patches.shape[0], -1))
            patches -= patches.mean(axis=1)[:,None]
        elif self.center_mode == 'channel': #center each colour channel individually
            patches -= patches.mean(axis=2).mean(axis=1).reshape((n_patches, 1, 1, 3))
            patches = patches.reshape((n_patches, -1))
        elif self.center_mode == 'none':
            patches = patches.reshape((n_patches, -1))
        else:
            assert False

        patches_var = np.sqrt((patches**2).mean(axis=1))
        if self.lcn_cap_divisor:
            min_divisor = (patches_var.min() + patches_var.mean()) / 2
        else:
            min_divisor = 0
        patches /= np.maximum(min_divisor, patches_var).reshape((patches.shape[0],1))
        return patches

    def fit(self, X):
        """Fit the feature extractor on a collection of 2D images"""
        X = self._check_images(X)

        # step 1: extract the patches
        offsets = [(o, o)
                   for o in range(0, self.patch_size - 1, self.step_size)]
        patch_size = (self.patch_size, self.patch_size)

        # this list of patches does not copy the memory allocated for raw
        # image data
        patches_by_offset = [extract_patches2d(
            X, self.image_size, patch_size, offsets=o) for o in offsets]

        # TODO: compute pca and kmeans taking other offsets into account to
        # step 2: whiten the patch space
        patches = patches_by_offset[0]
        n_patches = patches.shape[0]
        patches = self.local_contrast_normalization(patches)

        self.patches_= patches

        # step 3: compute the KMeans centers
        kmeans = KMeans(k=self.n_centers, init=self.kmeans_init_algo,
                        max_iter=self.max_iter, n_init=self.n_init, verbose=self.verbose)
        self.kmeans_ = kmeans
        # TODO: when whitening is enabled, implement curriculum learnin by
        # starting the kmeans on a the projection to the first singular
        # components and increase the number component with warm restarts by
        # padding previous centroids with zeros to keep up with the increasing
        # dim

        if self.whiten:
            self.pca = PCA(whiten=True, n_components=self.n_components)
            self.pca.fit(patches)
            patches_pca = self.pca.transform(patches)
            if self.n_drop_components:
                assert patches.shape[1] == self.pca.components_.shape[0]
                self.pca.components_[:,:self.n_drop_components] = 0
                patches_pca[:,:self.n_drop_components] = 0
            kmeans.fit(patches_pca)
            self.inertia_ = kmeans.inertia_

            # step 4: project back the centers in original, non-whitened space
            def unpca(X):
                return (np.dot(X* self.pca.components_coefs_**2,
                self.pca.components_.T) + self.pca.mean_)
            self.kernels_ = unpca(kmeans.cluster_centers_)

            self.patches_unpca_= unpca(patches_pca)
        else:
            self.pca = None
            kmeans.fit(patches)
            self.inertia_ = kmeans.inertia_
            self.kernels_ = kmeans.cluster_centers_
            self.patches_unpca_= patches

        return self

    def transform(self, X, rtype='float32'):
        """Map a collection of 2D images into the feature space"""
        #X = self._check_images(X)
        nX, nXrows, nXcols, nXchannels = X.shape
        if nXchannels != 3:
            raise NotImplementedError('should be checking channels throughout')
        out_features = np.zeros((X.shape[0], 2, 2, self.kernels_.shape[0]), dtype=rtype)
        for r in xrange(nXrows-self.patch_size+1):
            if self.verbose:
                print '%s::transform() row %i/%i'%(
                        self.__class__.__name__, r+1, nXrows-self.patch_size+1)
            for c in xrange(nXcols-self.patch_size+1):
                patches = X[:,r:r+self.patch_size,c:c+self.patch_size,:]
                n_patches = patches.shape[0]
                patches = self.local_contrast_normalization(patches)

                if self.whiten:
                    patches = self.pca.transform(patches)

                distances = np.sqrt(all_pairs_l2_distance_squared(
                        patches.reshape(patches.shape[0],-1),
                        self.kmeans_.cluster_centers_))

                if 1: #triangle features
                    features = np.maximum(0, distances.mean(axis=1)[:,None] - distances)
                elif 0: # hard assignment features
                    raise NotImplementedError()
                else: #posterior features
                    #TODO: include the cluster size as a component prior
                    features = np.exp(-0.5*distances)

                out_features[:,r//16,c//16,:] += features
        return out_features

    def tile_kernels(self, scale_each=False):
        shp = self.kernels_.shape[0], self.patch_size, self.patch_size, 3
        return tile_images(self.kernels_.reshape(shp),scale_each=scale_each)
    def tile_patches(self, N=256, scale_each=False):
        shp = self.patches_[:N].shape[0], self.patch_size, self.patch_size, 3
        return tile_images(self.patches_[:N].reshape(shp),scale_each=scale_each)
    def tile_patches_unpca(self, N=256, scale_each=False):
        shp = self.patches_unpca_[:N].shape[0], self.patch_size, self.patch_size, 3
        return tile_images(self.patches_unpca_[:N].reshape(shp),scale_each=scale_each)
    def tile_pca_components(self, scale_each=False):
        patches = self.pca.components_.T
        shp = patches.shape[0], self.patch_size, self.patch_size, 3
        return tile_images(patches.reshape(shp),scale_each=scale_each)
