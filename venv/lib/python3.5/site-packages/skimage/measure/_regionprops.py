# coding: utf-8
from __future__ import division
from math import sqrt, atan2, pi as PI
import itertools
from warnings import warn
import numpy as np
from scipy import ndimage as ndi

from ._label import label
from . import _moments


from functools import wraps

__all__ = ['regionprops', 'perimeter']


XY_TO_RC_DEPRECATION_MESSAGE = (
    'regionprops and image moments (including moments, normalized moments, '
    'central moments, and inertia tensor) of 2D images will change from xy '
    'coordinates to rc coordinates in version 0.16.\nSee '
    'http://scikit-image.org/docs/0.14.x/release_notes_and_installation.html#deprecations '
    'for details on how to avoid this message.'
)
STREL_4 = np.array([[0, 1, 0],
                    [1, 1, 1],
                    [0, 1, 0]], dtype=np.uint8)
STREL_8 = np.ones((3, 3), dtype=np.uint8)
STREL_26_3D = np.ones((3, 3, 3), dtype=np.uint8)
PROPS = {
    'Area': 'area',
    'BoundingBox': 'bbox',
    'BoundingBoxArea': 'bbox_area',
    'CentralMoments': 'moments_central',
    'Centroid': 'centroid',
    'ConvexArea': 'convex_area',
    # 'ConvexHull',
    'ConvexImage': 'convex_image',
    'Coordinates': 'coords',
    'Eccentricity': 'eccentricity',
    'EquivDiameter': 'equivalent_diameter',
    'EulerNumber': 'euler_number',
    'Extent': 'extent',
    # 'Extrema',
    'FilledArea': 'filled_area',
    'FilledImage': 'filled_image',
    'HuMoments': 'moments_hu',
    'Image': 'image',
    'Label': 'label',
    'MajorAxisLength': 'major_axis_length',
    'MaxIntensity': 'max_intensity',
    'MeanIntensity': 'mean_intensity',
    'MinIntensity': 'min_intensity',
    'MinorAxisLength': 'minor_axis_length',
    'Moments': 'moments',
    'NormalizedMoments': 'moments_normalized',
    'Orientation': 'orientation',
    'Perimeter': 'perimeter',
    # 'PixelIdxList',
    # 'PixelList',
    'Solidity': 'solidity',
    # 'SubarrayIdx'
    'WeightedCentralMoments': 'weighted_moments_central',
    'WeightedCentroid': 'weighted_centroid',
    'WeightedHuMoments': 'weighted_moments_hu',
    'WeightedMoments': 'weighted_moments',
    'WeightedNormalizedMoments': 'weighted_moments_normalized'
}

PROP_VALS = set(PROPS.values())


def _cached(f):
    @wraps(f)
    def wrapper(obj):
        cache = obj._cache
        prop = f.__name__

        if not ((prop in cache) and obj._cache_active):
            cache[prop] = f(obj)

        return cache[prop]

    return wrapper


def only2d(method):
    @wraps(method)
    def func2d(self, *args, **kwargs):
        if self._ndim > 2:
            raise NotImplementedError('Property %s is not implemented for '
                                      '3D images' % method.__name__)
        return method(self, *args, **kwargs)
    return func2d


class _RegionProperties(object):
    """Please refer to `skimage.measure.regionprops` for more information
    on the available region properties.
    """

    def __init__(self, slice, label, label_image, intensity_image,
                 cache_active, coordinates):

        if intensity_image is not None:
            if not intensity_image.shape == label_image.shape:
                raise ValueError('Label and intensity image must have the'
                                 'same shape.')

        self.label = label

        self._slice = slice
        self._label_image = label_image
        self._intensity_image = intensity_image

        self._cache_active = cache_active
        self._cache = {}
        self._ndim = label_image.ndim
        # Note: in PR 2603, we added support for nD moments in regionprops.
        # Many properties used xy coordinates, instead of rc. This attribute
        # helps with the deprecation process and should be removed in 0.16.
        if label_image.ndim > 2 or coordinates == 'rc':
            self._use_xy_warning = False
            self._transpose_moments = False
        elif coordinates == 'xy':
            self._use_xy_warning = False  # don't warn if 'xy' given explicitly
            self._transpose_moments = True
        elif coordinates is None:
            self._use_xy_warning = True
            self._transpose_moments = True
        else:
            raise ValueError('Incorrect value for regionprops coordinates: %s.'
                             ' Possible values are: "rc", "xy", or None')

    @_cached
    def area(self):
        return np.sum(self.image)

    def bbox(self):
        """
        Returns
        -------
        A tuple of the bounding box's start coordinates for each dimension,
        followed by the end coordinates for each dimension
        """
        return tuple([self._slice[i].start for i in range(self._ndim)] +
                     [self._slice[i].stop for i in range(self._ndim)])

    def bbox_area(self):
        return self.image.size

    def centroid(self):
        return tuple(self.coords.mean(axis=0))

    @_cached
    def convex_area(self):
        return np.sum(self.convex_image)

    @_cached
    def convex_image(self):
        from ..morphology.convex_hull import convex_hull_image
        return convex_hull_image(self.image)

    def coords(self):
        indices = np.nonzero(self.image)
        return np.vstack([indices[i] + self._slice[i].start
                          for i in range(self._ndim)]).T

    @only2d
    def eccentricity(self):
        l1, l2 = self.inertia_tensor_eigvals
        if l1 == 0:
            return 0
        return sqrt(1 - l2 / l1)

    def equivalent_diameter(self):
        if self._ndim == 2:
            return sqrt(4 * self.area / PI)
        elif self._ndim == 3:
            return (6 * self.area / PI) ** (1. / 3)

    def euler_number(self):
        euler_array = self.filled_image != self.image
        _, num = label(euler_array, connectivity=self._ndim, return_num=True,
                       background=0)
        return -num + 1

    def extent(self):
        return self.area / self.image.size

    def filled_area(self):
        return np.sum(self.filled_image)

    @_cached
    def filled_image(self):
        structure = np.ones((3,) * self._ndim)
        return ndi.binary_fill_holes(self.image, structure)

    @_cached
    def image(self):
        return self._label_image[self._slice] == self.label

    @_cached
    def inertia_tensor(self):
        mu = self.moments_central
        return _moments.inertia_tensor(self.image, mu)

    @_cached
    def inertia_tensor_eigvals(self):
        return _moments.inertia_tensor_eigvals(self.image,
                                               T=self.inertia_tensor)

    @_cached
    def intensity_image(self):
        if self._intensity_image is None:
            raise AttributeError('No intensity image specified.')
        return self._intensity_image[self._slice] * self.image

    def _intensity_image_double(self):
        return self.intensity_image.astype(np.double)

    def local_centroid(self):
        M = self.moments
        if self._transpose_moments:
            M = M.T
        return tuple(M[tuple(np.eye(self._ndim, dtype=int))] /
                     M[(0,) * self._ndim])

    def max_intensity(self):
        return np.max(self.intensity_image[self.image])

    def mean_intensity(self):
        return np.mean(self.intensity_image[self.image])

    def min_intensity(self):
        return np.min(self.intensity_image[self.image])

    def major_axis_length(self):
        l1 = self.inertia_tensor_eigvals[0]
        return 4 * sqrt(l1)

    def minor_axis_length(self):
        l2 = self.inertia_tensor_eigvals[-1]
        return 4 * sqrt(l2)

    @_cached
    def moments(self):
        M = _moments.moments(self.image.astype(np.uint8), 3)
        if self._use_xy_warning:
            warn(XY_TO_RC_DEPRECATION_MESSAGE)
        if self._transpose_moments:
            M = M.T
        return M

    @_cached
    def moments_central(self):
        mu = _moments.moments_central(self.image.astype(np.uint8),
                                      self.local_centroid, order=3)
        if self._use_xy_warning:
            warn(XY_TO_RC_DEPRECATION_MESSAGE)
        if self._transpose_moments:
            mu = mu.T
        return mu

    @only2d
    def moments_hu(self):
        return _moments.moments_hu(self.moments_normalized)

    @_cached
    def moments_normalized(self):
        return _moments.moments_normalized(self.moments_central, 3)

    @only2d
    def orientation(self):
        a, b, b, c = self.inertia_tensor.flat
        if a - c == 0:
            if b < 0:
                return -PI / 4.
            else:
                return PI / 4.
        else:
            return -0.5 * atan2(-2 * b, (a - c))

    @only2d
    def perimeter(self):
        return perimeter(self.image, 4)

    def solidity(self):
        return self.area / self.convex_area

    def weighted_centroid(self):
        ctr = self.weighted_local_centroid
        return tuple(idx + slc.start
                     for idx, slc in zip(ctr, self._slice))

    def weighted_local_centroid(self):
        M = self.weighted_moments
        return (M[tuple(np.eye(self._ndim, dtype=int))] /
                M[(0,) * self._ndim])

    @_cached
    def weighted_moments(self):
        return _moments.moments(self._intensity_image_double(), 3)

    @_cached
    def weighted_moments_central(self):
        ctr = self.weighted_local_centroid
        return _moments.moments_central(self._intensity_image_double(),
                                        center=ctr, order=3)

    @only2d
    def weighted_moments_hu(self):
        return _moments.moments_hu(self.weighted_moments_normalized)

    @_cached
    def weighted_moments_normalized(self):
        return _moments.moments_normalized(self.weighted_moments_central, 3)

    def __iter__(self):
        props = PROP_VALS

        if self._intensity_image is None:
            unavailable_props = ('intensity_image',
                                 'max_intensity',
                                 'mean_intensity',
                                 'min_intensity',
                                 'weighted_moments',
                                 'weighted_moments_central',
                                 'weighted_centroid',
                                 'weighted_local_centroid',
                                 'weighted_moments_hu',
                                 'weighted_moments_normalized')

            props = props.difference(unavailable_props)

        return iter(sorted(props))

    def __getitem__(self, key):
        value = getattr(self, key, None)
        if value is not None:
            return value
        else:  # backwards compatability
            return getattr(self, PROPS[key])

    def __eq__(self, other):
        if not isinstance(other, _RegionProperties):
            return False

        for key in PROP_VALS:
            try:
                # so that NaNs are equal
                np.testing.assert_equal(getattr(self, key, None),
                                        getattr(other, key, None))
            except AssertionError:
                return False

        return True


def regionprops(label_image, intensity_image=None, cache=True,
                coordinates=None):
    """Measure properties of labeled image regions.

    Parameters
    ----------
    label_image : (N, M) ndarray
        Labeled input image. Labels with value 0 are ignored.
    intensity_image : (N, M) ndarray, optional
        Intensity (i.e., input) image with same size as labeled image.
        Default is None.
    cache : bool, optional
        Determine whether to cache calculated properties. The computation is
        much faster for cached properties, whereas the memory consumption
        increases.
    coordinates : 'rc' or 'xy', optional
        Coordinate conventions for 2D images. (Only 'rc' coordinates are
        supported for 3D images.)

    Returns
    -------
    properties : list of RegionProperties
        Each item describes one labeled region, and can be accessed using the
        attributes listed below.

    Notes
    -----
    The following properties can be accessed as attributes or keys:

    **area** : int
        Number of pixels of region.
    **bbox** : tuple
        Bounding box ``(min_row, min_col, max_row, max_col)``.
        Pixels belonging to the bounding box are in the half-open interval
        ``[min_row; max_row)`` and ``[min_col; max_col)``.
    **bbox_area** : int
        Number of pixels of bounding box.
    **centroid** : array
        Centroid coordinate tuple ``(row, col)``.
    **convex_area** : int
        Number of pixels of convex hull image.
    **convex_image** : (H, J) ndarray
        Binary convex hull image which has the same size as bounding box.
    **coords** : (N, 2) ndarray
        Coordinate list ``(row, col)`` of the region.
    **eccentricity** : float
        Eccentricity of the ellipse that has the same second-moments as the
        region. The eccentricity is the ratio of the focal distance
        (distance between focal points) over the major axis length.
        The value is in the interval [0, 1).
        When it is 0, the ellipse becomes a circle.
    **equivalent_diameter** : float
        The diameter of a circle with the same area as the region.
    **euler_number** : int
        Euler characteristic of region. Computed as number of objects (= 1)
        subtracted by number of holes (8-connectivity).
    **extent** : float
        Ratio of pixels in the region to pixels in the total bounding box.
        Computed as ``area / (rows * cols)``
    **filled_area** : int
        Number of pixels of filled region.
    **filled_image** : (H, J) ndarray
        Binary region image with filled holes which has the same size as
        bounding box.
    **image** : (H, J) ndarray
        Sliced binary region image which has the same size as bounding box.
    **inertia_tensor** : (2, 2) ndarray
        Inertia tensor of the region for the rotation around its mass.
    **inertia_tensor_eigvals** : tuple
        The two eigen values of the inertia tensor in decreasing order.
    **intensity_image** : ndarray
        Image inside region bounding box.
    **label** : int
        The label in the labeled input image.
    **local_centroid** : array
        Centroid coordinate tuple ``(row, col)``, relative to region bounding
        box.
    **major_axis_length** : float
        The length of the major axis of the ellipse that has the same
        normalized second central moments as the region.
    **max_intensity** : float
        Value with the greatest intensity in the region.
    **mean_intensity** : float
        Value with the mean intensity in the region.
    **min_intensity** : float
        Value with the least intensity in the region.
    **minor_axis_length** : float
        The length of the minor axis of the ellipse that has the same
        normalized second central moments as the region.
    **moments** : (3, 3) ndarray
        Spatial moments up to 3rd order::

            m_ji = sum{ array(x, y) * x^j * y^i }

        where the sum is over the `x`, `y` coordinates of the region.
    **moments_central** : (3, 3) ndarray
        Central moments (translation invariant) up to 3rd order::

            mu_ji = sum{ array(x, y) * (x - x_c)^j * (y - y_c)^i }

        where the sum is over the `x`, `y` coordinates of the region,
        and `x_c` and `y_c` are the coordinates of the region's centroid.
    **moments_hu** : tuple
        Hu moments (translation, scale and rotation invariant).
    **moments_normalized** : (3, 3) ndarray
        Normalized moments (translation and scale invariant) up to 3rd order::

            nu_ji = mu_ji / m_00^[(i+j)/2 + 1]

        where `m_00` is the zeroth spatial moment.
    **orientation** : float
        Angle between the X-axis and the major axis of the ellipse that has
        the same second-moments as the region. Ranging from `-pi/2` to
        `pi/2` in counter-clockwise direction.
    **perimeter** : float
        Perimeter of object which approximates the contour as a line
        through the centers of border pixels using a 4-connectivity.
    **solidity** : float
        Ratio of pixels in the region to pixels of the convex hull image.
    **weighted_centroid** : array
        Centroid coordinate tuple ``(row, col)`` weighted with intensity
        image.
    **weighted_local_centroid** : array
        Centroid coordinate tuple ``(row, col)``, relative to region bounding
        box, weighted with intensity image.
    **weighted_moments** : (3, 3) ndarray
        Spatial moments of intensity image up to 3rd order::

            wm_ji = sum{ array(x, y) * x^j * y^i }

        where the sum is over the `x`, `y` coordinates of the region.
    **weighted_moments_central** : (3, 3) ndarray
        Central moments (translation invariant) of intensity image up to
        3rd order::

            wmu_ji = sum{ array(x, y) * (x - x_c)^j * (y - y_c)^i }

        where the sum is over the `x`, `y` coordinates of the region,
        and `x_c` and `y_c` are the coordinates of the region's weighted
        centroid.
    **weighted_moments_hu** : tuple
        Hu moments (translation, scale and rotation invariant) of intensity
        image.
    **weighted_moments_normalized** : (3, 3) ndarray
        Normalized moments (translation and scale invariant) of intensity
        image up to 3rd order::

            wnu_ji = wmu_ji / wm_00^[(i+j)/2 + 1]

        where ``wm_00`` is the zeroth spatial moment (intensity-weighted area).

    Each region also supports iteration, so that you can do::

      for prop in region:
          print(prop, region[prop])

    See Also
    --------
    label

    References
    ----------
    .. [1] Wilhelm Burger, Mark Burge. Principles of Digital Image Processing:
           Core Algorithms. Springer-Verlag, London, 2009.
    .. [2] B. Jähne. Digital Image Processing. Springer-Verlag,
           Berlin-Heidelberg, 6. edition, 2005.
    .. [3] T. H. Reiss. Recognizing Planar Objects Using Invariant Image
           Features, from Lecture notes in computer science, p. 676. Springer,
           Berlin, 1993.
    .. [4] http://en.wikipedia.org/wiki/Image_moment

    Examples
    --------
    >>> from skimage import data, util
    >>> from skimage.measure import label
    >>> img = util.img_as_ubyte(data.coins()) > 110
    >>> label_img = label(img, connectivity=img.ndim)
    >>> props = regionprops(label_img)
    >>> # centroid of first labeled object
    >>> props[0].centroid
    (22.729879860483141, 81.912285234465827)
    >>> # centroid of first labeled object
    >>> props[0]['centroid']
    (22.729879860483141, 81.912285234465827)

    """

    label_image = np.squeeze(label_image)

    if label_image.ndim not in (2, 3):
        raise TypeError('Only 2-D and 3-D images supported.')

    if not np.issubdtype(label_image.dtype, np.integer):
        raise TypeError('Label image must be of integer type.')

    regions = []

    objects = ndi.find_objects(label_image)
    for i, sl in enumerate(objects):
        if sl is None:
            continue

        label = i + 1

        props = _RegionProperties(sl, label, label_image, intensity_image,
                                  cache, coordinates=coordinates)
        regions.append(props)

    return regions


def perimeter(image, neighbourhood=4):
    """Calculate total perimeter of all objects in binary image.

    Parameters
    ----------
    image : array
        Binary image.
    neighbourhood : 4 or 8, optional
        Neighborhood connectivity for border pixel determination.

    Returns
    -------
    perimeter : float
        Total perimeter of all objects in binary image.

    References
    ----------
    .. [1] K. Benkrid, D. Crookes. Design and FPGA Implementation of
           a Perimeter Estimator. The Queen's University of Belfast.
           http://www.cs.qub.ac.uk/~d.crookes/webpubs/papers/perimeter.doc
    """
    if neighbourhood == 4:
        strel = STREL_4
    else:
        strel = STREL_8
    image = image.astype(np.uint8)
    eroded_image = ndi.binary_erosion(image, strel, border_value=0)
    border_image = image - eroded_image

    perimeter_weights = np.zeros(50, dtype=np.double)
    perimeter_weights[[5, 7, 15, 17, 25, 27]] = 1
    perimeter_weights[[21, 33]] = sqrt(2)
    perimeter_weights[[13, 23]] = (1 + sqrt(2)) / 2

    perimeter_image = ndi.convolve(border_image, np.array([[10, 2, 10],
                                                           [ 2, 1,  2],
                                                           [10, 2, 10]]),
                                   mode='constant', cval=0)

    # You can also write
    # return perimeter_weights[perimeter_image].sum()
    # but that was measured as taking much longer than bincount + np.dot (5x
    # as much time)
    perimeter_histogram = np.bincount(perimeter_image.ravel(), minlength=50)
    total_perimeter = np.dot(perimeter_histogram, perimeter_weights)
    return total_perimeter


def _parse_docs():
    import re
    import textwrap

    doc = regionprops.__doc__
    matches = re.finditer('\*\*(\w+)\*\* \:.*?\n(.*?)(?=\n    [\*\S]+)',
                          doc, flags=re.DOTALL)
    prop_doc = dict((m.group(1), textwrap.dedent(m.group(2))) for m in matches)

    return prop_doc


def _install_properties_docs():
    prop_doc = _parse_docs()

    for p in [member for member in dir(_RegionProperties)
              if not member.startswith('_')]:
        try:
            getattr(_RegionProperties, p).__doc__ = prop_doc[p]
        except AttributeError:
            # For Python 2.x
            getattr(_RegionProperties, p).im_func.__doc__ = prop_doc[p]

        setattr(_RegionProperties, p, property(getattr(_RegionProperties, p)))


_install_properties_docs()
