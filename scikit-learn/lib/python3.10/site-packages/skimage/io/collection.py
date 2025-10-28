"""Data structures to hold collections of images, with optional caching."""

import os
from glob import glob
import re
from collections.abc import Sequence
from copy import copy

import numpy as np
from PIL import Image

from tifffile import TiffFile


__all__ = [
    'MultiImage',
    'ImageCollection',
    'concatenate_images',
    'imread_collection_wrapper',
]


def concatenate_images(ic):
    """Concatenate all images in the image collection into an array.

    Parameters
    ----------
    ic : an iterable of images
        The images to be concatenated.

    Returns
    -------
    array_cat : ndarray
        An array having one more dimension than the images in `ic`.

    See Also
    --------
    ImageCollection.concatenate
    MultiImage.concatenate

    Raises
    ------
    ValueError
        If images in `ic` don't have identical shapes.

    Notes
    -----
    ``concatenate_images`` receives any iterable object containing images,
    including ImageCollection and MultiImage, and returns a NumPy array.
    """
    all_images = [image[np.newaxis, ...] for image in ic]
    try:
        array_cat = np.concatenate(all_images)
    except ValueError:
        raise ValueError('Image dimensions must agree.')
    return array_cat


def alphanumeric_key(s):
    """Convert string to list of strings and ints that gives intuitive sorting.

    Parameters
    ----------
    s : string

    Returns
    -------
    k : a list of strings and ints

    Examples
    --------
    >>> alphanumeric_key('z23a')
    ['z', 23, 'a']
    >>> filenames = ['f9.10.png', 'e10.png', 'f9.9.png', 'f10.10.png',
    ...              'f10.9.png']
    >>> sorted(filenames)
    ['e10.png', 'f10.10.png', 'f10.9.png', 'f9.10.png', 'f9.9.png']
    >>> sorted(filenames, key=alphanumeric_key)
    ['e10.png', 'f9.9.png', 'f9.10.png', 'f10.9.png', 'f10.10.png']
    """
    k = [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', s)]
    return k


def _is_multipattern(input_pattern):
    """Helping function. Returns True if pattern contains a tuple, list, or a
    string separated with os.pathsep."""
    # Conditions to be accepted by ImageCollection:
    has_str_ospathsep = isinstance(input_pattern, str) and os.pathsep in input_pattern
    not_a_string = not isinstance(input_pattern, str)
    has_iterable = isinstance(input_pattern, Sequence)
    has_strings = all(isinstance(pat, str) for pat in input_pattern)

    is_multipattern = has_str_ospathsep or (
        not_a_string and has_iterable and has_strings
    )
    return is_multipattern


class ImageCollection:
    """Load and manage a collection of image files.

    Parameters
    ----------
    load_pattern : str or list of str
        Pattern string or list of strings to load. The filename path can be
        absolute or relative.
    conserve_memory : bool, optional
        If True, :class:`skimage.io.ImageCollection` does not keep more than one in
        memory at a specific time. Otherwise, images will be cached once they are loaded.

    Other parameters
    ----------------
    load_func : callable
        ``imread`` by default. See Notes below.
    **load_func_kwargs : dict
        Any other keyword arguments are passed to `load_func`.

    Attributes
    ----------
    files : list of str
        If a pattern string is given for `load_pattern`, this attribute
        stores the expanded file list. Otherwise, this is equal to
        `load_pattern`.

    Notes
    -----
    Note that files are always returned in alphanumerical order. Also note that slicing
    returns a new :class:`skimage.io.ImageCollection`, *not* a view into the data.

    ImageCollection image loading can be customized through
    `load_func`. For an ImageCollection ``ic``, ``ic[5]`` calls
    ``load_func(load_pattern[5])`` to load that image.

    For example, here is an ImageCollection that, for each video provided,
    loads every second frame::

      import imageio.v3 as iio3
      import itertools

      def vidread_step(f, step):
          vid = iio3.imiter(f)
          return list(itertools.islice(vid, None, None, step)

      video_file = 'no_time_for_that_tiny.gif'
      ic = ImageCollection(video_file, load_func=vidread_step, step=2)

      ic  # is an ImageCollection object of length 1 because 1 video is provided

      x = ic[0]
      x[5]  # the 10th frame of the first video

    Alternatively, if `load_func` is provided and `load_pattern` is a
    sequence, an :class:`skimage.io.ImageCollection` of corresponding length will
    be created, and the individual images will be loaded by calling `load_func` with the
    matching element of the `load_pattern` as its first argument. In this
    case, the elements of the sequence do not need to be names of existing
    files (or strings at all). For example, to create an :class:`skimage.io.ImageCollection`
    containing 500 images from a video::

      class FrameReader:
          def __init__ (self, f):
              self.f = f
          def __call__ (self, index):
              return iio3.imread(self.f, index=index)

      ic = ImageCollection(range(500), load_func=FrameReader('movie.mp4'))

      ic  # is an ImageCollection object of length 500

    Another use of `load_func` would be to convert all images to ``uint8``::

      def imread_convert(f):
          return imread(f).astype(np.uint8)

      ic = ImageCollection('/tmp/*.png', load_func=imread_convert)

    Examples
    --------
    >>> import imageio.v3 as iio3
    >>> import skimage.io as io

    # Where your images are located
    >>> data_dir = os.path.join(os.path.dirname(__file__), '../data')

    >>> coll = io.ImageCollection(data_dir + '/chess*.png')
    >>> len(coll)
    2
    >>> coll[0].shape
    (200, 200)

    >>> image_col = io.ImageCollection([f'{data_dir}/*.png', '{data_dir}/*.jpg'])

    >>> class MultiReader:
    ...     def __init__ (self, f):
    ...         self.f = f
    ...     def __call__ (self, index):
    ...         return iio3.imread(self.f, index=index)
    ...
    >>> filename = data_dir + '/no_time_for_that_tiny.gif'
    >>> ic = io.ImageCollection(range(24), load_func=MultiReader(filename))
    >>> len(image_col)
    23
    >>> isinstance(ic[0], np.ndarray)
    True
    """

    def __init__(
        self, load_pattern, conserve_memory=True, load_func=None, **load_func_kwargs
    ):
        """Load and manage a collection of images."""
        self._files = []
        if _is_multipattern(load_pattern):
            if isinstance(load_pattern, str):
                load_pattern = load_pattern.split(os.pathsep)
            for pattern in load_pattern:
                self._files.extend(glob(pattern))
            self._files = sorted(self._files, key=alphanumeric_key)
        elif isinstance(load_pattern, str):
            self._files.extend(glob(load_pattern))
            self._files = sorted(self._files, key=alphanumeric_key)
        elif isinstance(load_pattern, Sequence) and load_func is not None:
            self._files = list(load_pattern)
        else:
            raise TypeError('Invalid pattern as input.')

        if load_func is None:
            from ._io import imread

            self.load_func = imread
            self._numframes = self._find_images()
        else:
            self.load_func = load_func
            self._numframes = len(self._files)
            self._frame_index = None

        if conserve_memory:
            memory_slots = 1
        else:
            memory_slots = self._numframes

        self._conserve_memory = conserve_memory
        self._cached = None

        self.load_func_kwargs = load_func_kwargs
        self.data = np.empty(memory_slots, dtype=object)

    @property
    def files(self):
        return self._files

    @property
    def conserve_memory(self):
        return self._conserve_memory

    def _find_images(self):
        index = []
        for fname in self._files:
            if fname.lower().endswith(('.tiff', '.tif')):
                with open(fname, 'rb') as f:
                    img = TiffFile(f)
                    index += [(fname, i) for i in range(len(img.pages))]
            else:
                try:
                    im = Image.open(fname)
                    im.seek(0)
                except OSError:
                    continue
                i = 0
                while True:
                    try:
                        im.seek(i)
                    except EOFError:
                        break
                    index.append((fname, i))
                    i += 1
                if hasattr(im, 'fp') and im.fp:
                    im.fp.close()
        self._frame_index = index
        return len(index)

    def __getitem__(self, n):
        """Return selected image(s) in the collection.

        Loading is done on demand.

        Parameters
        ----------
        n : int or slice
            The image number to be returned, or a slice selecting the images
            and ordering to be returned in a new ImageCollection.

        Returns
        -------
        img : ndarray or :class:`skimage.io.ImageCollection`
            The `n`-th image in the collection, or a new ImageCollection with
            the selected images.
        """
        if hasattr(n, '__index__'):
            n = n.__index__()

        if not isinstance(n, (int, slice)):
            raise TypeError('slicing must be with an int or slice object')

        if isinstance(n, int):
            n = self._check_imgnum(n)
            idx = n % len(self.data)

            if (self.conserve_memory and n != self._cached) or (self.data[idx] is None):
                kwargs = self.load_func_kwargs
                if self._frame_index:
                    fname, img_num = self._frame_index[n]
                    if img_num is not None:
                        kwargs['img_num'] = img_num
                    try:
                        self.data[idx] = self.load_func(fname, **kwargs)
                    # Account for functions that do not accept an img_num kwarg
                    except TypeError as e:
                        if "unexpected keyword argument 'img_num'" in str(e):
                            del kwargs['img_num']
                            self.data[idx] = self.load_func(fname, **kwargs)
                        else:
                            raise
                else:
                    self.data[idx] = self.load_func(self.files[n], **kwargs)
                self._cached = n

            return self.data[idx]
        else:
            # A slice object was provided, so create a new ImageCollection
            # object. Any loaded image data in the original ImageCollection
            # will be copied by reference to the new object.  Image data
            # loaded after this creation is not linked.
            fidx = range(self._numframes)[n]
            new_ic = copy(self)

            if self._frame_index:
                new_ic._files = [self._frame_index[i][0] for i in fidx]
                new_ic._frame_index = [self._frame_index[i] for i in fidx]
            else:
                new_ic._files = [self._files[i] for i in fidx]

            new_ic._numframes = len(fidx)

            if self.conserve_memory:
                if self._cached in fidx:
                    new_ic._cached = fidx.index(self._cached)
                    new_ic.data = np.copy(self.data)
                else:
                    new_ic.data = np.empty(1, dtype=object)
            else:
                new_ic.data = self.data[fidx]
            return new_ic

    def _check_imgnum(self, n):
        """Check that the given image number is valid."""
        num = self._numframes
        if -num <= n < num:
            n = n % num
        else:
            raise IndexError(f"There are only {num} images in the collection")
        return n

    def __iter__(self):
        """Iterate over the images."""
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        """Number of images in collection."""
        return self._numframes

    def __str__(self):
        return str(self.files)

    def reload(self, n=None):
        """Clear the image cache.

        Parameters
        ----------
        n : None or int
            Clear the cache for this image only. By default, the
            entire cache is erased.

        """
        self.data = np.empty_like(self.data)

    def concatenate(self):
        """Concatenate all images in the collection into an array.

        Returns
        -------
        ar : np.ndarray
            An array having one more dimension than the images in `self`.

        See Also
        --------
        skimage.io.concatenate_images

        Raises
        ------
        ValueError
            If images in the :class:`skimage.io.ImageCollection` do not have identical
            shapes.
        """
        return concatenate_images(self)


def imread_collection_wrapper(imread):
    def imread_collection(load_pattern, conserve_memory=True):
        """Return an `ImageCollection` from files matching the given pattern.

        Note that files are always stored in alphabetical order. Also note that
        slicing returns a new ImageCollection, *not* a view into the data.

        See `skimage.io.ImageCollection` for details.

        Parameters
        ----------
        load_pattern : str or list
            Pattern glob or filenames to load. The path can be absolute or
            relative.  Multiple patterns should be separated by a colon,
            e.g. ``/tmp/work/*.png:/tmp/other/*.jpg``.  Also see
            implementation notes below.
        conserve_memory : bool, optional
            If True, never keep more than one in memory at a specific
            time.  Otherwise, images will be cached once they are loaded.

        """
        return ImageCollection(
            load_pattern, conserve_memory=conserve_memory, load_func=imread
        )

    return imread_collection


class MultiImage(ImageCollection):
    """A class containing all frames from multi-frame TIFF images.

    Parameters
    ----------
    load_pattern : str or list of str
        Pattern glob or filenames to load. The path can be absolute or
        relative.
    conserve_memory : bool, optional
        Whether to conserve memory by only caching the frames of a single
        image. Default is True.

    Notes
    -----
    `MultiImage` returns a list of image-data arrays. In this
    regard, it is very similar to `ImageCollection`, but the two differ in
    their treatment of multi-frame images.

    For a TIFF image containing N frames of size WxH, `MultiImage` stores
    all frames of that image as a single element of shape `(N, W, H)` in the
    list. `ImageCollection` instead creates N elements of shape `(W, H)`.

    For an animated GIF image, `MultiImage` reads only the first frame, while
    `ImageCollection` reads all frames by default.

    Examples
    --------
    # Where your images are located
    >>> data_dir = os.path.join(os.path.dirname(__file__), '../data')

    >>> multipage_tiff = data_dir + '/multipage.tif'
    >>> multi_img = MultiImage(multipage_tiff)
    >>> len(multi_img)  # multi_img contains one element
    1
    >>> multi_img[0].shape  # this element is a two-frame image of shape:
    (2, 15, 10)

    >>> image_col = ImageCollection(multipage_tiff)
    >>> len(image_col)  # image_col contains two elements
    2
    >>> for frame in image_col:
    ...     print(frame.shape)  # each element is a frame of shape (15, 10)
    ...
    (15, 10)
    (15, 10)
    """

    def __init__(self, filename, conserve_memory=True, dtype=None, **imread_kwargs):
        """Load a multi-img."""
        from ._io import imread

        self._filename = filename
        super().__init__(filename, conserve_memory, load_func=imread, **imread_kwargs)

    @property
    def filename(self):
        return self._filename
