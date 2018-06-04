from __future__ import absolute_import, division, print_function

from glob import glob
import os

try:
    from skimage.io import imread as sk_imread
except ImportError:
    pass

from .core import Array
from ..base import tokenize


def add_leading_dimension(x):
    return x[None, ...]


def imread(filename, imread=None, preprocess=None):
    """ Read a stack of images into a dask array

    Parameters
    ----------

    filename: string
        A globstring like 'myfile.*.png'
    imread: function (optional)
        Optionally provide custom imread function.
        Function should expect a filename and produce a numpy array.
        Defaults to ``skimage.io.imread``.
    preprocess: function (optional)
        Optionally provide custom function to preprocess the image.
        Function should expect a numpy array for a single image.

    Examples
    --------

    >>> from dask.array.image import imread
    >>> im = imread('2015-*-*.png')  # doctest: +SKIP
    >>> im.shape  # doctest: +SKIP
    (365, 1000, 1000, 3)

    Returns
    -------

    Dask array of all images stacked along the first dimension.  All images
    will be treated as individual chunks
    """
    imread = imread or sk_imread
    filenames = sorted(glob(filename))
    if not filenames:
        raise ValueError("No files found under name %s" % filename)

    name = 'imread-%s' % tokenize(filenames, map(os.path.getmtime, filenames))

    sample = imread(filenames[0])
    if preprocess:
        sample = preprocess(sample)

    keys = [(name, i) + (0,) * len(sample.shape) for i in range(len(filenames))]
    if preprocess:
        values = [(add_leading_dimension, (preprocess, (imread, fn)))
                  for fn in filenames]
    else:
        values = [(add_leading_dimension, (imread, fn))
                  for fn in filenames]
    dsk = dict(zip(keys, values))

    chunks = ((1, ) * len(filenames), ) + tuple((d, ) for d in sample.shape)

    return Array(dsk, name, chunks, sample.dtype)
