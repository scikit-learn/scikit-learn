import numpy as np
import six

from ..io.manage_plugins import call_plugin
from ..color import rgb2gray
from .util import file_or_url_context
from ..exposure import is_low_contrast
from .._shared.utils import warn


__all__ = ['imread', 'imsave', 'imshow', 'show',
           'imread_collection', 'imshow_collection']


def imread(fname, as_gray=False, plugin=None, flatten=None,
           **plugin_args):
    """Load an image from file.

    Parameters
    ----------
    fname : string
        Image file name, e.g. ``test.jpg`` or URL.
    as_gray : bool, optional
        If True, convert color images to gray-scale (64-bit floats).
        Images that are already in gray-scale format are not converted.
    plugin : str, optional
        Name of plugin to use.  By default, the different plugins are
        tried (starting with the Python Imaging Library) until a suitable
        candidate is found.  If not given and fname is a tiff file, the
        tifffile plugin will be used.

    Other Parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.
    flatten : bool
        Backward compatible keyword, superseded by `as_gray`.

    Returns
    -------
    img_array : ndarray
        The different color bands/channels are stored in the
        third dimension, such that a gray-image is MxN, an
        RGB-image MxNx3 and an RGBA-image MxNx4.

    """
    if 'as_grey' in plugin_args.keys():
        as_gray = plugin_args.pop('as_grey', as_gray)
        warn('`as_grey` has been deprecated in favor of `as_gray`')

    # Backward compatibility
    if flatten is not None:
        as_gray = flatten
        warn('`flatten` has been deprecated in favor of `as_gray`'
             ' and will be removed in v0.16.')

    if plugin is None and hasattr(fname, 'lower'):
        if fname.lower().endswith(('.tiff', '.tif')):
            plugin = 'tifffile'

    with file_or_url_context(fname) as fname:
        img = call_plugin('imread', fname, plugin=plugin, **plugin_args)

    if not hasattr(img, 'ndim'):
        return img

    if img.ndim > 2:
        if img.shape[-1] not in (3, 4) and img.shape[-3] in (3, 4):
            img = np.swapaxes(img, -1, -3)
            img = np.swapaxes(img, -2, -3)

        if as_gray:
            img = rgb2gray(img)

    return img


def imread_collection(load_pattern, conserve_memory=True,
                      plugin=None, **plugin_args):
    """
    Load a collection of images.

    Parameters
    ----------
    load_pattern : str or list
        List of objects to load. These are usually filenames, but may
        vary depending on the currently active plugin.  See the docstring
        for ``ImageCollection`` for the default behaviour of this parameter.
    conserve_memory : bool, optional
        If True, never keep more than one in memory at a specific
        time.  Otherwise, images will be cached once they are loaded.

    Returns
    -------
    ic : ImageCollection
        Collection of images.

    Other parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    """
    return call_plugin('imread_collection', load_pattern, conserve_memory,
                       plugin=plugin, **plugin_args)


def imsave(fname, arr, plugin=None, **plugin_args):
    """Save an image to file.

    Parameters
    ----------
    fname : str
        Target filename.
    arr : ndarray of shape (M,N) or (M,N,3) or (M,N,4)
        Image data.
    plugin : str
        Name of plugin to use.  By default, the different plugins are
        tried (starting with the Python Imaging Library) until a suitable
        candidate is found.  If not given and fname is a tiff file, the
        tifffile plugin will be used.

    Other parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    Notes
    -----
    When saving a JPEG, the compression ratio may be controlled using the
    ``quality`` keyword argument which is an integer with values in [1, 100]
    where 1 is worst quality and smallest file size, and 100 is best quality
    and largest file size (default 75).  This is only available when using
    the PIL and imageio plugins.
    """
    if plugin is None and hasattr(fname, 'lower'):
        if fname.lower().endswith(('.tiff', '.tif')):
            plugin = 'tifffile'
    if is_low_contrast(arr):
        warn('%s is a low contrast image' % fname)
    if arr.dtype == bool:
        warn('%s is a boolean image: setting True to 1 and False to 0' % fname)
    return call_plugin('imsave', fname, arr, plugin=plugin, **plugin_args)


def imshow(arr, plugin=None, **plugin_args):
    """Display an image.

    Parameters
    ----------
    arr : ndarray or str
        Image data or name of image file.
    plugin : str
        Name of plugin to use.  By default, the different plugins are
        tried (starting with the Python Imaging Library) until a suitable
        candidate is found.

    Other parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    """
    if isinstance(arr, six.string_types):
        arr = call_plugin('imread', arr, plugin=plugin)
    return call_plugin('imshow', arr, plugin=plugin, **plugin_args)


def imshow_collection(ic, plugin=None, **plugin_args):
    """Display a collection of images.

    Parameters
    ----------
    ic : ImageCollection
        Collection to display.
    plugin : str
        Name of plugin to use.  By default, the different plugins are
        tried until a suitable candidate is found.

    Other parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    """
    return call_plugin('imshow_collection', ic, plugin=plugin, **plugin_args)


def show():
    '''Display pending images.

    Launch the event loop of the current gui plugin, and display all
    pending images, queued via `imshow`. This is required when using
    `imshow` from non-interactive scripts.

    A call to `show` will block execution of code until all windows
    have been closed.

    Examples
    --------
    >>> import skimage.io as io

    >>> for i in range(4):
    ...     ax_im = io.imshow(np.random.rand(50, 50))
    >>> io.show() # doctest: +SKIP

    '''
    return call_plugin('_app_show')
