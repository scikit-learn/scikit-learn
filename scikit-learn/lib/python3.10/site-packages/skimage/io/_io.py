import pathlib
import warnings

import numpy as np

from .._shared.utils import warn, deprecate_func, deprecate_parameter, DEPRECATED
from .._shared.version_requirements import require
from ..exposure import is_low_contrast
from ..color.colorconv import rgb2gray, rgba2rgb
from ..io.manage_plugins import call_plugin, _hide_plugin_deprecation_warnings
from .util import file_or_url_context

__all__ = [
    'imread',
    'imsave',
    'imshow',
    'show',
    'imread_collection',
    'imshow_collection',
]


_remove_plugin_param_template = (
    "The plugin infrastructure in `skimage.io` and the parameter "
    "`{deprecated_name}` are deprecated since version {deprecated_version} and "
    "will be removed in {changed_version} (or later). To avoid this warning, "
    "please do not use the parameter `{deprecated_name}`. Instead, use `imageio` "
    "or other I/O packages directly. See also `{func_name}`."
)


@deprecate_parameter(
    "plugin",
    start_version="0.25",
    stop_version="0.27",
    template=_remove_plugin_param_template,
)
def imread(fname, as_gray=False, plugin=DEPRECATED, **plugin_args):
    """Load an image from file.

    Parameters
    ----------
    fname : str or pathlib.Path
        Image file name, e.g. ``test.jpg`` or URL.
    as_gray : bool, optional
        If True, convert color images to gray-scale (64-bit floats).
        Images that are already in gray-scale format are not converted.

    Other Parameters
    ----------------
    plugin_args : DEPRECATED
        The plugin infrastructure is deprecated.

    Returns
    -------
    img_array : ndarray
        The different color bands/channels are stored in the
        third dimension, such that a gray-image is MxN, an
        RGB-image MxNx3 and an RGBA-image MxNx4.

    """
    if plugin is DEPRECATED:
        plugin = None
    if plugin_args:
        msg = (
            "The plugin infrastructure in `skimage.io` is deprecated since "
            "version 0.25 and will be removed in 0.27 (or later). To avoid "
            "this warning, please do not pass additional keyword arguments "
            "for plugins (`**plugin_args`). Instead, use `imageio` or other "
            "I/O packages directly. See also `skimage.io.imread`."
        )
        warnings.warn(msg, category=FutureWarning, stacklevel=3)

    if isinstance(fname, pathlib.Path):
        fname = str(fname.resolve())

    if plugin is None and hasattr(fname, 'lower'):
        if fname.lower().endswith(('.tiff', '.tif')):
            plugin = 'tifffile'

    with file_or_url_context(fname) as fname, _hide_plugin_deprecation_warnings():
        img = call_plugin('imread', fname, plugin=plugin, **plugin_args)

    if not hasattr(img, 'ndim'):
        return img

    if img.ndim > 2:
        if img.shape[-1] not in (3, 4) and img.shape[-3] in (3, 4):
            img = np.swapaxes(img, -1, -3)
            img = np.swapaxes(img, -2, -3)

        if as_gray:
            if img.shape[2] == 4:
                img = rgba2rgb(img)
            img = rgb2gray(img)

    return img


@deprecate_parameter(
    "plugin",
    start_version="0.25",
    stop_version="0.27",
    template=_remove_plugin_param_template,
)
def imread_collection(
    load_pattern, conserve_memory=True, plugin=DEPRECATED, **plugin_args
):
    """
    Load a collection of images.

    Parameters
    ----------
    load_pattern : str or list
        List of objects to load. These are usually filenames, but may
        vary depending on the currently active plugin. See :class:`ImageCollection`
        for the default behaviour of this parameter.
    conserve_memory : bool, optional
        If True, never keep more than one in memory at a specific
        time.  Otherwise, images will be cached once they are loaded.

    Returns
    -------
    ic : :class:`ImageCollection`
        Collection of images.

    Other Parameters
    ----------------
    plugin_args : DEPRECATED
        The plugin infrastructure is deprecated.

    """
    if plugin is DEPRECATED:
        plugin = None
    if plugin_args:
        msg = (
            "The plugin infrastructure in `skimage.io` is deprecated since "
            "version 0.25 and will be removed in 0.27 (or later). To avoid "
            "this warning, please do not pass additional keyword arguments "
            "for plugins (`**plugin_args`). Instead, use `imageio` or other "
            "I/O packages directly. See also `skimage.io.imread_collection`."
        )
        warnings.warn(msg, category=FutureWarning, stacklevel=3)
    with _hide_plugin_deprecation_warnings():
        return call_plugin(
            'imread_collection',
            load_pattern,
            conserve_memory,
            plugin=plugin,
            **plugin_args,
        )


@deprecate_parameter(
    "plugin",
    start_version="0.25",
    stop_version="0.27",
    template=_remove_plugin_param_template,
)
def imsave(fname, arr, plugin=DEPRECATED, *, check_contrast=True, **plugin_args):
    """Save an image to file.

    Parameters
    ----------
    fname : str or pathlib.Path
        Target filename.
    arr : ndarray of shape (M,N) or (M,N,3) or (M,N,4)
        Image data.
    check_contrast : bool, optional
        Check for low contrast and print warning (default: True).

    Other Parameters
    ----------------
    plugin_args : DEPRECATED
        The plugin infrastructure is deprecated.
    """
    if plugin is DEPRECATED:
        plugin = None
    if plugin_args:
        msg = (
            "The plugin infrastructure in `skimage.io` is deprecated since "
            "version 0.25 and will be removed in 0.27 (or later). To avoid "
            "this warning, please do not pass additional keyword arguments "
            "for plugins (`**plugin_args`). Instead, use `imageio` or other "
            "I/O packages directly. See also `skimage.io.imsave`."
        )
        warnings.warn(msg, category=FutureWarning, stacklevel=3)

    if isinstance(fname, pathlib.Path):
        fname = str(fname.resolve())
    if plugin is None and hasattr(fname, 'lower'):
        if fname.lower().endswith(('.tiff', '.tif')):
            plugin = 'tifffile'
    if arr.dtype == bool:
        warn(
            f'{fname} is a boolean image: setting True to 255 and False to 0. '
            'To silence this warning, please convert the image using '
            'img_as_ubyte.',
            stacklevel=3,
        )
        arr = arr.astype('uint8') * 255
    if check_contrast and is_low_contrast(arr):
        warn(f'{fname} is a low contrast image')

    with _hide_plugin_deprecation_warnings():
        return call_plugin('imsave', fname, arr, plugin=plugin, **plugin_args)


@deprecate_func(
    deprecated_version="0.25",
    removed_version="0.27",
    hint="Please use `matplotlib`, `napari`, etc. to visualize images.",
)
def imshow(arr, plugin=None, **plugin_args):
    """Display an image.

    Parameters
    ----------
    arr : ndarray or str
        Image data or name of image file.
    plugin : str
        Name of plugin to use.  By default, the different plugins are
        tried (starting with imageio) until a suitable candidate is found.

    Other Parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    """
    if isinstance(arr, str):
        arr = call_plugin('imread', arr, plugin=plugin)
    with _hide_plugin_deprecation_warnings():
        return call_plugin('imshow', arr, plugin=plugin, **plugin_args)


@deprecate_func(
    deprecated_version="0.25",
    removed_version="0.27",
    hint="Please use `matplotlib`, `napari`, etc. to visualize images.",
)
def imshow_collection(ic, plugin=None, **plugin_args):
    """Display a collection of images.

    Parameters
    ----------
    ic : :class:`ImageCollection`
        Collection to display.

    Other Parameters
    ----------------
    plugin_args : keywords
        Passed to the given plugin.

    """
    with _hide_plugin_deprecation_warnings():
        return call_plugin('imshow_collection', ic, plugin=plugin, **plugin_args)


@require("matplotlib", ">=3.3")
@deprecate_func(
    deprecated_version="0.25",
    removed_version="0.27",
    hint="Please use `matplotlib`, `napari`, etc. to visualize images.",
)
def show():
    """Display pending images.

    Launch the event loop of the current GUI plugin, and display all
    pending images, queued via `imshow`. This is required when using
    `imshow` from non-interactive scripts.

    A call to `show` will block execution of code until all windows
    have been closed.

    Examples
    --------
    .. testsetup::
        >>> import pytest; _ = pytest.importorskip('matplotlib')

    >>> import skimage.io as io
    >>> rng = np.random.default_rng()
    >>> for i in range(4):
    ...     ax_im = io.imshow(rng.random((50, 50)))  # doctest: +SKIP
    >>> io.show() # doctest: +SKIP

    """
    with _hide_plugin_deprecation_warnings():
        return call_plugin('_app_show')
