# -*- coding: utf-8 -*-
# imageio is distributed under the terms of the (new) BSD License.

import re
import warnings
from numbers import Number
from pathlib import Path
from typing import Dict

import numpy as np

from imageio.core.legacy_plugin_wrapper import LegacyPlugin
from imageio.core.util import Array
from imageio.core.v3_plugin_api import PluginV3

from . import formats
from .config import known_extensions, known_plugins
from .core import RETURN_BYTES
from .core.imopen import imopen

MEMTEST_DEFAULT_MIM = "256MB"
MEMTEST_DEFAULT_MVOL = "1GB"


mem_re = re.compile(r"^(\d+\.?\d*)\s*([kKMGTPEZY]?i?)B?$")
sizes = {"": 1, None: 1}
for i, si in enumerate([""] + list("kMGTPEZY")):
    sizes[si] = 1000**i
    if si:
        sizes[si.upper() + "i"] = 1024**i


def to_nbytes(arg, default=None):
    if not arg:
        arg = float("inf")

    if arg is True:
        arg = default

    if isinstance(arg, Number):
        return arg

    match = mem_re.match(arg)
    if match is None:
        raise ValueError(
            "Memory size could not be parsed "
            "(is your capitalisation correct?): {}".format(arg)
        )

    num, unit = match.groups()

    try:
        return float(num) * sizes[unit]
    except KeyError:  # pragma: no cover
        # Note: I don't think we can reach this
        raise ValueError(
            "Memory size unit not recognised "
            "(is your capitalisation correct?): {}".format(unit)
        )


def help(name=None):
    """help(name=None)

    Print the documentation of the format specified by name, or a list
    of supported formats if name is omitted.

    Parameters
    ----------
    name : str
        Can be the name of a format, a filename extension, or a full
        filename. See also the :doc:`formats page <../formats/index>`.
    """
    if not name:
        print(formats)
    else:
        print(formats[name])


def decypher_format_arg(format_name: str) -> Dict[str, str]:
    """Split format into plugin and format

    The V2 API aliases plugins and supported formats. This function
    splits these so that they can be fed separately to `iio.imopen`.

    """

    plugin = None
    extension = None

    if format_name is None:
        pass  # nothing to do
    elif Path(format_name).suffix.lower() in known_extensions:
        extension = Path(format_name).suffix.lower()
    elif format_name in known_plugins:
        plugin = format_name
    elif format_name.upper() in known_plugins:
        plugin = format_name.upper()
    elif format_name.lower() in known_extensions:
        extension = format_name.lower()
    elif "." + format_name.lower() in known_extensions:
        extension = "." + format_name.lower()
    else:
        raise IndexError(f"No format known by name `{plugin}`.")

    return {"plugin": plugin, "extension": extension}


class LegacyReader:
    def __init__(self, plugin_instance: PluginV3, **kwargs):
        self.instance = plugin_instance
        self.last_index = 0
        self.closed = False

        if (
            type(self.instance).__name__ == "PillowPlugin"
            and kwargs.get("pilmode") is not None
        ):
            kwargs["mode"] = kwargs["pilmode"]
            del kwargs["pilmode"]

        self.read_args = kwargs

    def close(self):
        if not self.closed:
            self.instance.close()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    @property
    def request(self):
        return self.instance.request

    @property
    def format(self):
        raise TypeError("V3 Plugins don't have a format.")

    def get_length(self):
        return self.instance.properties(index=...).n_images

    def get_data(self, index):
        self.last_index = index
        img = self.instance.read(index=index, **self.read_args)
        metadata = self.instance.metadata(index=index, exclude_applied=False)
        return Array(img, metadata)

    def get_next_data(self):
        return self.get_data(self.last_index + 1)

    def set_image_index(self, index):
        self.last_index = index - 1

    def get_meta_data(self, index=None):
        return self.instance.metadata(index=index, exclude_applied=False)

    def iter_data(self):
        for idx, img in enumerate(self.instance.iter()):
            metadata = self.instance.metadata(index=idx, exclude_applied=False)
            yield Array(img, metadata)

    def __iter__(self):
        return self.iter_data()

    def __len__(self):
        return self.get_length()


class LegacyWriter:
    def __init__(self, plugin_instance: PluginV3, **kwargs):
        self.instance = plugin_instance
        self.last_index = 0
        self.closed = False

        if type(self.instance).__name__ == "PillowPlugin" and "pilmode" in kwargs:
            kwargs["mode"] = kwargs["pilmode"]
            del kwargs["pilmode"]

        self.write_args = kwargs

    def close(self):
        if not self.closed:
            self.instance.close()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    @property
    def request(self):
        return self.instance.request

    @property
    def format(self):
        raise TypeError("V3 Plugins don't have a format.")

    def append_data(self, im, meta=None):
        # TODO: write metadata in the future; there is currently no
        # generic way to do this with v3 plugins :(
        if meta is not None:
            warnings.warn(
                "V3 Plugins currently don't have a uniform way to"
                " write metadata, so any metadata is ignored."
            )

        # total_meta = dict()
        # if meta is None:
        #     meta = {}
        # if hasattr(im, "meta") and isinstance(im.meta, dict):
        #     total_meta.update(im.meta)
        # total_meta.update(meta)

        return self.instance.write(im, **self.write_args)

    def set_meta_data(self, meta):
        # TODO: write metadata
        raise NotImplementedError(
            "V3 Plugins don't have a uniform way to write metadata (yet)."
        )


def is_batch(ndimage):
    if isinstance(ndimage, (list, tuple)):
        return True

    ndimage = np.asarray(ndimage)
    if ndimage.ndim <= 2:
        return False
    elif ndimage.ndim == 3 and ndimage.shape[2] < 5:
        return False

    return True


def is_volume(ndimage):
    ndimage = np.asarray(ndimage)
    if not is_batch(ndimage):
        return False

    if ndimage.ndim == 3 and ndimage.shape[2] >= 5:
        return True
    elif ndimage.ndim == 4 and ndimage.shape[3] < 5:
        return True
    else:
        return False


# Base functions that return a reader/writer


def get_reader(uri, format=None, mode="?", **kwargs):
    """get_reader(uri, format=None, mode='?', **kwargs)

    Returns a :class:`.Reader` object which can be used to read data
    and meta data from the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    mode : {'i', 'I', 'v', 'V', '?'}
        Used to give the reader a hint on what the user expects (default "?"):
        "i" for an image, "I" for multiple images, "v" for a volume,
        "V" for multiple volumes, "?" for don't care.
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True

    image_file = imopen(uri, "r" + mode, **imopen_args)

    if isinstance(image_file, LegacyPlugin):
        return image_file.legacy_get_reader(**kwargs)
    else:
        return LegacyReader(image_file, **kwargs)


def get_writer(uri, format=None, mode="?", **kwargs):
    """get_writer(uri, format=None, mode='?', **kwargs)

    Returns a :class:`.Writer` object which can be used to write data
    and meta data to the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, file}
        The resource to write the image to, e.g. a filename, pathlib.Path
        or file object, see the docs for more info.
    format : str
        The format to use to write the file. By default imageio selects
        the appropriate for you based on the filename.
    mode : {'i', 'I', 'v', 'V', '?'}
        Used to give the writer a hint on what the user expects (default '?'):
        "i" for an image, "I" for multiple images, "v" for a volume,
        "V" for multiple volumes, "?" for don't care.
    kwargs : ...
        Further keyword arguments are passed to the writer. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True

    image_file = imopen(uri, "w" + mode, **imopen_args)
    if isinstance(image_file, LegacyPlugin):
        return image_file.legacy_get_writer(**kwargs)
    else:
        return LegacyWriter(image_file, **kwargs)


# Images


def imread(uri, format=None, **kwargs):
    """imread(uri, format=None, **kwargs)

    Reads an image from the specified file. Returns a numpy array, which
    comes with a dict of meta data at its 'meta' attribute.

    Note that the image data is returned as-is, and may not always have
    a dtype of uint8 (and thus may differ from what e.g. PIL returns).

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True

    with imopen(uri, "ri", **imopen_args) as file:
        result = file.read(index=0, **kwargs)

    return result


def imwrite(uri, im, format=None, **kwargs):
    """imwrite(uri, im, format=None, **kwargs)

    Write an image to the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, file}
        The resource to write the image to, e.g. a filename, pathlib.Path
        or file object, see the docs for more info.
    im : numpy.ndarray
        The image data. Must be NxM, NxMx3 or NxMx4.
    format : str
        The format to use to write the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the writer. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    # Test image
    imt = type(im)
    im = np.asarray(im)
    if not np.issubdtype(im.dtype, np.number):
        raise ValueError("Image is not numeric, but {}.".format(imt.__name__))

    if is_batch(im) or im.ndim < 2:
        raise ValueError("Image must be 2D (grayscale, RGB, or RGBA).")

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "wi", **imopen_args) as file:
        return file.write(im, **kwargs)


# Multiple images


def mimread(uri, format=None, memtest=MEMTEST_DEFAULT_MIM, **kwargs):
    """mimread(uri, format=None, memtest="256MB", **kwargs)

    Reads multiple images from the specified file. Returns a list of
    numpy arrays, each with a dict of meta data at its 'meta' attribute.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the images from, e.g. a filename,pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    memtest : {bool, int, float, str}
        If truthy, this function will raise an error if the resulting
        list of images consumes greater than the amount of memory specified.
        This is to protect the system from using so much memory that it needs
        to resort to swapping, and thereby stall the computer. E.g.
        ``mimread('hunger_games.avi')``.

        If the argument is a number, that will be used as the threshold number
        of bytes.

        If the argument is a string, it will be interpreted as a number of bytes with
        SI/IEC prefixed units (e.g. '1kB', '250MiB', '80.3YB').

        - Units are case sensitive
        - k, M etc. represent a 1000-fold change, where Ki, Mi etc. represent 1024-fold
        - The "B" is optional, but if present, must be capitalised

        If the argument is True, the default will be used, for compatibility reasons.

        Default: '256MB'
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    # used for mimread and mvolread
    nbyte_limit = to_nbytes(memtest, MEMTEST_DEFAULT_MIM)

    images = list()
    nbytes = 0

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "rI", **imopen_args) as file:
        for image in file.iter(**kwargs):
            images.append(image)
            nbytes += image.nbytes
            if nbytes > nbyte_limit:
                raise RuntimeError(
                    "imageio.mimread() has read over {}B of "
                    "image data.\nStopped to avoid memory problems."
                    " Use imageio.get_reader(), increase threshold, or memtest=False".format(
                        int(nbyte_limit)
                    )
                )

    if len(images) == 1 and is_batch(images[0]):
        images = [*images[0]]

    return images


def mimwrite(uri, ims, format=None, **kwargs):
    """mimwrite(uri, ims, format=None, **kwargs)

    Write multiple images to the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, file}
        The resource to write the images to, e.g. a filename, pathlib.Path
        or file object, see the docs for more info.
    ims : sequence of numpy arrays
        The image data. Each array must be NxM, NxMx3 or NxMx4.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the writer. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    if not is_batch(ims):
        raise ValueError("Image data must be a sequence of ndimages.")

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "wI", **imopen_args) as file:
        return file.write(ims, is_batch=True, **kwargs)


# Volumes


def volread(uri, format=None, **kwargs):
    """volread(uri, format=None, **kwargs)

    Reads a volume from the specified file. Returns a numpy array, which
    comes with a dict of meta data at its 'meta' attribute.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the volume from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "rv", **imopen_args) as file:
        return file.read(index=0, **kwargs)


def volwrite(uri, im, format=None, **kwargs):
    """volwrite(uri, vol, format=None, **kwargs)

    Write a volume to the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, file}
        The resource to write the image to, e.g. a filename, pathlib.Path
        or file object, see the docs for more info.
    vol : numpy.ndarray
        The image data. Must be NxMxL (or NxMxLxK if each voxel is a tuple).
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the writer. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    # Test image
    im = np.asarray(im)
    if not is_volume(im):
        raise ValueError("Image must be 3D, or 4D if each voxel is a tuple.")

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True

    with imopen(uri, "wv", **imopen_args) as file:
        return file.write(im, is_batch=False, **kwargs)


# Multiple volumes


def mvolread(uri, format=None, memtest=MEMTEST_DEFAULT_MVOL, **kwargs):
    """mvolread(uri, format=None, memtest='1GB', **kwargs)

    Reads multiple volumes from the specified file. Returns a list of
    numpy arrays, each with a dict of meta data at its 'meta' attribute.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the volumes from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    memtest : {bool, int, float, str}
        If truthy, this function will raise an error if the resulting
        list of images consumes greater than the amount of memory specified.
        This is to protect the system from using so much memory that it needs
        to resort to swapping, and thereby stall the computer. E.g.
        ``mimread('hunger_games.avi')``.

        If the argument is a number, that will be used as the threshold number
        of bytes.

        If the argument is a string, it will be interpreted as a number of bytes with
        SI/IEC prefixed units (e.g. '1kB', '250MiB', '80.3YB').

        - Units are case sensitive
        - k, M etc. represent a 1000-fold change, where Ki, Mi etc. represent 1024-fold
        - The "B" is optional, but if present, must be capitalised

        If the argument is True, the default will be used, for compatibility reasons.

        Default: '1GB'
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    # used for mimread and mvolread
    nbyte_limit = to_nbytes(memtest, MEMTEST_DEFAULT_MVOL)

    images = list()
    nbytes = 0
    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "rV", **imopen_args) as file:
        for image in file.iter(**kwargs):
            images.append(image)
            nbytes += image.nbytes
            if nbytes > nbyte_limit:
                raise RuntimeError(
                    "imageio.mimread() has read over {}B of "
                    "image data.\nStopped to avoid memory problems."
                    " Use imageio.get_reader(), increase threshold, or memtest=False".format(
                        int(nbyte_limit)
                    )
                )

    return images


def mvolwrite(uri, ims, format=None, **kwargs):
    """mvolwrite(uri, vols, format=None, **kwargs)

    Write multiple volumes to the specified file.

    Parameters
    ----------
    uri : {str, pathlib.Path, file}
        The resource to write the volumes to, e.g. a filename, pathlib.Path
        or file object, see the docs for more info.
    ims : sequence of numpy arrays
        The image data. Each array must be NxMxL (or NxMxLxK if each
        voxel is a tuple).
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the writer. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    for im in ims:
        if not is_volume(im):
            raise ValueError("Image must be 3D, or 4D if each voxel is a tuple.")

    imopen_args = decypher_format_arg(format)
    imopen_args["legacy_mode"] = True
    with imopen(uri, "wV", **imopen_args) as file:
        return file.write(ims, is_batch=True, **kwargs)


# aliases
read = get_reader
save = get_writer
imsave = imwrite
mimsave = mimwrite
volsave = volwrite
mvolsave = mvolwrite

__all__ = [
    "imread",
    "mimread",
    "volread",
    "mvolread",
    "imwrite",
    "mimwrite",
    "volwrite",
    "mvolwrite",
    # misc
    "help",
    "get_reader",
    "get_writer",
    "RETURN_BYTES",
]
