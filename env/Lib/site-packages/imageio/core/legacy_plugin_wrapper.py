from pathlib import Path

import numpy as np

from ..config import known_extensions
from .request import InitializationError, IOMode
from .v3_plugin_api import ImageProperties, PluginV3


def _legacy_default_index(format):
    if format._name == "FFMPEG":
        index = Ellipsis
    elif format._name == "GIF-PIL":
        index = Ellipsis
    else:
        index = 0

    return index


class LegacyPlugin(PluginV3):
    """A plugin to  make old (v2.9) plugins compatible with v3.0

    .. depreciated:: 2.9
        `legacy_get_reader` will be removed in a future version of imageio.
        `legacy_get_writer` will be removed in a future version of imageio.

    This plugin is a wrapper around the old FormatManager class and exposes
    all the old plugins via the new API. On top of this it has
    ``legacy_get_reader`` and ``legacy_get_writer`` methods to allow using
    it with the v2.9 API.

    Methods
    -------
    read(index=None, **kwargs)
        Read the image at position ``index``.
    write(image, **kwargs)
        Write image to the URI.
    iter(**kwargs)
        Iteratively yield images from the given URI.
    get_meta(index=None)
        Return the metadata for the image at position ``index``.
    legacy_get_reader(**kwargs)
        Returns the v2.9 image reader. (depreciated)
    legacy_get_writer(**kwargs)
        Returns the v2.9 image writer. (depreciated)

    Examples
    --------

    >>> import imageio.v3 as iio
    >>> with iio.imopen("/path/to/image.tiff", "r", legacy_mode=True) as file:
    >>>     reader = file.legacy_get_reader()  # depreciated
    >>>     for im in file.iter():
    >>>         print(im.shape)

    """

    def __init__(self, request, legacy_plugin):
        """Instantiate a new Legacy Plugin

        Parameters
        ----------
        uri : {str, pathlib.Path, bytes, file}
            The resource to load the image from, e.g. a filename, pathlib.Path,
            http address or file object, see the docs for more info.
        legacy_plugin : Format
            The (legacy) format to use to interface with the URI.

        """
        self._request = request
        self._format = legacy_plugin

        source = (
            "<bytes>"
            if isinstance(self._request.raw_uri, bytes)
            else self._request.raw_uri
        )
        if self._request.mode.io_mode == IOMode.read:
            if not self._format.can_read(request):
                raise InitializationError(
                    f"`{self._format.name}`" f" can not read `{source}`."
                )
        else:
            if not self._format.can_write(request):
                raise InitializationError(
                    f"`{self._format.name}`" f" can not write to `{source}`."
                )

    def legacy_get_reader(self, **kwargs):
        """legacy_get_reader(**kwargs)

        a utility method to provide support vor the V2.9 API

        Parameters
        ----------
        kwargs : ...
            Further keyword arguments are passed to the reader. See :func:`.help`
            to see what arguments are available for a particular format.
        """

        # Note: this will break thread-safety
        self._request._kwargs = kwargs

        # safeguard for DICOM plugin reading from folders
        try:
            assert Path(self._request.filename).is_dir()
        except OSError:
            pass  # not a valid path on this OS
        except AssertionError:
            pass  # not a folder
        else:
            return self._format.get_reader(self._request)

        self._request.get_file().seek(0)
        return self._format.get_reader(self._request)

    def read(self, *, index=None, **kwargs):
        """
        Parses the given URI and creates a ndarray from it.

        Parameters
        ----------
        index : {integer, None}
            If the URI contains a list of ndimages return the index-th
            image. If None, stack all images into an ndimage along the
            0-th dimension (equivalent to np.stack(imgs, axis=0)).
        kwargs : ...
            Further keyword arguments are passed to the reader. See
            :func:`.help` to see what arguments are available for a particular
            format.

        Returns
        -------
        ndimage : np.ndarray
            A numpy array containing the decoded image data.

        """

        if index is None:
            index = _legacy_default_index(self._format)

        if index is Ellipsis:
            img = np.stack([im for im in self.iter(**kwargs)])
            return img

        reader = self.legacy_get_reader(**kwargs)
        return reader.get_data(index)

    def legacy_get_writer(self, **kwargs):
        """legacy_get_writer(**kwargs)

        Returns a :class:`.Writer` object which can be used to write data
        and meta data to the specified file.

        Parameters
        ----------
        kwargs : ...
            Further keyword arguments are passed to the writer. See :func:`.help`
            to see what arguments are available for a particular format.
        """

        # Note: this will break thread-safety
        self._request._kwargs = kwargs
        return self._format.get_writer(self._request)

    def write(self, ndimage, *, is_batch=None, metadata=None, **kwargs):
        """
        Write an ndimage to the URI specified in path.

        If the URI points to a file on the current host and the file does not
        yet exist it will be created. If the file exists already, it will be
        appended if possible; otherwise, it will be replaced.

        Parameters
        ----------
        ndimage : numpy.ndarray
            The ndimage or list of ndimages to write.
        is_batch : bool
            If True, treat the supplied ndimage as a batch of images. If False,
            treat the supplied ndimage as a single image. If None, try to
            determine ``is_batch`` from the ndimage's shape and ndim.
        metadata : dict
            The metadata passed to write alongside the image.
        kwargs : ...
            Further keyword arguments are passed to the writer. See
            :func:`.help` to see what arguments are available for a
            particular format.


        Returns
        -------
        buffer : bytes
            When writing to the special target "<bytes>", this function will
            return the encoded image data as a bytes string. Otherwise it
            returns None.

        Notes
        -----
        Automatically determining ``is_batch`` may fail for some images due to
        shape aliasing. For example, it may classify a channel-first color image
        as a batch of gray images. In most cases this automatic deduction works
        fine (it has for almost a decade), but if you do have one of those edge
        cases (or are worried that you might) consider explicitly setting
        ``is_batch``.

        """

        if is_batch or isinstance(ndimage, (list, tuple)):
            pass  # ndimage is list of images
        elif is_batch is False:
            ndimage = [ndimage]
        else:
            # Write the largest possible block by guessing the meaning of each
            # dimension from the shape/ndim and then checking if any batch
            # dimensions are left.
            ndimage = np.asanyarray(ndimage)
            batch_dims = ndimage.ndim

            # two spatial dimensions
            batch_dims = max(batch_dims - 2, 0)

            # packed (channel-last) image
            if ndimage.ndim >= 3 and ndimage.shape[-1] < 5:
                batch_dims = max(batch_dims - 1, 0)

            # format supports volumetric images
            ext_infos = known_extensions.get(self._request.extension, list())
            for ext_info in ext_infos:
                if self._format.name in ext_info.priority and ext_info.volume_support:
                    batch_dims = max(batch_dims - 1, 0)
                    break

            if batch_dims == 0:
                ndimage = [ndimage]

        with self.legacy_get_writer(**kwargs) as writer:
            for image in ndimage:
                image = np.asanyarray(image)

                if image.ndim < 2:
                    raise ValueError(
                        "The image must have at least two spatial dimensions."
                    )

                if not np.issubdtype(image.dtype, np.number) and not np.issubdtype(
                    image.dtype, bool
                ):
                    raise ValueError(
                        f"All images have to be numeric, and not `{image.dtype}`."
                    )

                writer.append_data(image, metadata)

        return writer.request.get_result()

    def iter(self, **kwargs):
        """Iterate over a list of ndimages given by the URI

        Parameters
        ----------
        kwargs : ...
            Further keyword arguments are passed to the reader. See
            :func:`.help` to see what arguments are available for a particular
            format.
        """

        reader = self.legacy_get_reader(**kwargs)
        for image in reader:
            yield image

    def properties(self, index=None):
        """Standardized ndimage metadata.

        Parameters
        ----------
        index : int
            The index of the ndimage for which to return properties. If the
            index is out of bounds a ``ValueError`` is raised. If ``None``,
            return the properties for the ndimage stack. If this is impossible,
            e.g., due to shape mismatch, an exception will be raised.

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.

        """

        if index is None:
            index = _legacy_default_index(self._format)

        # for backwards compatibility ... actually reads pixel data :(
        if index is Ellipsis:
            image = self.read(index=0)
            n_images = self.legacy_get_reader().get_length()
            return ImageProperties(
                shape=(n_images, *image.shape),
                dtype=image.dtype,
                n_images=n_images,
                is_batch=True,
            )

        image = self.read(index=index)
        return ImageProperties(
            shape=image.shape,
            dtype=image.dtype,
            is_batch=False,
        )

    def get_meta(self, *, index=None):
        """Read ndimage metadata from the URI

        Parameters
        ----------
        index : {integer, None}
            If the URI contains a list of ndimages return the metadata
            corresponding to the index-th image. If None, behavior depends on
            the used api

            Legacy-style API: return metadata of the first element (index=0)
            New-style API: Behavior depends on the used Plugin.

        Returns
        -------
        metadata : dict
            A dictionary of metadata.

        """

        return self.metadata(index=index, exclude_applied=False)

    def metadata(self, index=None, exclude_applied: bool = True):
        """Format-Specific ndimage metadata.

        Parameters
        ----------
        index : int
            The index of the ndimage to read. If the index is out of bounds a
            ``ValueError`` is raised. If ``None``, global metadata is returned.
        exclude_applied : bool
            This parameter exists for compatibility and has no effect. Legacy
            plugins always report all metadata they find.

        Returns
        -------
        metadata : dict
            A dictionary filled with format-specific metadata fields and their
            values.

        """

        if index is None:
            index = _legacy_default_index(self._format)

        return self.legacy_get_reader().get_meta_data(index=index)

    def __del__(self) -> None:
        pass
        # turns out we can't close the file here for LegacyPlugin
        # because it would break backwards compatibility
        # with legacy_get_writer and legacy_get_reader
        # self._request.finish()
