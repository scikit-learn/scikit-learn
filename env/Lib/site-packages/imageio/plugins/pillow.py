# -*- coding: utf-8 -*-
# imageio is distributed under the terms of the (new) BSD License.

""" Read/Write images using Pillow/PIL.

Backend Library: `Pillow <https://pillow.readthedocs.io/en/stable/>`_

Plugin that wraps the the Pillow library. Pillow is a friendly fork of PIL
(Python Image Library) and supports reading and writing of common formats (jpg,
png, gif, tiff, ...). For, the complete list of features and supported formats
please refer to pillows official docs (see the Backend Library link).

Parameters
----------
request : Request
    A request object representing the resource to be operated on.

Methods
-------

.. autosummary::
    :toctree: _plugins/pillow

    PillowPlugin.read
    PillowPlugin.write
    PillowPlugin.iter
    PillowPlugin.get_meta

"""

import sys
import warnings
from io import BytesIO
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple, Union, cast

import numpy as np
from PIL import ExifTags, GifImagePlugin, Image, ImageSequence, UnidentifiedImageError
from PIL import __version__ as pil_version  # type: ignore

from ..core.request import URI_BYTES, InitializationError, IOMode, Request
from ..core.v3_plugin_api import ImageProperties, PluginV3
from ..typing import ArrayLike


def pillow_version() -> Tuple[int]:
    return tuple(int(x) for x in pil_version.split("."))


def _exif_orientation_transform(orientation: int, mode: str) -> Callable:
    # get transformation that transforms an image from a
    # given EXIF orientation into the standard orientation

    # -1 if the mode has color channel, 0 otherwise
    axis = -2 if Image.getmodebands(mode) > 1 else -1

    EXIF_ORIENTATION = {
        1: lambda x: x,
        2: lambda x: np.flip(x, axis=axis),
        3: lambda x: np.rot90(x, k=2),
        4: lambda x: np.flip(x, axis=axis - 1),
        5: lambda x: np.flip(np.rot90(x, k=3), axis=axis),
        6: lambda x: np.rot90(x, k=3),
        7: lambda x: np.flip(np.rot90(x, k=1), axis=axis),
        8: lambda x: np.rot90(x, k=1),
    }

    return EXIF_ORIENTATION[orientation]


class PillowPlugin(PluginV3):
    def __init__(self, request: Request) -> None:
        """Instantiate a new Pillow Plugin Object

        Parameters
        ----------
        request : {Request}
            A request object representing the resource to be operated on.

        """

        super().__init__(request)

        # Register HEIF opener for Pillow
        try:
            from pillow_heif import register_heif_opener
        except ImportError:
            pass
        else:
            register_heif_opener()

        # Register AVIF opener for Pillow
        try:
            from pillow_heif import register_avif_opener
        except ImportError:
            pass
        else:
            register_avif_opener()

        self._image: Image = None
        self.images_to_write = []

        if request.mode.io_mode == IOMode.read:
            try:
                with Image.open(request.get_file()):
                    # Check if it is generally possible to read the image.
                    # This will not read any data and merely try to find a
                    # compatible pillow plugin (ref: the pillow docs).
                    pass
            except UnidentifiedImageError:
                if request._uri_type == URI_BYTES:
                    raise InitializationError(
                        "Pillow can not read the provided bytes."
                    ) from None
                else:
                    raise InitializationError(
                        f"Pillow can not read {request.raw_uri}."
                    ) from None

            self._image = Image.open(self._request.get_file())
        else:
            self.save_args = {}

            extension = self.request.extension or self.request.format_hint
            if extension is None:
                warnings.warn(
                    "Can't determine file format to write as. You _must_"
                    " set `format` during write or the call will fail. Use "
                    "`extension` to supress this warning. ",
                    UserWarning,
                )
                return

            tirage = [Image.preinit, Image.init]
            for format_loader in tirage:
                format_loader()
                if extension in Image.registered_extensions().keys():
                    return

            raise InitializationError(
                f"Pillow can not write `{extension}` files."
            ) from None

    def close(self) -> None:
        self._flush_writer()

        if self._image:
            self._image.close()

        self._request.finish()

    def read(
        self,
        *,
        index: int = None,
        mode: str = None,
        rotate: bool = False,
        apply_gamma: bool = False,
        writeable_output: bool = True,
        pilmode: str = None,
        exifrotate: bool = None,
        as_gray: bool = None,
    ) -> np.ndarray:
        """
        Parses the given URI and creates a ndarray from it.

        Parameters
        ----------
        index : int
            If the ImageResource contains multiple ndimages, and index is an
            integer, select the index-th ndimage from among them and return it.
            If index is an ellipsis (...), read all ndimages in the file and
            stack them along a new batch dimension and return them. If index is
            None, this plugin reads the first image of the file (index=0) unless
            the image is a GIF or APNG, in which case all images are read
            (index=...).
        mode : str
            Convert the image to the given mode before returning it. If None,
            the mode will be left unchanged. Possible modes can be found at:
            https://pillow.readthedocs.io/en/stable/handbook/concepts.html#modes
        rotate : bool
            If True and the image contains an EXIF orientation tag,
            apply the orientation before returning the ndimage.
        apply_gamma : bool
            If True and the image contains metadata about gamma, apply gamma
            correction to the image.
        writable_output : bool
            If True, ensure that the image is writable before returning it to
            the user. This incurs a full copy of the pixel data if the data
            served by pillow is read-only. Consequentially, setting this flag to
            False improves performance for some images.
        pilmode : str
            Deprecated, use `mode` instead.
        exifrotate : bool
            Deprecated, use `rotate` instead.
        as_gray : bool
            Deprecated. Exists to raise a constructive error message.

        Returns
        -------
        ndimage : ndarray
            A numpy array containing the loaded image data

        Notes
        -----
        If you read a paletted image (e.g. GIF) then the plugin will apply the
        palette by default. Should you wish to read the palette indices of each
        pixel use ``mode="P"``. The coresponding color pallete can be found in
        the image's metadata using the ``palette`` key when metadata is
        extracted using the ``exclude_applied=False`` kwarg. The latter is
        needed, as palettes are applied by default and hence excluded by default
        to keep metadata and pixel data consistent.

        """

        if pilmode is not None:
            warnings.warn(
                "`pilmode` is deprecated. Use `mode` instead.", DeprecationWarning
            )
            mode = pilmode

        if exifrotate is not None:
            warnings.warn(
                "`exifrotate` is deprecated. Use `rotate` instead.", DeprecationWarning
            )
            rotate = exifrotate

        if as_gray is not None:
            raise TypeError(
                "The keyword `as_gray` is no longer supported."
                "Use `mode='F'` for a backward-compatible result, or "
                " `mode='L'` for an integer-valued result."
            )

        if self._image.format == "GIF":
            # Converting GIF P frames to RGB
            # https://github.com/python-pillow/Pillow/pull/6150
            GifImagePlugin.LOADING_STRATEGY = (
                GifImagePlugin.LoadingStrategy.RGB_AFTER_DIFFERENT_PALETTE_ONLY
            )

        if index is None:
            if self._image.format == "GIF":
                index = Ellipsis
            elif self._image.custom_mimetype == "image/apng":
                index = Ellipsis
            else:
                index = 0

        if isinstance(index, int):
            # will raise IO error if index >= number of frames in image
            self._image.seek(index)
            image = self._apply_transforms(
                self._image, mode, rotate, apply_gamma, writeable_output
            )
        else:
            iterator = self.iter(
                mode=mode,
                rotate=rotate,
                apply_gamma=apply_gamma,
                writeable_output=writeable_output,
            )
            image = np.stack([im for im in iterator], axis=0)

        return image

    def iter(
        self,
        *,
        mode: str = None,
        rotate: bool = False,
        apply_gamma: bool = False,
        writeable_output: bool = True,
    ) -> Iterator[np.ndarray]:
        """
        Iterate over all ndimages/frames in the URI

        Parameters
        ----------
        mode : {str, None}
            Convert the image to the given mode before returning it. If None,
            the mode will be left unchanged. Possible modes can be found at:
            https://pillow.readthedocs.io/en/stable/handbook/concepts.html#modes
        rotate : {bool}
            If set to ``True`` and the image contains an EXIF orientation tag,
            apply the orientation before returning the ndimage.
        apply_gamma : {bool}
            If ``True`` and the image contains metadata about gamma, apply gamma
            correction to the image.
        writable_output : bool
            If True, ensure that the image is writable before returning it to
            the user. This incurs a full copy of the pixel data if the data
            served by pillow is read-only. Consequentially, setting this flag to
            False improves performance for some images.
        """

        for im in ImageSequence.Iterator(self._image):
            yield self._apply_transforms(
                im, mode, rotate, apply_gamma, writeable_output
            )

    def _apply_transforms(
        self, image, mode, rotate, apply_gamma, writeable_output
    ) -> np.ndarray:
        if mode is not None:
            image = image.convert(mode)
        elif image.mode == "P":
            # adjust for pillow9 changes
            # see: https://github.com/python-pillow/Pillow/issues/5929
            image = image.convert(image.palette.mode)
        elif image.format == "PNG" and image.mode == "I":
            major, minor, patch = pillow_version()

            if sys.byteorder == "little":
                desired_mode = "I;16"
            else:  # pragma: no cover
                # can't test big-endian in GH-Actions
                desired_mode = "I;16B"

            if major < 10:  # pragma: no cover
                warnings.warn(
                    "Loading 16-bit (uint16) PNG as int32 due to limitations "
                    "in pillow's PNG decoder. This will be fixed in a future "
                    "version of pillow which will make this warning dissapear.",
                    UserWarning,
                )
            elif minor < 1:  # pragma: no cover
                # pillow<10.1.0 can directly decode into 16-bit grayscale
                image.mode = desired_mode
            else:
                # pillow >= 10.1.0
                image = image.convert(desired_mode)

        image = np.asarray(image)

        meta = self.metadata(index=self._image.tell(), exclude_applied=False)
        if rotate and "Orientation" in meta:
            transformation = _exif_orientation_transform(
                meta["Orientation"], self._image.mode
            )
            image = transformation(image)

        if apply_gamma and "gamma" in meta:
            gamma = float(meta["gamma"])
            scale = float(65536 if image.dtype == np.uint16 else 255)
            gain = 1.0
            image = ((image / scale) ** gamma) * scale * gain + 0.4999
            image = np.round(image).astype(np.uint8)

        if writeable_output and not image.flags["WRITEABLE"]:
            image = np.array(image)

        return image

    def write(
        self,
        ndimage: Union[ArrayLike, List[ArrayLike]],
        *,
        mode: str = None,
        format: str = None,
        is_batch: bool = None,
        **kwargs,
    ) -> Optional[bytes]:
        """
        Write an ndimage to the URI specified in path.

        If the URI points to a file on the current host and the file does not
        yet exist it will be created. If the file exists already, it will be
        appended if possible; otherwise, it will be replaced.

        If necessary, the image is broken down along the leading dimension to
        fit into individual frames of the chosen format. If the format doesn't
        support multiple frames, and IOError is raised.

        Parameters
        ----------
        image : ndarray or list
            The ndimage to write. If a list is given each element is expected to
            be an ndimage.
        mode : str
            Specify the image's color format. If None (default), the mode is
            inferred from the array's shape and dtype. Possible modes can be
            found at:
            https://pillow.readthedocs.io/en/stable/handbook/concepts.html#modes
        format : str
            Optional format override. If omitted, the format to use is
            determined from the filename extension. If a file object was used
            instead of a filename, this parameter must always be used.
        is_batch : bool
            Explicitly tell the writer that ``image`` is a batch of images
            (True) or not (False). If None, the writer will guess this from the
            provided ``mode`` or ``image.shape``. While the latter often works,
            it may cause problems for small images due to aliasing of spatial
            and color-channel axes.
        kwargs : ...
            Extra arguments to pass to pillow. If a writer doesn't recognise an
            option, it is silently ignored. The available options are described
            in pillow's `image format documentation
            <https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html>`_
            for each writer.

        Notes
        -----
        When writing batches of very narrow (2-4 pixels wide) gray images set
        the ``mode`` explicitly to avoid the batch being identified as a colored
        image.

        """
        if "fps" in kwargs:
            warnings.warn(
                "The keyword `fps` is no longer supported. Use `duration`"
                "(in ms) instead, e.g. `fps=50` == `duration=20` (1000 * 1/50).",
                DeprecationWarning,
            )
            kwargs["duration"] = 1000 * 1 / kwargs.get("fps")

        if isinstance(ndimage, list):
            ndimage = np.stack(ndimage, axis=0)
            is_batch = True
        else:
            ndimage = np.asarray(ndimage)

        # check if ndimage is a batch of frames/pages (e.g. for writing GIF)
        # if mode is given, use it; otherwise fall back to image.ndim only
        if is_batch is not None:
            pass
        elif mode is not None:
            is_batch = (
                ndimage.ndim > 3 if Image.getmodebands(mode) > 1 else ndimage.ndim > 2
            )
        elif ndimage.ndim == 2:
            is_batch = False
        elif ndimage.ndim == 3 and ndimage.shape[-1] == 1:
            raise ValueError("Can't write images with one color channel.")
        elif ndimage.ndim == 3 and ndimage.shape[-1] in [2, 3, 4]:
            # Note: this makes a channel-last assumption
            is_batch = False
        else:
            is_batch = True

        if not is_batch:
            ndimage = ndimage[None, ...]

        for frame in ndimage:
            pil_frame = Image.fromarray(frame, mode=mode)
            if "bits" in kwargs:
                pil_frame = pil_frame.quantize(colors=2 ** kwargs["bits"])
            self.images_to_write.append(pil_frame)

        if (
            format is not None
            and "format" in self.save_args
            and self.save_args["format"] != format
        ):
            old_format = self.save_args["format"]
            warnings.warn(
                "Changing the output format during incremental"
                " writes is strongly discouraged."
                f" Was `{old_format}`, is now `{format}`.",
                UserWarning,
            )

        extension = self.request.extension or self.request.format_hint
        self.save_args["format"] = format or Image.registered_extensions()[extension]
        self.save_args.update(kwargs)

        # when writing to `bytes` we flush instantly
        result = None
        if self._request._uri_type == URI_BYTES:
            self._flush_writer()
            file = cast(BytesIO, self._request.get_file())
            result = file.getvalue()

        return result

    def _flush_writer(self):
        if len(self.images_to_write) == 0:
            return

        primary_image = self.images_to_write.pop(0)

        if len(self.images_to_write) > 0:
            self.save_args["save_all"] = True
            self.save_args["append_images"] = self.images_to_write

        primary_image.save(self._request.get_file(), **self.save_args)
        self.images_to_write.clear()
        self.save_args.clear()

    def get_meta(self, *, index=0) -> Dict[str, Any]:
        return self.metadata(index=index, exclude_applied=False)

    def metadata(
        self, index: int = None, exclude_applied: bool = True
    ) -> Dict[str, Any]:
        """Read ndimage metadata.

        Parameters
        ----------
        index : {integer, None}
            If the ImageResource contains multiple ndimages, and index is an
            integer, select the index-th ndimage from among them and return its
            metadata. If index is an ellipsis (...), read and return global
            metadata. If index is None, this plugin reads metadata from the
            first image of the file (index=0) unless the image is a GIF or APNG,
            in which case global metadata is read (index=...).
        exclude_applied : bool
            If True, exclude metadata fields that are applied to the image while
            reading. For example, if the binary data contains a rotation flag,
            the image is rotated by default and the rotation flag is excluded
            from the metadata to avoid confusion.

        Returns
        -------
        metadata : dict
            A dictionary of format-specific metadata.

        """

        if index is None:
            if self._image.format == "GIF":
                index = Ellipsis
            elif self._image.custom_mimetype == "image/apng":
                index = Ellipsis
            else:
                index = 0

        if isinstance(index, int) and self._image.tell() != index:
            self._image.seek(index)

        metadata = self._image.info.copy()
        metadata["mode"] = self._image.mode
        metadata["shape"] = self._image.size

        if self._image.mode == "P" and not exclude_applied:
            metadata["palette"] = np.asarray(tuple(self._image.palette.colors.keys()))

        if self._image.getexif():
            exif_data = {
                ExifTags.TAGS.get(key, "unknown"): value
                for key, value in dict(self._image.getexif()).items()
            }
            exif_data.pop("unknown", None)
            metadata.update(exif_data)

        if exclude_applied:
            metadata.pop("Orientation", None)

        return metadata

    def properties(self, index: int = None) -> ImageProperties:
        """Standardized ndimage metadata
        Parameters
        ----------
        index : int
            If the ImageResource contains multiple ndimages, and index is an
            integer, select the index-th ndimage from among them and return its
            properties. If index is an ellipsis (...), read and return the
            properties of all ndimages in the file stacked along a new batch
            dimension. If index is None, this plugin reads and returns the
            properties of the first image (index=0) unless the image is a GIF or
            APNG, in which case it reads and returns the properties all images
            (index=...).

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.

        Notes
        -----
        This does not decode pixel data and is fast for large images.

        """

        if index is None:
            if self._image.format == "GIF":
                index = Ellipsis
            elif self._image.custom_mimetype == "image/apng":
                index = Ellipsis
            else:
                index = 0

        if index is Ellipsis:
            self._image.seek(0)
        else:
            self._image.seek(index)

        if self._image.mode == "P":
            # mode of palette images is determined by their palette
            mode = self._image.palette.mode
        else:
            mode = self._image.mode

        width: int = self._image.width
        height: int = self._image.height
        shape: Tuple[int, ...] = (height, width)

        n_frames: Optional[int] = None
        if index is ...:
            n_frames = getattr(self._image, "n_frames", 1)
            shape = (n_frames, *shape)

        dummy = np.asarray(Image.new(mode, (1, 1)))
        pil_shape: Tuple[int, ...] = dummy.shape
        if len(pil_shape) > 2:
            shape = (*shape, *pil_shape[2:])

        return ImageProperties(
            shape=shape,
            dtype=dummy.dtype,
            n_images=n_frames,
            is_batch=index is Ellipsis,
        )
