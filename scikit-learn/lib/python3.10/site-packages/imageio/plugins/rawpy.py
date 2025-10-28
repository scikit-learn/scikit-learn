""" Read/Write images using rawpy.

rawpy is an easy-to-use Python wrapper for the LibRaw library.
It also contains some extra functionality for finding and repairing hot/dead pixels.
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple, Union
import rawpy
import numpy as np

from ..core.request import URI_BYTES, InitializationError, IOMode, Request
from ..core.v3_plugin_api import ImageProperties, PluginV3
from ..typing import ArrayLike


class RawPyPlugin(PluginV3):
    """A class representing the rawpy plugin.

    Methods
    -------

    .. autosummary::
    :toctree: _plugins/rawpy

    RawPyPlugin.read
    """

    def __init__(self, request: Request) -> None:
        """Instantiates a new rawpy plugin object

        Parameters
        ----------
        request: Request
            A request object representing the resource to be operated on.
        """

        super().__init__(request)

        self._image_file = None

        if request.mode.io_mode == IOMode.read:
            try:
                self._image_file = rawpy.imread(request.get_file())
            except (
                rawpy.NotSupportedError,
                rawpy.LibRawFileUnsupportedError,
                rawpy.LibRawIOError,
            ):
                if request._uri_type == URI_BYTES:
                    raise InitializationError(
                        "RawPy can not read the provided bytes."
                    ) from None
                else:
                    raise InitializationError(
                        f"RawPy can not read {request.raw_uri}."
                    ) from None
        elif request.mode.io_mode == IOMode.write:
            raise InitializationError("RawPy does not support writing.") from None

    def close(self) -> None:
        if self._image_file:
            self._image_file.close()

        self._request.finish()

    def read(self, *, index: int = 0, **kwargs) -> np.ndarray:
        """Read Raw Image.

        Returns
        -------
        nd_image: ndarray
            The image data
        """

        nd_image: np.ndarray

        try:
            nd_image = self._image_file.postprocess(**kwargs)
        except Exception:
            pass

        if index is Ellipsis:
            nd_image = nd_image[None, ...]

        return nd_image

    def write(self, ndimage: Union[ArrayLike, List[ArrayLike]]) -> Optional[bytes]:
        """RawPy does not support writing."""
        raise NotImplementedError()

    def iter(self) -> Iterator[np.ndarray]:
        """Load the image.

        Returns
        -------
        nd_image: ndarray
            The image data
        """

        try:
            yield self.read()
        except Exception:
            pass

    def metadata(
        self, index: int = None, exclude_applied: bool = True
    ) -> Dict[str, Any]:
        """Read ndimage metadata.

        Parameters
        ----------
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

        metadata = {}

        image_size = self._image_file.sizes

        metadata["black_level_per_channel"] = self._image_file.black_level_per_channel
        metadata["camera_white_level_per_channel"] = (
            self._image_file.camera_white_level_per_channel
        )
        metadata["color_desc"] = self._image_file.color_desc
        metadata["color_matrix"] = self._image_file.color_matrix
        metadata["daylight_whitebalance"] = self._image_file.daylight_whitebalance
        metadata["dtype"] = self._image_file.raw_image.dtype
        metadata["flip"] = image_size.flip
        metadata["num_colors"] = self._image_file.num_colors
        metadata["tone_curve"] = self._image_file.tone_curve
        metadata["width"] = image_size.width
        metadata["height"] = image_size.height
        metadata["raw_width"] = image_size.raw_width
        metadata["raw_height"] = image_size.raw_height
        metadata["raw_shape"] = self._image_file.raw_image.shape
        metadata["iwidth"] = image_size.iwidth
        metadata["iheight"] = image_size.iheight
        metadata["pixel_aspect"] = image_size.pixel_aspect
        metadata["white_level"] = self._image_file.white_level

        if exclude_applied:
            metadata.pop("black_level_per_channel", None)
            metadata.pop("camera_white_level_per_channel", None)
            metadata.pop("color_desc", None)
            metadata.pop("color_matrix", None)
            metadata.pop("daylight_whitebalance", None)
            metadata.pop("dtype", None)
            metadata.pop("flip", None)
            metadata.pop("num_colors", None)
            metadata.pop("tone_curve", None)
            metadata.pop("raw_width", None)
            metadata.pop("raw_height", None)
            metadata.pop("raw_shape", None)
            metadata.pop("iwidth", None)
            metadata.pop("iheight", None)
            metadata.pop("white_level", None)

        return metadata

    def properties(self, index: int = None) -> ImageProperties:
        """Standardized ndimage metadata

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.

        Notes
        -----
        This does not decode pixel data and is fast for large images.

        """

        ImageSize = self._image_file.sizes

        width: int = ImageSize.width
        height: int = ImageSize.height
        shape: Tuple[int, ...] = (height, width)

        dtype = self._image_file.raw_image.dtype

        return ImageProperties(shape=shape, dtype=dtype)
