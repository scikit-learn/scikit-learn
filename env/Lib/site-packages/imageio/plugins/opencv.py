"""Read/Write images using OpenCV.

Backend Library: `OpenCV <https://opencv.org/>`_

This plugin wraps OpenCV (also known as ``cv2``), a popular image processing
library. Currently, it exposes OpenCVs image reading capability (no video or GIF
support yet); however, this may be added in future releases.

Methods
-------
.. note::
    Check the respective function for a list of supported kwargs and their
    documentation.

.. autosummary::
    :toctree:

    OpenCVPlugin.read
    OpenCVPlugin.iter
    OpenCVPlugin.write
    OpenCVPlugin.properties
    OpenCVPlugin.metadata

Pixel Formats (Colorspaces)
---------------------------

OpenCV is known to process images in BGR; however, most of the python ecosystem
(in particular matplotlib and other pydata libraries) use the RGB. As such,
images are converted to RGB, RGBA, or grayscale (where applicable) by default.

"""

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import cv2
import numpy as np

from ..core import Request
from ..core.request import URI_BYTES, InitializationError, IOMode
from ..core.v3_plugin_api import ImageProperties, PluginV3
from ..typing import ArrayLike


class OpenCVPlugin(PluginV3):
    def __init__(self, request: Request) -> None:
        super().__init__(request)

        self.file_handle = request.get_local_filename()
        if request._uri_type is URI_BYTES:
            self.filename = "<bytes>"
        else:
            self.filename = request.raw_uri

        mode = request.mode.io_mode
        if mode == IOMode.read and not cv2.haveImageReader(self.file_handle):
            raise InitializationError(f"OpenCV can't read `{self.filename}`.")
        elif mode == IOMode.write and not cv2.haveImageWriter(self.file_handle):
            raise InitializationError(f"OpenCV can't write to `{self.filename}`.")

    def read(
        self,
        *,
        index: int = None,
        colorspace: Union[int, str] = None,
        flags: int = cv2.IMREAD_COLOR,
    ) -> np.ndarray:
        """Read an image from the ImageResource.

        Parameters
        ----------
        index : int, Ellipsis
            If int, read the index-th image from the ImageResource. If ``...``,
            read all images from the ImageResource and stack them along a new,
            prepended, batch dimension. If None (default), use ``index=0`` if
            the image contains exactly one image and ``index=...`` otherwise.
        colorspace : str, int
            The colorspace to convert into after loading and before returning
            the image. If None (default) keep grayscale images as is, convert
            images with an alpha channel to ``RGBA`` and all other images to
            ``RGB``. If int, interpret ``colorspace`` as one of OpenCVs
            `conversion flags
            <https://docs.opencv.org/4.x/d8/d01/group__imgproc__color__conversions.html>`_
            and use it for conversion. If str, convert the image into the given
            colorspace. Possible string values are: ``"RGB"``, ``"BGR"``,
            ``"RGBA"``, ``"BGRA"``, ``"GRAY"``, ``"HSV"``, or ``"LAB"``.
        flags : int
            The OpenCV flag(s) to pass to the reader. Refer to the `OpenCV docs
            <https://docs.opencv.org/4.x/d4/da8/group__imgcodecs.html#ga288b8b3da0892bd651fce07b3bbd3a56>`_
            for details.

        Returns
        -------
        ndimage : np.ndarray
            The decoded image as a numpy array.

        """

        if index is None:
            n_images = cv2.imcount(self.file_handle, flags)
            index = 0 if n_images == 1 else ...

        if index is ...:
            retval, img = cv2.imreadmulti(self.file_handle, flags=flags)
            is_batch = True
        else:
            retval, img = cv2.imreadmulti(self.file_handle, index, 1, flags=flags)
            is_batch = False

        if retval is False:
            raise ValueError(f"Could not read index `{index}` from `{self.filename}`.")

        if img[0].ndim == 2:
            in_colorspace = "GRAY"
            out_colorspace = colorspace or "GRAY"
        elif img[0].shape[-1] == 4:
            in_colorspace = "BGRA"
            out_colorspace = colorspace or "RGBA"
        else:
            in_colorspace = "BGR"
            out_colorspace = colorspace or "RGB"

        if isinstance(colorspace, int):
            cvt_space = colorspace
        elif in_colorspace == out_colorspace.upper():
            cvt_space = None
        else:
            out_colorspace = out_colorspace.upper()
            cvt_space = getattr(cv2, f"COLOR_{in_colorspace}2{out_colorspace}")

        if cvt_space is not None:
            img = np.stack([cv2.cvtColor(x, cvt_space) for x in img])
        else:
            img = np.stack(img)

        return img if is_batch else img[0]

    def iter(
        self,
        colorspace: Union[int, str] = None,
        flags: int = cv2.IMREAD_COLOR,
    ) -> np.ndarray:
        """Yield images from the ImageResource.

        Parameters
        ----------
        colorspace : str, int
            The colorspace to convert into after loading and before returning
            the image. If None (default) keep grayscale images as is, convert
            images with an alpha channel to ``RGBA`` and all other images to
            ``RGB``. If int, interpret ``colorspace`` as one of OpenCVs
            `conversion flags
            <https://docs.opencv.org/4.x/d8/d01/group__imgproc__color__conversions.html>`_
            and use it for conversion. If str, convert the image into the given
            colorspace. Possible string values are: ``"RGB"``, ``"BGR"``,
            ``"RGBA"``, ``"BGRA"``, ``"GRAY"``, ``"HSV"``, or ``"LAB"``.
        flags : int
            The OpenCV flag(s) to pass to the reader. Refer to the `OpenCV docs
            <https://docs.opencv.org/4.x/d4/da8/group__imgcodecs.html#ga288b8b3da0892bd651fce07b3bbd3a56>`_
            for details.

        Yields
        ------
        ndimage : np.ndarray
            The decoded image as a numpy array.

        """
        for idx in range(cv2.imcount(self.file_handle)):
            yield self.read(index=idx, flags=flags, colorspace=colorspace)

    def write(
        self,
        ndimage: Union[ArrayLike, List[ArrayLike]],
        is_batch: bool = False,
        params: List[int] = None,
    ) -> Optional[bytes]:
        """Save an ndimage in the ImageResource.

        Parameters
        ----------
        ndimage : ArrayLike, List[ArrayLike]
            The image data that will be written to the file. It is either a
            single image, a batch of images, or a list of images.
        is_batch : bool
            If True, the provided ndimage is a batch of images. If False (default), the
            provided ndimage is a single image. If the provided ndimage is a list of images,
            this parameter has no effect.
        params : List[int]
            A list of parameters that will be passed to OpenCVs imwrite or
            imwritemulti functions. Possible values are documented in the
            `OpenCV documentation
            <https://docs.opencv.org/4.x/d4/da8/group__imgcodecs.html#gabbc7ef1aa2edfaa87772f1202d67e0ce>`_.

        Returns
        -------
        encoded_image : bytes, None
            If the ImageResource is ``"<bytes>"`` the call to write returns the
            encoded image as a bytes string. Otherwise it returns None.

        """

        if isinstance(ndimage, list):
            ndimage = np.stack(ndimage, axis=0)
        elif not is_batch:
            ndimage = ndimage[None, ...]

        if ndimage[0].ndim == 2:
            n_channels = 1
        else:
            n_channels = ndimage[0].shape[-1]

        if n_channels == 1:
            ndimage_cv2 = [x for x in ndimage]
        elif n_channels == 4:
            ndimage_cv2 = [cv2.cvtColor(x, cv2.COLOR_RGBA2BGRA) for x in ndimage]
        else:
            ndimage_cv2 = [cv2.cvtColor(x, cv2.COLOR_RGB2BGR) for x in ndimage]

        retval = cv2.imwritemulti(self.file_handle, ndimage_cv2, params)

        if retval is False:
            # not sure what scenario would trigger this, but
            # it can occur theoretically.
            raise IOError("OpenCV failed to write.")  # pragma: no cover

        if self.request._uri_type == URI_BYTES:
            return Path(self.file_handle).read_bytes()

    def properties(
        self,
        index: int = None,
        colorspace: Union[int, str] = None,
        flags: int = cv2.IMREAD_COLOR,
    ) -> ImageProperties:
        """Standardized image metadata.

        Parameters
        ----------
        index : int, Ellipsis
            If int, get the properties of the index-th image in the
            ImageResource. If ``...``, get the properties of the image stack
            that contains all images. If None (default), use ``index=0`` if the
            image contains exactly one image and ``index=...`` otherwise.
        colorspace : str, int
            The colorspace to convert into after loading and before returning
            the image. If None (default) keep grayscale images as is, convert
            images with an alpha channel to ``RGBA`` and all other images to
            ``RGB``. If int, interpret ``colorspace`` as one of OpenCVs
            `conversion flags
            <https://docs.opencv.org/4.x/d8/d01/group__imgproc__color__conversions.html>`_
            and use it for conversion. If str, convert the image into the given
            colorspace. Possible string values are: ``"RGB"``, ``"BGR"``,
            ``"RGBA"``, ``"BGRA"``, ``"GRAY"``, ``"HSV"``, or ``"LAB"``.
        flags : int
            The OpenCV flag(s) to pass to the reader. Refer to the `OpenCV docs
            <https://docs.opencv.org/4.x/d4/da8/group__imgcodecs.html#ga288b8b3da0892bd651fce07b3bbd3a56>`_
            for details.

        Returns
        -------
        props : ImageProperties
            A dataclass filled with standardized image metadata.

        Notes
        -----
        Reading properties with OpenCV involves decoding pixel data, because
        OpenCV doesn't provide a direct way to access metadata.

        """

        if index is None:
            n_images = cv2.imcount(self.file_handle, flags)
            is_batch = n_images > 1
        elif index is Ellipsis:
            n_images = cv2.imcount(self.file_handle, flags)
            is_batch = True
        else:
            is_batch = False

        # unfortunately, OpenCV doesn't allow reading shape without reading pixel data
        if is_batch:
            img = self.read(index=0, flags=flags, colorspace=colorspace)
            return ImageProperties(
                shape=(n_images, *img.shape),
                dtype=img.dtype,
                n_images=n_images,
                is_batch=True,
            )

        img = self.read(index=index, flags=flags, colorspace=colorspace)
        return ImageProperties(shape=img.shape, dtype=img.dtype, is_batch=False)

    def metadata(
        self, index: int = None, exclude_applied: bool = True
    ) -> Dict[str, Any]:
        """Format-specific metadata.

        .. warning::
            OpenCV does not support reading metadata. When called, this function
            will raise a ``NotImplementedError``.

        Parameters
        ----------
        index : int
            This parameter has no effect.
        exclude_applied : bool
            This parameter has no effect.

        """

        warnings.warn("OpenCV does not support reading metadata.", UserWarning)
        return dict()
