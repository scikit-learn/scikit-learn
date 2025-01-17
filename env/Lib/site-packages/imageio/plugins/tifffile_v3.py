"""Read/Write TIFF files using tifffile.

.. note::
    To use this plugin you need to have `tifffile
    <https://github.com/cgohlke/tifffile>`_ installed::

        pip install tifffile

This plugin wraps tifffile, a powerful library to manipulate TIFF files. It
superseeds our previous tifffile plugin and aims to expose all the features of
tifffile.

The plugin treats individual TIFF series as ndimages. A series is a sequence of
TIFF pages that, when combined describe a meaningful unit, e.g., a volumetric
image (where each slice is stored on an individual page) or a multi-color
staining picture (where each stain is stored on an individual page). Different
TIFF flavors/variants use series in different ways and, as such, the resulting
reading behavior may vary depending on the program used while creating a
particular TIFF file.

Methods
-------
.. note::
    Check the respective function for a list of supported kwargs and detailed
    documentation.

.. autosummary::
    :toctree:

    TifffilePlugin.read
    TifffilePlugin.iter
    TifffilePlugin.write
    TifffilePlugin.properties
    TifffilePlugin.metadata

Additional methods available inside the :func:`imopen <imageio.v3.imopen>`
context:

.. autosummary::
    :toctree:

    TifffilePlugin.iter_pages

"""

from io import BytesIO
from typing import Any, Dict, Optional, cast
import warnings

import numpy as np
import tifffile

from ..core.request import URI_BYTES, InitializationError, Request
from ..core.v3_plugin_api import ImageProperties, PluginV3
from ..typing import ArrayLike


def _get_resolution(page: tifffile.TiffPage) -> Dict[str, Any]:
    metadata = {}

    try:
        metadata["resolution_unit"] = page.tags[296].value.value
    except KeyError:
        # tag 296 missing
        return metadata

    try:
        resolution_x = page.tags[282].value
        resolution_y = page.tags[283].value

        metadata["resolution"] = (
            resolution_x[0] / resolution_x[1],
            resolution_y[0] / resolution_y[1],
        )
    except KeyError:
        # tag 282 or 283 missing
        pass
    except ZeroDivisionError:
        warnings.warn(
            "Ignoring resolution metadata because at least one direction has a 0 "
            "denominator.",
            RuntimeWarning,
        )

    return metadata


class TifffilePlugin(PluginV3):
    """Support for tifffile as backend.

    Parameters
    ----------
    request : iio.Request
        A request object that represents the users intent. It provides a
        standard interface for a plugin to access the various ImageResources.
        Check the docs for details.
    kwargs : Any
        Additional kwargs are forwarded to tifffile's constructor, i.e.
        to ``TiffFile`` for reading or ``TiffWriter`` for writing.

    """

    def __init__(self, request: Request, **kwargs) -> None:
        super().__init__(request)
        self._fh = None

        if request.mode.io_mode == "r":
            try:
                self._fh = tifffile.TiffFile(request.get_file(), **kwargs)
            except tifffile.tifffile.TiffFileError:
                raise InitializationError("Tifffile can not read this file.")
        else:
            self._fh = tifffile.TiffWriter(request.get_file(), **kwargs)

    # ---------------------
    # Standard V3 Interface
    # ---------------------

    def read(self, *, index: int = None, page: int = None, **kwargs) -> np.ndarray:
        """Read a ndimage or page.

        The ndimage returned depends on the value of both ``index`` and
        ``page``. ``index`` selects the series to read and ``page`` allows
        selecting a single page from the selected series. If ``index=None``,
        ``page`` is understood as a flat index, i.e., the selection ignores
        individual series inside the file. If both ``index`` and ``page`` are
        ``None``, then all the series are read and returned as a batch.

        Parameters
        ----------
        index : int
            If ``int``, select the ndimage (series) located at that index inside
            the file and return ``page`` from it. If ``None`` and ``page`` is
            ``int`` read the page located at that (flat) index inside the file.
            If ``None`` and ``page=None``, read all ndimages from the file and
            return them as a batch.
        page : int
            If ``None`` return the full selected ndimage. If ``int``, read the
            page at the selected index and return it.
        kwargs : Any
            Additional kwargs are forwarded to TiffFile's ``as_array`` method.

        Returns
        -------
        ndarray : np.ndarray
            The decoded ndimage or page.
        """

        if "key" not in kwargs:
            kwargs["key"] = page
        elif page is not None:
            raise ValueError("Can't use `page` and `key` at the same time.")

        # set plugin default for ``index``
        if index is not None and "series" in kwargs:
            raise ValueError("Can't use `series` and `index` at the same time.")
        elif "series" in kwargs:
            index = kwargs.pop("series")
        elif index is not None:
            pass
        else:
            index = 0

        if index is Ellipsis and page is None:
            # read all series in the file and return them as a batch
            ndimage = np.stack([x for x in self.iter(**kwargs)])
        else:
            index = None if index is Ellipsis else index
            ndimage = self._fh.asarray(series=index, **kwargs)

        return ndimage

    def iter(self, **kwargs) -> np.ndarray:
        """Yield ndimages from the TIFF.

        Parameters
        ----------
        kwargs : Any
            Additional kwargs are forwarded to the TiffPageSeries' ``as_array``
            method.

        Yields
        ------
        ndimage : np.ndarray
            A decoded ndimage.
        """

        for sequence in self._fh.series:
            yield sequence.asarray(**kwargs)

    def write(
        self, ndimage: ArrayLike, *, is_batch: bool = False, **kwargs
    ) -> Optional[bytes]:
        """Save a ndimage as TIFF.

        Parameters
        ----------
        ndimage : ArrayLike
            The ndimage to encode and write to the ImageResource.
        is_batch : bool
            If True, the first dimension of the given ndimage is treated as a
            batch dimension and each element will create a new series.
        kwargs : Any
            Additional kwargs are forwarded to TiffWriter's ``write`` method.

        Returns
        -------
        encoded_image : bytes
            If the ImageResource is ``"<bytes>"``, return the encoded bytes.
            Otherwise write returns None.

        Notes
        -----
        Incremental writing is supported. Subsequent calls to ``write`` will
        create new series unless ``contiguous=True`` is used, in which case the
        call to write will append to the current series.

        """

        if not is_batch:
            ndimage = np.asarray(ndimage)[None, :]

        for image in ndimage:
            self._fh.write(image, **kwargs)

        if self._request._uri_type == URI_BYTES:
            self._fh.close()
            file = cast(BytesIO, self._request.get_file())
            return file.getvalue()

    def metadata(
        self, *, index: int = Ellipsis, page: int = None, exclude_applied: bool = True
    ) -> Dict[str, Any]:
        """Format-Specific TIFF metadata.

        The metadata returned depends on the value of both ``index`` and
        ``page``. ``index`` selects a series and ``page`` allows selecting a
        single page from the selected series. If ``index=Ellipsis``, ``page`` is
        understood as a flat index, i.e., the selection ignores individual
        series inside the file. If ``index=Ellipsis`` and ``page=None`` then
        global (file-level) metadata is returned.

        Parameters
        ----------
        index : int
            Select the series of which to extract metadata from. If Ellipsis, treat
            page as a flat index into the file's pages.
        page : int
            If not None, select the page of which to extract metadata from. If
            None, read series-level metadata or, if ``index=...`` global,
            file-level metadata.
        exclude_applied : bool
            For API compatibility. Currently ignored.

        Returns
        -------
        metadata : dict
            A dictionary with information regarding the tiff flavor (file-level)
            or tiff tags (page-level).
        """

        if index is not Ellipsis and page is not None:
            target = self._fh.series[index].pages[page]
        elif index is not Ellipsis and page is None:
            # This is based on my understanding that series-level metadata is
            # stored in the first TIFF page.
            target = self._fh.series[index].pages[0]
        elif index is Ellipsis and page is not None:
            target = self._fh.pages[page]
        else:
            target = None

        metadata = {}
        if target is None:
            # return file-level metadata
            metadata["byteorder"] = self._fh.byteorder

            for flag in tifffile.TIFF.FILE_FLAGS:
                flag_value = getattr(self._fh, "is_" + flag)
                metadata["is_" + flag] = flag_value

                if flag_value and hasattr(self._fh, flag + "_metadata"):
                    flavor_metadata = getattr(self._fh, flag + "_metadata")
                    if isinstance(flavor_metadata, tuple):
                        metadata.update(flavor_metadata[0])
                    else:
                        metadata.update(flavor_metadata)
        else:
            # tifffile may return a TiffFrame instead of a page
            target = target.keyframe

            metadata.update({tag.name: tag.value for tag in target.tags})
            metadata.update(
                {
                    "planar_configuration": target.planarconfig,
                    "compression": target.compression,
                    "predictor": target.predictor,
                    "orientation": None,  # TODO
                    "description1": target.description1,
                    "description": target.description,
                    "software": target.software,
                    **_get_resolution(target),
                    "datetime": target.datetime,
                }
            )

        return metadata

    def properties(self, *, index: int = None, page: int = None) -> ImageProperties:
        """Standardized metadata.

        The properties returned depend on the value of both ``index`` and
        ``page``. ``index`` selects a series and ``page`` allows selecting a
        single page from the selected series. If ``index=Ellipsis``, ``page`` is
        understood as a flat index, i.e., the selection ignores individual
        series inside the file. If ``index=Ellipsis`` and ``page=None`` then
        global (file-level) properties are returned. If ``index=Ellipsis``
        and ``page=...``, file-level properties for the flattened index are
        returned.

        Parameters
        ----------
        index : int
            If ``int``, select the ndimage (series) located at that index inside
            the file. If ``Ellipsis`` and ``page`` is ``int`` extract the
            properties of the page located at that (flat) index inside the file.
            If ``Ellipsis`` and ``page=None``, return the properties for the
            batch of all ndimages in the file.
        page : int
            If ``None`` return the properties of the full ndimage. If ``...``
            return the properties of the flattened index. If ``int``,
            return the properties of the page at the selected index only.

        Returns
        -------
        image_properties : ImageProperties
            The standardized metadata (properties) of the selected ndimage or series.

        """
        index = index or 0
        page_idx = 0 if page in (None, Ellipsis) else page

        if index is Ellipsis:
            target_page = self._fh.pages[page_idx]
        else:
            target_page = self._fh.series[index].pages[page_idx]

        if index is Ellipsis and page is None:
            n_series = len(self._fh.series)
            props = ImageProperties(
                shape=(n_series, *target_page.shape),
                dtype=target_page.dtype,
                n_images=n_series,
                is_batch=True,
                spacing=_get_resolution(target_page).get("resolution"),
            )
        elif index is Ellipsis and page is Ellipsis:
            n_pages = len(self._fh.pages)
            props = ImageProperties(
                shape=(n_pages, *target_page.shape),
                dtype=target_page.dtype,
                n_images=n_pages,
                is_batch=True,
                spacing=_get_resolution(target_page).get("resolution"),
            )
        else:
            props = ImageProperties(
                shape=target_page.shape,
                dtype=target_page.dtype,
                is_batch=False,
                spacing=_get_resolution(target_page).get("resolution"),
            )

        return props

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()

        super().close()

    # ------------------------------
    # Add-on Interface inside imopen
    # ------------------------------

    def iter_pages(self, index=..., **kwargs):
        """Yield pages from a TIFF file.

        This generator walks over the flat index of the pages inside an
        ImageResource and yields them in order.

        Parameters
        ----------
        index : int
            The index of the series to yield pages from. If Ellipsis, walk over
            the file's flat index (and ignore individual series).
        kwargs : Any
            Additional kwargs are passed to TiffPage's ``as_array`` method.

        Yields
        ------
        page : np.ndarray
            A page stored inside the TIFF file.

        """

        if index is Ellipsis:
            pages = self._fh.pages
        else:
            pages = self._fh.series[index]

        for page in pages:
            yield page.asarray(**kwargs)
