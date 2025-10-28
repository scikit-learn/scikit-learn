import numpy as np

from .core.imopen import imopen


def imread(uri, *, index=None, plugin=None, extension=None, format_hint=None, **kwargs):
    """Read an ndimage from a URI.

    Opens the given URI and reads an ndimage from it. The exact behavior
    depends on both the file type and plugin used to open the file. To learn
    about the exact behavior, check the documentation of the relevant plugin.
    Typically, imread attempts to read all data stored in the URI.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    index : {int, Ellipsis, None}
        If the ImageResource contains multiple ndimages, and index is an
        integer, select the index-th ndimage from among them and return it. If
        index is an ellipsis (...), read all ndimages in the file and stack them
        along a new batch dimension. If index is None, let the plugin decide.
    plugin : {str, None}
        The plugin to use. If set to None (default) imread will perform a
        search for a matching plugin. If not None, this takes priority over
        the provided format hint  (if present).
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered.
    format_hint : str
        Deprecated. Use `extension` instead.
    **kwargs :
        Additional keyword arguments will be passed to the plugin's read call.

    Returns
    -------
    image : ndimage
        The ndimage located at the given URI.
    """

    plugin_kwargs = {
        "legacy_mode": False,
        "plugin": plugin,
        "format_hint": format_hint,
        "extension": extension,
    }

    call_kwargs = kwargs
    if index is not None:
        call_kwargs["index"] = index

    with imopen(uri, "r", **plugin_kwargs) as img_file:
        return np.asarray(img_file.read(**call_kwargs))


def imiter(uri, *, plugin=None, extension=None, format_hint=None, **kwargs):
    """Read a sequence of ndimages from a URI.

    Returns an iterable that yields ndimages from the given URI. The exact
    behavior depends on both, the file type and plugin used to open the file.
    To learn about the exact behavior, check the documentation of the relevant
    plugin.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    plugin : {str, None}
        The plugin to use. If set to None (default) imiter will perform a
        search for a matching plugin. If not None, this takes priority over
        the provided format hint (if present).
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered.
    format_hint : str
        Deprecated. Use `extension` instead.
    **kwargs :
        Additional keyword arguments will be passed to the plugin's ``iter``
        call.

    Yields
    ------
    image : ndimage
        The next ndimage located at the given URI.

    """

    with imopen(
        uri,
        "r",
        legacy_mode=False,
        plugin=plugin,
        format_hint=format_hint,
        extension=extension,
    ) as img_file:
        for image in img_file.iter(**kwargs):
            # Note: casting to ndarray here to ensure compatibility
            # with the v2.9 API
            yield np.asarray(image)


def imwrite(uri, image, *, plugin=None, extension=None, format_hint=None, **kwargs):
    """Write an ndimage to the given URI.

    The exact behavior depends on the file type and plugin used. To learn about
    the exact behavior, check the documentation of the relevant plugin.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to save the image to, e.g. a filename, pathlib.Path,
        http address or file object, check the docs for more info.
    image : np.ndarray
        The image to write to disk.
    plugin : {str, None}
        The plugin to use. If set to None (default) imwrite will perform a
        search for a matching plugin. If not None, this takes priority over
        the provided format hint (if present).
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered, and
        may also influence the format used when encoding.
    format_hint : str
        Deprecated. Use `extension` instead.
    **kwargs :
        Additional keyword arguments will be passed to the plugin's ``write``
        call.

    Returns
    -------
    encoded_image : None or Bytes
        Returns ``None`` in all cases, except when ``uri`` is set to ``<bytes>``.
        In this case it returns the encoded ndimage as a bytes string.

    """

    with imopen(
        uri,
        "w",
        legacy_mode=False,
        plugin=plugin,
        format_hint=format_hint,
        extension=extension,
    ) as img_file:
        encoded = img_file.write(image, **kwargs)

    return encoded


def improps(uri, *, index=None, plugin=None, extension=None, **kwargs):
    """Read standardized metadata.

    Opens the given URI and reads the properties of an ndimage from it. The
    properties represent standardized metadata. This means that they will have
    the same name regardless of the format being read or plugin/backend being
    used. Further, any field will be, where possible, populated with a sensible
    default (may be `None`) if the ImageResource does not declare a value in its
    metadata.

    Parameters
    ----------
    index : int
        If the ImageResource contains multiple ndimages, and index is an
        integer, select the index-th ndimage from among them and return its
        properties. If index is an ellipsis (...), read all ndimages in the file
        and stack them along a new batch dimension and return their properties.
        If index is None, let the plugin decide.
    plugin : {str, None}
        The plugin to be used. If None, performs a search for a matching
        plugin.
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered.
    **kwargs :
        Additional keyword arguments will be passed to the plugin's ``properties``
        call.

    Returns
    -------
    properties : ImageProperties
        A dataclass filled with standardized image metadata.

    Notes
    -----
    Where possible, this will avoid loading pixel data.

    See Also
    --------
    imageio.core.v3_plugin_api.ImageProperties

    """

    plugin_kwargs = {"legacy_mode": False, "plugin": plugin, "extension": extension}

    call_kwargs = kwargs
    if index is not None:
        call_kwargs["index"] = index

    with imopen(uri, "r", **plugin_kwargs) as img_file:
        properties = img_file.properties(**call_kwargs)

    return properties


def immeta(
    uri, *, index=None, plugin=None, extension=None, exclude_applied=True, **kwargs
):
    """Read format-specific metadata.

    Opens the given URI and reads metadata for an ndimage from it. The contents
    of the returned metadata dictionary is specific to both the image format and
    plugin used to open the ImageResource. To learn about the exact behavior,
    check the documentation of the relevant plugin. Typically, immeta returns a
    dictionary specific to the image format, where keys match metadata field
    names and values are a field's contents.

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path, http
        address or file object, see the docs for more info.
    index : {int, None}
        If the ImageResource contains multiple ndimages, and index is an
        integer, select the index-th ndimage from among them and return its
        metadata. If index is an ellipsis (...), return global metadata. If
        index is None, let the plugin decide the default.
    plugin : {str, None}
        The plugin to be used. If None (default), performs a search for a
        matching plugin.
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered.
    **kwargs :
        Additional keyword arguments will be passed to the plugin's metadata
        method.

    Returns
    -------
    image : ndimage
        The ndimage located at the given URI.

    """

    plugin_kwargs = {"legacy_mode": False, "plugin": plugin, "extension": extension}

    call_kwargs = kwargs
    call_kwargs["exclude_applied"] = exclude_applied
    if index is not None:
        call_kwargs["index"] = index

    with imopen(uri, "r", **plugin_kwargs) as img_file:
        metadata = img_file.metadata(**call_kwargs)

    return metadata


__all__ = ["imopen", "imread", "imwrite", "imiter", "improps", "immeta"]
