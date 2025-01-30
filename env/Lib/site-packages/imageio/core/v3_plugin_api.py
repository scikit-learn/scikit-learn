from . import Request
from ..typing import ArrayLike
import numpy as np
from typing import Optional, Dict, Any, Tuple, Union, List, Iterator
from dataclasses import dataclass


@dataclass
class ImageProperties:
    """Standardized Metadata

    ImageProperties represent a set of standardized metadata that is available
    under the same name for every supported format. If the ImageResource (or
    format) does not specify the value, a sensible default value is chosen
    instead.

    Attributes
    ----------
    shape : Tuple[int, ...]
        The shape of the loaded ndimage.
    dtype : np.dtype
        The dtype of the loaded ndimage.
    n_images : int
        Number of images in the file if ``index=...``, `None` for single images.
    is_batch : bool
        If True, the first dimension of the ndimage represents a batch dimension
        along which several images are stacked.
    spacing : Tuple
        A tuple describing the spacing between pixels along each axis of the
        ndimage. If the spacing is uniform along an axis the value corresponding
        to that axis is a single float. If the spacing is non-uniform, the value
        corresponding to that axis is a tuple in which the i-th element
        indicates the spacing between the i-th and (i+1)-th pixel along that
        axis.

    """

    shape: Tuple[int, ...]
    dtype: np.dtype
    n_images: Optional[int] = None
    is_batch: bool = False
    spacing: Optional[tuple] = None


class PluginV3:
    """A ImageIO Plugin.

    This is an abstract plugin that documents the v3 plugin API interface. A
    plugin is an adapter/wrapper around a backend that converts a request from
    iio.core (e.g., read an image from file) into a sequence of instructions for
    the backend that fulfill the request.

    Plugin authors may choose to subclass this class when implementing a new
    plugin, but aren't obliged to do so. As long as the plugin class implements
    the interface (methods) described below the ImageIO core will treat it just
    like any other plugin.


    Parameters
    ----------
    request : iio.Request
        A request object that represents the users intent. It provides a
        standard interface to access the various ImageResources and serves them
        to the plugin as a file object (or file). Check the docs for details.
    **kwargs : Any
        Additional configuration arguments for the plugin or backend. Usually
        these match the configuration arguments available on the backend and
        are forwarded to it.


    Raises
    ------
    InitializationError
        During ``__init__`` the plugin tests if it can fulfill the request. If
        it can't, e.g., because the request points to a file in the wrong
        format, then it should raise an ``InitializationError`` and provide a
        reason for failure. This reason may be reported to the user.
    ImportError
        Plugins will be imported dynamically when listed in
        ``iio.config.known_plugins`` to fulfill requests. This way, users only
        have to load plugins/backends they actually use. If this plugin's backend
        is not installed, it should raise an ``ImportError`` either during
        module import or during class construction.

    Notes
    -----
    Upon successful construction the plugin takes ownership of the provided
    request. This means that it is the plugin's responsibility to call
    request.finish() to close the resource when it is no longer needed.

    Plugins _must_ implement a context manager that closes and cleans any
    resources held by the plugin upon exit.

    """

    def __init__(self, request: Request) -> None:
        """Initialize a new Plugin Instance.

        See Plugin's docstring for detailed documentation.

        Notes
        -----
        The implementation here stores the request as a local variable that is
        exposed using a @property below. If you inherit from PluginV3, remember
        to call ``super().__init__(request)``.

        """

        self._request = request

    def read(self, *, index: int = 0) -> np.ndarray:
        """Read a ndimage.

        The ``read`` method loads a (single) ndimage, located at ``index`` from
        the requested ImageResource.

        It is at the plugin's descretion to decide (and document) what
        constitutes a single ndimage. A sensible way to make this decision is to
        choose based on the ImageResource's format and on what users will expect
        from such a format. For example, a sensible choice for a TIFF file
        produced by an ImageJ hyperstack is to read it as a volumetric ndimage
        (1 color dimension followed by 3 spatial dimensions). On the other hand,
        a sensible choice for a MP4 file produced by Davinci Resolve is to treat
        each frame as a ndimage (2 spatial dimensions followed by 1 color
        dimension).

        The value ``index=None`` is special. It requests the plugin to load all
        ndimages in the file and stack them along a new first axis. For example,
        if a MP4 file is read with ``index=None`` and the plugin identifies
        single frames as ndimages, then the plugin should read all frames and
        stack them into a new ndimage which now contains a time axis as its
        first axis. If a PNG file (single image format) is read with
        ``index=None`` the plugin does a very similar thing: It loads all
        ndimages in the file (here it's just one) and stacks them along a new
        first axis, effectively prepending an axis with size 1 to the image. If
        a plugin does not wish to support ``index=None`` it should set a more
        sensible default and raise a ``ValueError`` when requested to read using
        ``index=None``.

        Parameters
        ----------
        index : int
            If the ImageResource contains multiple ndimages, and index is an
            integer, select the index-th ndimage from among them and return it.
            If index is an ellipsis (...), read all ndimages in the file and
            stack them along a new batch dimension. If index is None, let the
            plugin decide. If the index is out of bounds a ``ValueError`` is
            raised.
        **kwargs : Any
            The read method may accept any number of plugin-specific keyword
            arguments to further customize the read behavior. Usually these
            match the arguments available on the backend and are forwarded to
            it.

        Returns
        -------
        ndimage : np.ndarray
            A ndimage containing decoded pixel data (sometimes called bitmap).

        Notes
        -----
        The ImageResource from which the plugin should read is managed by the
        provided request object. Directly accessing the managed ImageResource is
        _not_ permitted. Instead, you can get FileLike access to the
        ImageResource via request.get_file().

        If the backend doesn't support reading from FileLike objects, you can
        request a temporary file to pass to the backend via
        ``request.get_local_filename()``. This is, however, not very performant
        (involves copying the Request's content into a temporary file), so you
        should avoid doing this whenever possible. Consider it a fallback method
        in case all else fails.

        """
        raise NotImplementedError()

    def write(self, ndimage: Union[ArrayLike, List[ArrayLike]]) -> Optional[bytes]:
        """Write a ndimage to a ImageResource.

        The ``write`` method encodes the given ndimage into the format handled
        by the backend and writes it to the ImageResource. It overwrites
        any content that may have been previously stored in the file.

        If the backend supports only a single format then it must check if
        the ImageResource matches that format and raise an exception if not.
        Typically, this should be done during initialization in the form of a
        ``InitializationError``.

        If the backend supports more than one format it must determine the
        requested/desired format. Usually this can be done by inspecting the
        ImageResource (e.g., by checking ``request.extension``), or by providing
        a mechanism to explicitly set the format (perhaps with a - sensible -
        default value). If the plugin can not determine the desired format, it
        **must not** write to the ImageResource, but raise an exception instead.

        If the backend supports at least one format that can hold multiple
        ndimages it should be capable of handling ndimage batches and lists of
        ndimages. If the ``ndimage`` input is a list of ndimages, the plugin
        should not assume that the ndimages are not stackable, i.e., ndimages
        may have different shapes. Otherwise, the ``ndimage`` may be a batch of
        multiple ndimages stacked along the first axis of the array. The plugin
        must be able to discover this, either automatically or via additional
        `kwargs`. If there is ambiguity in the process, the plugin must clearly
        document what happens in such cases and, if possible, describe how to
        resolve this ambiguity.

        Parameters
        ----------
        ndimage : ArrayLike
            The ndimage to encode and write to the current ImageResource.
        **kwargs : Any
            The write method may accept any number of plugin-specific keyword
            arguments to customize the writing behavior. Usually these match the
            arguments available on the backend and are forwarded to it.

        Returns
        -------
        encoded_image : bytes or None
            If the chosen ImageResource is the special target ``"<bytes>"`` then
            write should return a byte string containing the encoded image data.
            Otherwise, it returns None.

        Notes
        -----
        The ImageResource to which the plugin should write to is managed by the
        provided request object. Directly accessing the managed ImageResource is
        _not_ permitted. Instead, you can get FileLike access to the
        ImageResource via request.get_file().

        If the backend doesn't support writing to FileLike objects, you can
        request a temporary file to pass to the backend via
        ``request.get_local_filename()``. This is, however, not very performant
        (involves copying the Request's content from a temporary file), so you
        should avoid doing this whenever possible. Consider it a fallback method
        in case all else fails.

        """
        raise NotImplementedError()

    def iter(self) -> Iterator[np.ndarray]:
        """Iterate the ImageResource.

        This method returns a generator that yields ndimages in the order in which
        they appear in the file. This is roughly equivalent to::

            idx = 0
            while True:
                try:
                    yield self.read(index=idx)
                except ValueError:
                    break

        It works very similar to ``read``, and you can consult the documentation
        of that method for additional information on desired behavior.

        Parameters
        ----------
        **kwargs : Any
            The iter method may accept any number of plugin-specific keyword
            arguments to further customize the reading/iteration behavior.
            Usually these match the arguments available on the backend and are
            forwarded to it.

        Yields
        ------
        ndimage : np.ndarray
            A ndimage containing decoded pixel data (sometimes called bitmap).

        See Also
        --------
        PluginV3.read

        """
        raise NotImplementedError()

    def properties(self, index: int = 0) -> ImageProperties:
        """Standardized ndimage metadata.

        Parameters
        ----------
        index : int
            If the ImageResource contains multiple ndimages, and index is an
        integer, select the index-th ndimage from among them and return its
        properties. If index is an ellipsis (...), read all ndimages in the file
        and stack them along a new batch dimension and return their properties.
        If index is None, the plugin decides the default.

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.

        """
        raise NotImplementedError()

    def metadata(self, index: int = 0, exclude_applied: bool = True) -> Dict[str, Any]:
        """Format-Specific ndimage metadata.

        The method reads metadata stored in the ImageResource and returns it as
        a python dict. The plugin is free to choose which name to give a piece
        of metadata; however, if possible, it should match the name given by the
        format. There is no requirement regarding the fields a plugin must
        expose; however, if a plugin does expose any,``exclude_applied`` applies
        to these fields.

        If the plugin does return metadata items, it must check the value of
        ``exclude_applied`` before returning them. If ``exclude applied`` is
        True, then any metadata item that would be applied to an ndimage
        returned by ``read`` (or ``iter``) must not be returned. This is done to
        avoid confusion; for example, if an ImageResource defines the ExIF
        rotation tag, and the plugin applies the rotation to the data before
        returning it, then ``exclude_applied`` prevents confusion on whether the
        tag was already applied or not.

        The `kwarg` ``index`` behaves similar to its counterpart in ``read``
        with one exception: If the ``index`` is None, then global metadata is
        returned instead of returning a combination of all metadata items. If
        there is no global metadata, the Plugin should return an empty dict or
        raise an exception.

        Parameters
        ----------
        index : int
            If the ImageResource contains multiple ndimages, and index is an
            integer, select the index-th ndimage from among them and return its
            metadata. If index is an ellipsis (...), return global metadata. If
            index is None, the plugin decides the default.
        exclude_applied : bool
            If True (default), do not report metadata fields that the plugin
            would apply/consume while reading the image.

        Returns
        -------
        metadata : dict
            A dictionary filled with format-specific metadata fields and their
            values.

        """
        raise NotImplementedError()

    def close(self) -> None:
        """Close the ImageResource.

        This method allows a plugin to behave similar to the python built-in ``open``::

            image_file = my_plugin(Request, "r")
            ...
            image_file.close()

        It is used by the context manager and deconstructor below to avoid leaking
        ImageResources. If the plugin has no other cleanup to do it doesn't have
        to overwrite this method itself and can rely on the implementation
        below.

        """

        self.request.finish()

    @property
    def request(self) -> Request:
        return self._request

    def __enter__(self) -> "PluginV3":
        return self

    def __exit__(self, type, value, traceback) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()
