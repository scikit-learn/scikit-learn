from pathlib import Path
import warnings

from ..config import known_plugins
from ..config.extensions import known_extensions
from .request import (
    SPECIAL_READ_URIS,
    URI_FILENAME,
    InitializationError,
    IOMode,
    Request,
)


def imopen(
    uri,
    io_mode,
    *,
    plugin=None,
    extension=None,
    format_hint=None,
    legacy_mode=False,
    **kwargs,
):
    """Open an ImageResource.

    .. warning::
        This warning is for pypy users. If you are not using a context manager,
        remember to deconstruct the returned plugin to avoid leaking the file
        handle to an unclosed file.

    Parameters
    ----------
    uri : str or pathlib.Path or bytes or file or Request
        The :doc:`ImageResource <../../user_guide/requests>` to load the
        image from.
    io_mode : str
        The mode in which the file is opened. Possible values are::

            ``r`` - open the file for reading
            ``w`` - open the file for writing

        Depreciated since v2.9:
        A second character can be added to give the reader a hint on what
        the user expects. This will be ignored by new plugins and will
        only have an effect on legacy plugins. Possible values are::

            ``i`` for a single image,
            ``I`` for multiple images,
            ``v`` for a single volume,
            ``V`` for multiple volumes,
            ``?`` for don't care

    plugin : str, Plugin, or None
        The plugin to use. If set to None imopen will perform a
        search for a matching plugin. If not None, this takes priority over
        the provided format hint.
    extension : str
        If not None, treat the provided ImageResource as if it had the given
        extension. This affects the order in which backends are considered, and
        when writing this may also influence the format used when encoding.
    format_hint : str
        Deprecated. Use `extension` instead.
    legacy_mode : bool
        If true use the v2 behavior when searching for a suitable
        plugin. This will ignore v3 plugins and will check ``plugin``
        against known extensions if no plugin with the given name can be found.
    **kwargs : Any
        Additional keyword arguments will be passed to the plugin upon
        construction.

    Notes
    -----
    Registered plugins are controlled via the ``known_plugins`` dict in
    ``imageio.config``.

    Passing a ``Request`` as the uri is only supported if ``legacy_mode``
    is ``True``. In this case ``io_mode`` is ignored.

    Using the kwarg ``format_hint`` does not enforce the given format. It merely
    provides a `hint` to the selection process and plugin. The selection
    processes uses this hint for optimization; however, a plugin's decision how
    to read a ImageResource will - typically - still be based on the content of
    the resource.


    Examples
    --------

    >>> import imageio.v3 as iio
    >>> with iio.imopen("/path/to/image.png", "r") as file:
    >>>     im = file.read()

    >>> with iio.imopen("/path/to/output.jpg", "w") as file:
    >>>     file.write(im)

    """

    if isinstance(uri, Request) and legacy_mode:
        warnings.warn(
            "`iio.core.Request` is a low-level object and using it"
            " directly as input to `imopen` is discouraged. This will raise"
            " an exception in ImageIO v3.",
            DeprecationWarning,
            stacklevel=2,
        )

        request = uri
        uri = request.raw_uri
        io_mode = request.mode.io_mode
        request.format_hint = format_hint
    else:
        request = Request(uri, io_mode, format_hint=format_hint, extension=extension)

    source = "<bytes>" if isinstance(uri, bytes) else uri

    # fast-path based on plugin
    # (except in legacy mode)
    if plugin is not None:
        if isinstance(plugin, str):
            try:
                config = known_plugins[plugin]
            except KeyError:
                request.finish()
                raise ValueError(
                    f"`{plugin}` is not a registered plugin name."
                ) from None

            def loader(request, **kwargs):
                return config.plugin_class(request, **kwargs)

        else:

            def loader(request, **kwargs):
                return plugin(request, **kwargs)

        try:
            return loader(request, **kwargs)
        except InitializationError as class_specific:
            err_from = class_specific
            err_type = RuntimeError if legacy_mode else IOError
            err_msg = f"`{plugin}` can not handle the given uri."
        except ImportError:
            err_from = None
            err_type = ImportError
            err_msg = (
                f"The `{config.name}` plugin is not installed. "
                f"Use `pip install imageio[{config.install_name}]` to install it."
            )
        except Exception as generic_error:
            err_from = generic_error
            err_type = IOError
            err_msg = f"An unknown error occurred while initializing plugin `{plugin}`."

        request.finish()
        raise err_type(err_msg) from err_from

    # fast-path based on format_hint
    if request.format_hint is not None:
        for candidate_format in known_extensions[format_hint]:
            for plugin_name in candidate_format.priority:
                config = known_plugins[plugin_name]

                try:
                    candidate_plugin = config.plugin_class
                except ImportError:
                    # not installed
                    continue

                try:
                    plugin_instance = candidate_plugin(request, **kwargs)
                except InitializationError:
                    # file extension doesn't match file type
                    continue

                return plugin_instance
        else:
            resource = (
                "<bytes>" if isinstance(request.raw_uri, bytes) else request.raw_uri
            )
            warnings.warn(f"`{resource}` can not be opened as a `{format_hint}` file.")

    # fast-path based on file extension
    if request.extension in known_extensions:
        for candidate_format in known_extensions[request.extension]:
            for plugin_name in candidate_format.priority:
                config = known_plugins[plugin_name]

                try:
                    candidate_plugin = config.plugin_class
                except ImportError:
                    # not installed
                    continue

                try:
                    plugin_instance = candidate_plugin(request, **kwargs)
                except InitializationError:
                    # file extension doesn't match file type
                    continue

                return plugin_instance

    # error out for read-only special targets
    # this is hacky; can we come up with a better solution for this?
    if request.mode.io_mode == IOMode.write:
        if isinstance(uri, str) and uri.startswith(SPECIAL_READ_URIS):
            request.finish()
            err_type = ValueError if legacy_mode else IOError
            err_msg = f"`{source}` is read-only."
            raise err_type(err_msg)

    # error out for directories
    # this is a bit hacky and should be cleaned once we decide
    # how to gracefully handle DICOM
    if request._uri_type == URI_FILENAME and Path(request.raw_uri).is_dir():
        request.finish()
        err_type = ValueError if legacy_mode else IOError
        err_msg = (
            "ImageIO does not generally support reading folders. "
            "Limited support may be available via specific plugins. "
            "Specify the plugin explicitly using the `plugin` kwarg, e.g. `plugin='DICOM'`"
        )
        raise err_type(err_msg)

    # close the current request here and use fresh/new ones while trying each
    # plugin This is slow (means potentially reopening a resource several
    # times), but should only happen rarely because this is the fallback if all
    # else fails.
    request.finish()

    # fallback option: try all plugins
    for config in known_plugins.values():
        # each plugin gets its own request
        request = Request(uri, io_mode, format_hint=format_hint)

        try:
            plugin_instance = config.plugin_class(request, **kwargs)
        except InitializationError:
            continue
        except ImportError:
            continue
        else:
            return plugin_instance

    err_type = ValueError if legacy_mode else IOError
    err_msg = f"Could not find a backend to open `{source}`` with iomode `{io_mode}`."

    # check if a missing plugin could help
    if request.extension in known_extensions:
        missing_plugins = list()

        formats = known_extensions[request.extension]
        plugin_names = [
            plugin for file_format in formats for plugin in file_format.priority
        ]
        for name in plugin_names:
            config = known_plugins[name]

            try:
                config.plugin_class
                continue
            except ImportError:
                missing_plugins.append(config)

        if len(missing_plugins) > 0:
            install_candidates = "\n".join(
                [
                    (
                        f"  {config.name}:  "
                        f"pip install imageio[{config.install_name}]"
                    )
                    for config in missing_plugins
                ]
            )
            err_msg += (
                "\nBased on the extension, the following plugins might add capable backends:\n"
                f"{install_candidates}"
            )

    request.finish()
    raise err_type(err_msg)
