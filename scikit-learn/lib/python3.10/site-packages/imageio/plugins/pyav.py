"""Read/Write Videos (and images) using PyAV.

.. note::
    To use this plugin you need to have `PyAV <https://pyav.org/docs/stable/>`_
    installed::

        pip install av

This plugin wraps pyAV, a pythonic binding for the FFMPEG library. It is similar
to our FFMPEG plugin, has improved performance, features a robust interface, and
aims to supersede the FFMPEG plugin in the future.


Methods
-------
.. note::
    Check the respective function for a list of supported kwargs and detailed
    documentation.

.. autosummary::
    :toctree:

    PyAVPlugin.read
    PyAVPlugin.iter
    PyAVPlugin.write
    PyAVPlugin.properties
    PyAVPlugin.metadata

Additional methods available inside the :func:`imopen <imageio.v3.imopen>`
context:

.. autosummary::
    :toctree:

    PyAVPlugin.init_video_stream
    PyAVPlugin.write_frame
    PyAVPlugin.set_video_filter
    PyAVPlugin.container_metadata
    PyAVPlugin.video_stream_metadata

Advanced API
------------

In addition to the default ImageIO v3 API this plugin exposes custom functions
that are specific to reading/writing video and its metadata. These are available
inside the :func:`imopen <imageio.v3.imopen>` context and allow fine-grained
control over how the video is processed. The functions are documented above and
below you can find a usage example::

    import imageio.v3 as iio

    with iio.imopen("test.mp4", "w", plugin="pyav") as file:
        file.init_video_stream("libx264")
        file.container_metadata["comment"] = "This video was created using ImageIO."

        for _ in range(5):
            for frame in iio.imiter("imageio:newtonscradle.gif"):
                file.write_frame(frame)

    meta = iio.immeta("test.mp4", plugin="pyav")
    assert meta["comment"] == "This video was created using ImageIO."



Pixel Formats (Colorspaces)
---------------------------

By default, this plugin converts the video into 8-bit RGB (called ``rgb24`` in
ffmpeg). This is a useful behavior for many use-cases, but sometimes you may
want to use the video's native colorspace or you may wish to convert the video
into an entirely different colorspace. This is controlled using the ``format``
kwarg. You can use ``format=None`` to leave the image in its native colorspace
or specify any colorspace supported by FFMPEG as long as it is stridable, i.e.,
as long as it can be represented by a single numpy array. Some useful choices
include:

- rgb24 (default; 8-bit RGB)
- rgb48le (16-bit lower-endian RGB)
- bgr24 (8-bit BGR; openCVs default colorspace)
- gray (8-bit grayscale)
- yuv444p (8-bit channel-first YUV)

Further, FFMPEG maintains a list of available formats, albeit not as part of the
narrative docs. It can be `found here
<https://ffmpeg.org/doxygen/trunk/pixfmt_8h_source.html>`_ (warning: C source
code).

Filters
-------

On top of providing basic read/write functionality, this plugin allows you to
use the full collection of `video filters available in FFMPEG
<https://ffmpeg.org/ffmpeg-filters.html#Video-Filters>`_. This means that you
can apply excessive preprocessing to your video before retrieving it as a numpy
array or apply excessive post-processing before you encode your data.

Filters come in two forms: sequences or graphs. Filter sequences are, as the
name suggests, sequences of filters that are applied one after the other. They
are specified using the ``filter_sequence`` kwarg. Filter graphs, on the other
hand, come in the form of a directed graph and are specified using the
``filter_graph`` kwarg.

.. note::
    All filters are either sequences or graphs. If all you want is to apply a
    single filter, you can do this by specifying a filter sequence with a single
    entry.

A ``filter_sequence`` is a list of filters, each defined through a 2-element
tuple of the form ``(filter_name, filter_parameters)``. The first element of the
tuple is the name of the filter. The second element are the filter parameters,
which can be given either as a string or a dict. The string matches the same
format that you would use when specifying the filter using the ffmpeg
command-line tool and the dict has entries of the form ``parameter:value``. For
example::

    import imageio.v3 as iio

    # using a filter_parameters str
    img1 = iio.imread(
        "imageio:cockatoo.mp4",
        plugin="pyav",
        filter_sequence=[
            ("rotate", "45*PI/180")
        ]
    )

    # using a filter_parameters dict
    img2 = iio.imread(
        "imageio:cockatoo.mp4",
        plugin="pyav",
        filter_sequence=[
            ("rotate", {"angle":"45*PI/180", "fillcolor":"AliceBlue"})
        ]
    )

A ``filter_graph``, on the other hand, is specified using a ``(nodes, edges)``
tuple. It is best explained using an example::

    img = iio.imread(
        "imageio:cockatoo.mp4",
        plugin="pyav",
        filter_graph=(
            {
                "split": ("split", ""),
                "scale_overlay":("scale", "512:-1"),
                "overlay":("overlay", "x=25:y=25:enable='between(t,1,8)'"),
            },
            [
                ("video_in", "split", 0, 0),
                ("split", "overlay", 0, 0),
                ("split", "scale_overlay", 1, 0),
                ("scale_overlay", "overlay", 0, 1),
                ("overlay", "video_out", 0, 0),
            ]
        )
    )

The above transforms the video to have picture-in-picture of itself in the top
left corner. As you can see, nodes are specified using a dict which has names as
its keys and filter tuples as values; the same tuples as the ones used when
defining a filter sequence. Edges are a list of a 4-tuples of the form
``(node_out, node_in, output_idx, input_idx)`` and specify which two filters are
connected and which inputs/outputs should be used for this.

Further, there are two special nodes in a filter graph: ``video_in`` and
``video_out``, which represent the graph's input and output respectively. These
names can not be chosen for other nodes (those nodes would simply be
overwritten), and for a graph to be valid there must be a path from the input to
the output and all nodes in the graph must be connected.

While most graphs are quite simple, they can become very complex and we
recommend that you read through the `FFMPEG documentation
<https://ffmpeg.org/ffmpeg-filters.html#Filtergraph-description>`_ and their
examples to better understand how to use them.

"""

from fractions import Fraction
from math import ceil
from typing import Any, Dict, Generator, List, Optional, Tuple, Union

import av
import av.filter
import numpy as np
from av.codec.context import Flags
from numpy.lib.stride_tricks import as_strided

from ..core import Request
from ..core.request import URI_BYTES, InitializationError, IOMode
from ..core.v3_plugin_api import ImageProperties, PluginV3


def _format_to_dtype(format: av.VideoFormat) -> np.dtype:
    """Convert a pyAV video format into a numpy dtype"""

    if len(format.components) == 0:
        # fake format
        raise ValueError(
            f"Can't determine dtype from format `{format.name}`. It has no channels."
        )

    endian = ">" if format.is_big_endian else "<"
    dtype = "f" if "f32" in format.name else "u"
    bits_per_channel = [x.bits for x in format.components]
    n_bytes = str(int(ceil(bits_per_channel[0] / 8)))

    return np.dtype(endian + dtype + n_bytes)


def _get_frame_shape(frame: av.VideoFrame) -> Tuple[int, ...]:
    """Compute the frame's array shape

    Parameters
    ----------
    frame : av.VideoFrame
        A frame for which the resulting shape should be computed.

    Returns
    -------
    shape : Tuple[int, ...]
        A tuple describing the shape of the image data in the frame.

    """

    widths = [component.width for component in frame.format.components]
    heights = [component.height for component in frame.format.components]
    bits = np.array([component.bits for component in frame.format.components])
    line_sizes = [plane.line_size for plane in frame.planes]

    subsampled_width = widths[:-1] != widths[1:]
    subsampled_height = heights[:-1] != heights[1:]
    unaligned_components = np.any(bits % 8 != 0) or (line_sizes[:-1] != line_sizes[1:])
    if subsampled_width or subsampled_height or unaligned_components:
        raise IOError(
            f"{frame.format.name} can't be expressed as a strided array."
            "Use `format=` to select a format to convert into."
        )

    shape = [frame.height, frame.width]

    # ffmpeg doesn't have a notion of channel-first or channel-last formats
    # instead it stores frames in one or more planes which contain individual
    # components of a pixel depending on the pixel format. For channel-first
    # formats each component lives on a separate plane (n_planes) and for
    # channel-last formats all components are packed on a single plane
    # (n_channels)
    n_planes = max([component.plane for component in frame.format.components]) + 1
    if n_planes > 1:
        shape = [n_planes] + shape

    channels_per_plane = [0] * n_planes
    for component in frame.format.components:
        channels_per_plane[component.plane] += 1
    n_channels = max(channels_per_plane)

    if n_channels > 1:
        shape = shape + [n_channels]

    return tuple(shape)


def _get_frame_type(picture_type: int) -> str:
    """Return a human-readable name for provided picture type

    Parameters
    ----------
    picture_type : int
        The picture type extracted from Frame.pict_type

    Returns
    -------
    picture_name : str
        A human readable name of the picture type

    """

    if not isinstance(picture_type, int):
        # old pyAV versions send an enum, not an int
        return picture_type.name

    picture_types = [
        "NONE",
        "I",
        "P",
        "B",
        "S",
        "SI",
        "SP",
        "BI",
    ]

    return picture_types[picture_type]


class PyAVPlugin(PluginV3):
    """Support for pyAV as backend.

    Parameters
    ----------
    request : iio.Request
        A request object that represents the users intent. It provides a
        standard interface to access various the various ImageResources and
        serves them to the plugin as a file object (or file). Check the docs for
        details.
    container : str
        Only used during `iio_mode="w"`! If not None, overwrite the default container
        format chosen by pyav.
    kwargs : Any
        Additional kwargs are forwarded to PyAV's constructor.

    """

    def __init__(self, request: Request, *, container: str = None, **kwargs) -> None:
        """Initialize a new Plugin Instance.

        See Plugin's docstring for detailed documentation.

        Notes
        -----
        The implementation here stores the request as a local variable that is
        exposed using a @property below. If you inherit from PluginV3, remember
        to call ``super().__init__(request)``.

        """

        super().__init__(request)

        self._container = None
        self._video_stream = None
        self._video_filter = None

        if request.mode.io_mode == IOMode.read:
            self._next_idx = 0
            try:
                if request._uri_type == 5:  # 5 is the value of URI_HTTP
                    # pyav should read from HTTP by itself. This enables reading
                    # HTTP-based streams like DASH. Note that solving streams
                    # like this is temporary until the new request object gets
                    # implemented.
                    self._container = av.open(request.raw_uri, **kwargs)
                else:
                    self._container = av.open(request.get_file(), **kwargs)
                self._video_stream = self._container.streams.video[0]
                self._decoder = self._container.decode(video=0)
            except av.FFmpegError:
                if isinstance(request.raw_uri, bytes):
                    msg = "PyAV does not support these `<bytes>`"
                else:
                    msg = f"PyAV does not support `{request.raw_uri}`"
                raise InitializationError(msg) from None
        else:
            self.frames_written = 0
            file_handle = self.request.get_file()
            filename = getattr(file_handle, "name", None)
            extension = self.request.extension or self.request.format_hint
            if extension is None:
                raise InitializationError("Can't determine output container to use.")

            # hacky, but beats running our own format selection logic
            # (since av_guess_format is not exposed)
            try:
                setattr(file_handle, "name", filename or "tmp" + extension)
            except AttributeError:
                pass  # read-only, nothing we can do

            try:
                self._container = av.open(
                    file_handle, mode="w", format=container, **kwargs
                )
            except ValueError:
                raise InitializationError(
                    f"PyAV can not write to `{self.request.raw_uri}`"
                )

    # ---------------------
    # Standard V3 Interface
    # ---------------------

    def read(
        self,
        *,
        index: int = ...,
        format: str = "rgb24",
        filter_sequence: List[Tuple[str, Union[str, dict]]] = None,
        filter_graph: Tuple[dict, List] = None,
        constant_framerate: bool = None,
        thread_count: int = 0,
        thread_type: str = None,
    ) -> np.ndarray:
        """Read frames from the video.

        If ``index`` is an integer, this function reads the index-th frame from
        the file. If ``index`` is ... (Ellipsis), this function reads all frames
        from the video, stacks them along the first dimension, and returns a
        batch of frames.

        Parameters
        ----------
        index : int
            The index of the frame to read, e.g. ``index=5`` reads the 5th
            frame. If ``...``, read all the frames in the video and stack them
            along a new, prepended, batch dimension.
        format : str
            Set the returned colorspace. If not None (default: rgb24), convert
            the data into the given format before returning it. If ``None``
            return the data in the encoded format if it can be expressed as a
            strided array; otherwise raise an Exception.
        filter_sequence : List[str, str, dict]
            If not None, apply the given sequence of FFmpeg filters to each
            ndimage. Check the (module-level) plugin docs for details and
            examples.
        filter_graph : (dict, List)
            If not None, apply the given graph of FFmpeg filters to each
            ndimage. The graph is given as a tuple of two dicts. The first dict
            contains a (named) set of nodes, and the second dict contains a set
            of edges between nodes of the previous dict. Check the (module-level)
            plugin docs for details and examples.
        constant_framerate : bool
            If True assume the video's framerate is constant. This allows for
            faster seeking inside the file. If False, the video is reset before
            each read and searched from the beginning. If None (default), this
            value will be read from the container format.
        thread_count : int
            How many threads to use when decoding a frame. The default is 0,
            which will set the number using ffmpeg's default, which is based on
            the codec, number of available cores, threadding model, and other
            considerations.
        thread_type : str
            The threading model to be used. One of

            - `"SLICE"`: threads assemble parts of the current frame
            - `"FRAME"`: threads may assemble future frames
            - None (default): Uses ``"FRAME"`` if ``index=...`` and ffmpeg's
              default otherwise.


        Returns
        -------
        frame : np.ndarray
            A numpy array containing loaded frame data.

        Notes
        -----
        Accessing random frames repeatedly is costly (O(k), where k is the
        average distance between two keyframes). You should do so only sparingly
        if possible. In some cases, it can be faster to bulk-read the video (if
        it fits into memory) and to then access the returned ndarray randomly.

        The current implementation may cause problems for b-frames, i.e.,
        bidirectionaly predicted pictures. I lack test videos to write unit
        tests for this case.

        Reading from an index other than ``...``, i.e. reading a single frame,
        currently doesn't support filters that introduce delays.

        """

        if index is ...:
            props = self.properties(format=format)
            uses_filter = (
                self._video_filter is not None
                or filter_graph is not None
                or filter_sequence is not None
            )

            self._container.seek(0)
            if not uses_filter and props.shape[0] != 0:
                frames = np.empty(props.shape, dtype=props.dtype)
                for idx, frame in enumerate(
                    self.iter(
                        format=format,
                        filter_sequence=filter_sequence,
                        filter_graph=filter_graph,
                        thread_count=thread_count,
                        thread_type=thread_type or "FRAME",
                    )
                ):
                    frames[idx] = frame
            else:
                frames = np.stack(
                    [
                        x
                        for x in self.iter(
                            format=format,
                            filter_sequence=filter_sequence,
                            filter_graph=filter_graph,
                            thread_count=thread_count,
                            thread_type=thread_type or "FRAME",
                        )
                    ]
                )

            # reset stream container, because threading model can't change after
            # first access
            self._video_stream = self._container.streams.video[0]

            return frames

        if thread_type is not None and not (
            self._video_stream.thread_type == thread_type
            or self._video_stream.thread_type.name == thread_type
        ):
            self._video_stream.thread_type = thread_type

        if (
            thread_count != 0
            and thread_count != self._video_stream.codec_context.thread_count
        ):
            # in FFMPEG thread_count == 0 means use the default count, which we
            # change to mean don't change the thread count.
            self._video_stream.codec_context.thread_count = thread_count

        if constant_framerate is None:
            # "variable_fps" is now a flag (handle got removed). Full list at
            # https://pyav.org/docs/stable/api/container.html#module-av.format
            variable_fps = bool(self._container.format.flags & 0x400)
            constant_framerate = not variable_fps

        # note: cheap for contigous incremental reads
        self._seek(index, constant_framerate=constant_framerate)
        desired_frame = next(self._decoder)
        self._next_idx += 1

        self.set_video_filter(filter_sequence, filter_graph)
        if self._video_filter is not None:
            desired_frame = self._video_filter.send(desired_frame)

        return self._unpack_frame(desired_frame, format=format)

    def iter(
        self,
        *,
        format: str = "rgb24",
        filter_sequence: List[Tuple[str, Union[str, dict]]] = None,
        filter_graph: Tuple[dict, List] = None,
        thread_count: int = 0,
        thread_type: str = None,
    ) -> np.ndarray:
        """Yield frames from the video.

        Parameters
        ----------
        frame : np.ndarray
            A numpy array containing loaded frame data.
        format : str
            Convert the data into the given format before returning it. If None,
            return the data in the encoded format if it can be expressed as a
            strided array; otherwise raise an Exception.
        filter_sequence : List[str, str, dict]
            Set the returned colorspace. If not None (default: rgb24), convert
            the data into the given format before returning it. If ``None``
            return the data in the encoded format if it can be expressed as a
            strided array; otherwise raise an Exception.
        filter_graph : (dict, List)
            If not None, apply the given graph of FFmpeg filters to each
            ndimage. The graph is given as a tuple of two dicts. The first dict
            contains a (named) set of nodes, and the second dict contains a set
            of edges between nodes of the previous dict. Check the (module-level)
            plugin docs for details and examples.
        thread_count : int
            How many threads to use when decoding a frame. The default is 0,
            which will set the number using ffmpeg's default, which is based on
            the codec, number of available cores, threadding model, and other
            considerations.
        thread_type : str
            The threading model to be used. One of

            - `"SLICE"` (default): threads assemble parts of the current frame
            - `"FRAME"`: threads may assemble future frames (faster for bulk reading)


        Yields
        ------
        frame : np.ndarray
            A (decoded) video frame.


        """

        self._video_stream.thread_type = thread_type or "SLICE"
        self._video_stream.codec_context.thread_count = thread_count

        self.set_video_filter(filter_sequence, filter_graph)

        for frame in self._decoder:
            self._next_idx += 1

            if self._video_filter is not None:
                try:
                    frame = self._video_filter.send(frame)
                except StopIteration:
                    break

            if frame is None:
                continue

            yield self._unpack_frame(frame, format=format)

        if self._video_filter is not None:
            for frame in self._video_filter:
                yield self._unpack_frame(frame, format=format)

    def write(
        self,
        ndimage: Union[np.ndarray, List[np.ndarray]],
        *,
        codec: str = None,
        is_batch: bool = True,
        fps: int = 24,
        in_pixel_format: str = "rgb24",
        out_pixel_format: str = None,
        filter_sequence: List[Tuple[str, Union[str, dict]]] = None,
        filter_graph: Tuple[dict, List] = None,
    ) -> Optional[bytes]:
        """Save a ndimage as a video.

        Given a batch of frames (stacked along the first axis) or a list of
        frames, encode them and add the result to the ImageResource.

        Parameters
        ----------
        ndimage : ArrayLike, List[ArrayLike]
            The ndimage to encode and write to the ImageResource.
        codec : str
            The codec to use when encoding frames. Only needed on first write
            and ignored on subsequent writes.
        is_batch : bool
            If True (default), the ndimage is a batch of images, otherwise it is
            a single image. This parameter has no effect on lists of ndimages.
        fps : str
            The resulting videos frames per second.
        in_pixel_format : str
            The pixel format of the incoming ndarray. Defaults to "rgb24" and can
            be any stridable pix_fmt supported by FFmpeg.
        out_pixel_format : str
            The pixel format to use while encoding frames. If None (default)
            use the codec's default.
        filter_sequence : List[str, str, dict]
            If not None, apply the given sequence of FFmpeg filters to each
            ndimage. Check the (module-level) plugin docs for details and
            examples.
        filter_graph : (dict, List)
            If not None, apply the given graph of FFmpeg filters to each
            ndimage. The graph is given as a tuple of two dicts. The first dict
            contains a (named) set of nodes, and the second dict contains a set
            of edges between nodes of the previous dict. Check the (module-level)
            plugin docs for details and examples.

        Returns
        -------
        encoded_image : bytes or None
            If the chosen ImageResource is the special target ``"<bytes>"`` then
            write will return a byte string containing the encoded image data.
            Otherwise, it returns None.

        Notes
        -----
        When writing ``<bytes>``, the video is finalized immediately after the
        first write call and calling write multiple times to append frames is
        not possible.

        """

        if isinstance(ndimage, list):
            # frames shapes must agree for video
            if any(f.shape != ndimage[0].shape for f in ndimage):
                raise ValueError("All frames should have the same shape")
        elif not is_batch:
            ndimage = np.asarray(ndimage)[None, ...]
        else:
            ndimage = np.asarray(ndimage)

        if self._video_stream is None:
            self.init_video_stream(codec, fps=fps, pixel_format=out_pixel_format)

        self.set_video_filter(filter_sequence, filter_graph)

        for img in ndimage:
            self.write_frame(img, pixel_format=in_pixel_format)

        if self.request._uri_type == URI_BYTES:
            # bytes are immutuable, so we have to flush immediately
            # and can't support appending
            self._flush_writer()
            self._container.close()

            return self.request.get_file().getvalue()

    def properties(self, index: int = ..., *, format: str = "rgb24") -> ImageProperties:
        """Standardized ndimage metadata.

        Parameters
        ----------
        index : int
            The index of the ndimage for which to return properties. If ``...``
            (Ellipsis, default), return the properties for the resulting batch
            of frames.
        format : str
            If not None (default: rgb24), convert the data into the given format
            before returning it. If None return the data in the encoded format
            if that can be expressed as a strided array; otherwise raise an
            Exception.

        Returns
        -------
        properties : ImageProperties
            A dataclass filled with standardized image metadata.

        Notes
        -----
        This function is efficient and won't process any pixel data.

        The provided metadata does not include modifications by any filters
        (through ``filter_sequence`` or ``filter_graph``).

        """

        video_width = self._video_stream.codec_context.width
        video_height = self._video_stream.codec_context.height
        pix_format = format or self._video_stream.codec_context.pix_fmt
        frame_template = av.VideoFrame(video_width, video_height, pix_format)

        shape = _get_frame_shape(frame_template)
        if index is ...:
            n_frames = self._video_stream.frames
            shape = (n_frames,) + shape

        return ImageProperties(
            shape=tuple(shape),
            dtype=_format_to_dtype(frame_template.format),
            n_images=shape[0] if index is ... else None,
            is_batch=index is ...,
        )

    def metadata(
        self,
        index: int = ...,
        exclude_applied: bool = True,
        constant_framerate: bool = None,
    ) -> Dict[str, Any]:
        """Format-specific metadata.

        Returns a dictionary filled with metadata that is either stored in the
        container, the video stream, or the frame's side-data.

        Parameters
        ----------
        index : int
            If ... (Ellipsis, default) return global metadata (the metadata
            stored in the container and video stream). If not ..., return the
            side data stored in the frame at the given index.
        exclude_applied : bool
            Currently, this parameter has no effect. It exists for compliance with
            the ImageIO v3 API.
        constant_framerate : bool
            If True assume the video's framerate is constant. This allows for
            faster seeking inside the file. If False, the video is reset before
            each read and searched from the beginning. If None (default), this
            value will be read from the container format.

        Returns
        -------
        metadata : dict
            A dictionary filled with format-specific metadata fields and their
            values.

        """

        metadata = dict()

        if index is ...:
            # useful flags defined on the container and/or video stream
            metadata.update(
                {
                    "video_format": self._video_stream.codec_context.pix_fmt,
                    "codec": self._video_stream.codec.name,
                    "long_codec": self._video_stream.codec.long_name,
                    "profile": self._video_stream.profile,
                    "fps": float(self._video_stream.guessed_rate),
                }
            )
            if self._video_stream.duration is not None:
                duration = float(
                    self._video_stream.duration * self._video_stream.time_base
                )
                metadata.update({"duration": duration})

            metadata.update(self.container_metadata)
            metadata.update(self.video_stream_metadata)
            return metadata

        if constant_framerate is None:
            # "variable_fps" is now a flag (handle got removed). Full list at
            # https://pyav.org/docs/stable/api/container.html#module-av.format
            variable_fps = bool(self._container.format.flags & 0x400)
            constant_framerate = not variable_fps

        self._seek(index, constant_framerate=constant_framerate)
        desired_frame = next(self._decoder)
        self._next_idx += 1

        # useful flags defined on the frame
        metadata.update(
            {
                "key_frame": bool(desired_frame.key_frame),
                "time": desired_frame.time,
                "interlaced_frame": bool(desired_frame.interlaced_frame),
                "frame_type": _get_frame_type(desired_frame.pict_type),
            }
        )

        # side data
        metadata.update(
            {item.type.name: bytes(item) for item in desired_frame.side_data}
        )

        return metadata

    def close(self) -> None:
        """Close the Video."""

        is_write = self.request.mode.io_mode == IOMode.write
        if is_write and self._video_stream is not None:
            self._flush_writer()

        if self._video_stream is not None:
            self._video_stream = None

        if self._container is not None:
            self._container.close()

        self.request.finish()

    def __enter__(self) -> "PyAVPlugin":
        return super().__enter__()

    # ------------------------------
    # Add-on Interface inside imopen
    # ------------------------------

    def init_video_stream(
        self,
        codec: str,
        *,
        fps: float = 24,
        pixel_format: str = None,
        max_keyframe_interval: int = None,
        force_keyframes: bool = None,
    ) -> None:
        """Initialize a new video stream.

        This function adds a new video stream to the ImageResource using the
        selected encoder (codec), framerate, and colorspace.

        Parameters
        ----------
        codec : str
            The codec to use, e.g. ``"h264"`` or ``"vp9"``.
        fps : float
            The desired framerate of the video stream (frames per second).
        pixel_format : str
            The pixel format to use while encoding frames. If None (default) use
            the codec's default.
        max_keyframe_interval : int
            The maximum distance between two intra frames (I-frames). Also known
            as GOP size. If unspecified use the codec's default. Note that not
            every I-frame is a keyframe; see the notes for details.
        force_keyframes : bool
            If True, limit inter frames dependency to frames within the current
            keyframe interval (GOP), i.e., force every I-frame to be a keyframe.
            If unspecified, use the codec's default.

        Notes
        -----
        You can usually leave ``max_keyframe_interval`` and ``force_keyframes``
        at their default values, unless you try to generate seek-optimized video
        or have a similar specialist use-case. In this case, ``force_keyframes``
        controls the ability to seek to _every_ I-frame, and
        ``max_keyframe_interval`` controls how close to a random frame you can
        seek. Low values allow more fine-grained seek at the expense of
        file-size (and thus I/O performance).

        """

        fps = Fraction.from_float(fps)
        stream = self._container.add_stream(codec, fps)
        stream.time_base = Fraction(1 / fps).limit_denominator(int(2**16 - 1))
        if pixel_format is not None:
            stream.pix_fmt = pixel_format
        if max_keyframe_interval is not None:
            stream.gop_size = max_keyframe_interval
        if force_keyframes is not None:
            if force_keyframes:
                stream.codec_context.flags |= Flags.closed_gop
            else:
                stream.codec_context.flags &= ~Flags.closed_gop

        self._video_stream = stream

    def write_frame(self, frame: np.ndarray, *, pixel_format: str = "rgb24") -> None:
        """Add a frame to the video stream.

        This function appends a new frame to the video. It assumes that the
        stream previously has been initialized. I.e., ``init_video_stream`` has
        to be called before calling this function for the write to succeed.

        Parameters
        ----------
        frame : np.ndarray
            The image to be appended/written to the video stream.
        pixel_format : str
            The colorspace (pixel format) of the incoming frame.

        Notes
        -----
        Frames may be held in a buffer, e.g., by the filter pipeline used during
        writing or by FFMPEG to batch them prior to encoding. Make sure to
        ``.close()`` the plugin or to use a context manager to ensure that all
        frames are written to the ImageResource.

        """

        # manual packing of ndarray into frame
        # (this should live in pyAV, but it doesn't support all the formats we
        # want and PRs there are slow)
        pixel_format = av.VideoFormat(pixel_format)
        img_dtype = _format_to_dtype(pixel_format)
        width = frame.shape[2 if pixel_format.is_planar else 1]
        height = frame.shape[1 if pixel_format.is_planar else 0]
        av_frame = av.VideoFrame(width, height, pixel_format.name)
        if pixel_format.is_planar:
            for idx, plane in enumerate(av_frame.planes):
                plane_array = np.frombuffer(plane, dtype=img_dtype)
                plane_array = as_strided(
                    plane_array,
                    shape=(plane.height, plane.width),
                    strides=(plane.line_size, img_dtype.itemsize),
                )
                plane_array[...] = frame[idx]
        else:
            if pixel_format.name.startswith("bayer_"):
                # ffmpeg doesn't describe bayer formats correctly
                # see https://github.com/imageio/imageio/issues/761#issuecomment-1059318851
                # and following for details.
                n_channels = 1
            else:
                n_channels = len(pixel_format.components)

            plane = av_frame.planes[0]
            plane_shape = (plane.height, plane.width)
            plane_strides = (plane.line_size, n_channels * img_dtype.itemsize)
            if n_channels > 1:
                plane_shape += (n_channels,)
                plane_strides += (img_dtype.itemsize,)

            plane_array = as_strided(
                np.frombuffer(plane, dtype=img_dtype),
                shape=plane_shape,
                strides=plane_strides,
            )
            plane_array[...] = frame

        stream = self._video_stream
        av_frame.time_base = stream.codec_context.time_base
        av_frame.pts = self.frames_written
        self.frames_written += 1

        if self._video_filter is not None:
            av_frame = self._video_filter.send(av_frame)
            if av_frame is None:
                return

        if stream.frames == 0:
            stream.width = av_frame.width
            stream.height = av_frame.height

        for packet in stream.encode(av_frame):
            self._container.mux(packet)

    def set_video_filter(
        self,
        filter_sequence: List[Tuple[str, Union[str, dict]]] = None,
        filter_graph: Tuple[dict, List] = None,
    ) -> None:
        """Set the filter(s) to use.

        This function creates a new FFMPEG filter graph to use when reading or
        writing video. In the case of reading, frames are passed through the
        filter graph before begin returned and, in case of writing, frames are
        passed through the filter before being written to the video.

        Parameters
        ----------
        filter_sequence : List[str, str, dict]
            If not None, apply the given sequence of FFmpeg filters to each
            ndimage. Check the (module-level) plugin docs for details and
            examples.
        filter_graph : (dict, List)
            If not None, apply the given graph of FFmpeg filters to each
            ndimage. The graph is given as a tuple of two dicts. The first dict
            contains a (named) set of nodes, and the second dict contains a set
            of edges between nodes of the previous dict. Check the
            (module-level) plugin docs for details and examples.

        Notes
        -----
        Changing a filter graph with lag during reading or writing will
        currently cause frames in the filter queue to be lost.

        """

        if filter_sequence is None and filter_graph is None:
            self._video_filter = None
            return

        if filter_sequence is None:
            filter_sequence = list()

        node_descriptors: Dict[str, Tuple[str, Union[str, Dict]]]
        edges: List[Tuple[str, str, int, int]]
        if filter_graph is None:
            node_descriptors, edges = dict(), [("video_in", "video_out", 0, 0)]
        else:
            node_descriptors, edges = filter_graph

        graph = av.filter.Graph()

        previous_node = graph.add_buffer(template=self._video_stream)
        for filter_name, argument in filter_sequence:
            if isinstance(argument, str):
                current_node = graph.add(filter_name, argument)
            else:
                current_node = graph.add(filter_name, **argument)
            previous_node.link_to(current_node)
            previous_node = current_node

        nodes = dict()
        nodes["video_in"] = previous_node
        nodes["video_out"] = graph.add("buffersink")
        for name, (filter_name, arguments) in node_descriptors.items():
            if isinstance(arguments, str):
                nodes[name] = graph.add(filter_name, arguments)
            else:
                nodes[name] = graph.add(filter_name, **arguments)

        for from_note, to_node, out_idx, in_idx in edges:
            nodes[from_note].link_to(nodes[to_node], out_idx, in_idx)

        graph.configure()

        def video_filter():
            # this starts a co-routine
            # send frames using graph.send()
            frame = yield None

            # send and receive frames in "parallel"
            while frame is not None:
                graph.push(frame)
                try:
                    frame = yield graph.pull()
                except av.error.BlockingIOError:
                    # filter has lag and needs more frames
                    frame = yield None
                except av.error.EOFError:
                    break

            try:
                # send EOF in av>=9.0
                graph.push(None)
            except ValueError:  # pragma: no cover
                # handle av<9.0
                pass

            # all frames have been sent, empty the filter
            while True:
                try:
                    yield graph.pull()
                except av.error.EOFError:
                    break  # EOF
                except av.error.BlockingIOError:  # pragma: no cover
                    # handle av<9.0
                    break

        self._video_filter = video_filter()
        self._video_filter.send(None)

    @property
    def container_metadata(self):
        """Container-specific metadata.

        A dictionary containing metadata stored at the container level.

        """
        return self._container.metadata

    @property
    def video_stream_metadata(self):
        """Stream-specific metadata.

        A dictionary containing metadata stored at the stream level.

        """
        return self._video_stream.metadata

    # -------------------------------
    # Internals and private functions
    # -------------------------------

    def _unpack_frame(self, frame: av.VideoFrame, *, format: str = None) -> np.ndarray:
        """Convert a av.VideoFrame into a ndarray

        Parameters
        ----------
        frame : av.VideoFrame
            The frame to unpack.
        format : str
            If not None, convert the frame to the given format before unpacking.

        """

        if format is not None:
            frame = frame.reformat(format=format)

        dtype = _format_to_dtype(frame.format)
        shape = _get_frame_shape(frame)

        planes = list()
        for idx in range(len(frame.planes)):
            n_channels = sum(
                [
                    x.bits // (dtype.itemsize * 8)
                    for x in frame.format.components
                    if x.plane == idx
                ]
            )
            av_plane = frame.planes[idx]
            plane_shape = (av_plane.height, av_plane.width)
            plane_strides = (av_plane.line_size, n_channels * dtype.itemsize)
            if n_channels > 1:
                plane_shape += (n_channels,)
                plane_strides += (dtype.itemsize,)

            np_plane = as_strided(
                np.frombuffer(av_plane, dtype=dtype),
                shape=plane_shape,
                strides=plane_strides,
            )
            planes.append(np_plane)

        if len(planes) > 1:
            # Note: the planes *should* exist inside a contigous memory block
            # somewhere inside av.Frame however pyAV does not appear to expose this,
            # so we are forced to copy the planes individually instead of wrapping
            # them :(
            out = np.concatenate(planes).reshape(shape)
        else:
            out = planes[0]

        return out

    def _seek(self, index, *, constant_framerate: bool = True) -> Generator:
        """Seeks to the frame at the given index."""

        if index == self._next_idx:
            return  # fast path :)

        # we must decode at least once before we seek otherwise the
        # returned frames become corrupt.
        if self._next_idx == 0:
            next(self._decoder)
            self._next_idx += 1

            if index == self._next_idx:
                return  # fast path :)

        # remove this branch until I find a way to efficiently find the next
        # keyframe. keeping this as a reminder
        # if self._next_idx < index and index < self._next_keyframe_idx:
        #     frames_to_yield = index - self._next_idx
        if not constant_framerate and index > self._next_idx:
            frames_to_yield = index - self._next_idx
        elif not constant_framerate:
            # seek backwards and can't link idx and pts
            self._container.seek(0)
            self._decoder = self._container.decode(video=0)
            self._next_idx = 0

            frames_to_yield = index
        else:
            # we know that the time between consecutive frames is constant
            # hence we can link index and pts

            # how many pts lie between two frames
            sec_delta = 1 / self._video_stream.guessed_rate
            pts_delta = sec_delta / self._video_stream.time_base

            index_pts = int(index * pts_delta)

            # this only seeks to the closed (preceeding) keyframe
            self._container.seek(index_pts, stream=self._video_stream)
            self._decoder = self._container.decode(video=0)

            # this may be made faster if we could get the keyframe's time without
            # decoding it
            keyframe = next(self._decoder)
            keyframe_time = keyframe.pts * keyframe.time_base
            keyframe_pts = int(keyframe_time / self._video_stream.time_base)
            keyframe_index = keyframe_pts // pts_delta

            self._container.seek(index_pts, stream=self._video_stream)
            self._next_idx = keyframe_index

            frames_to_yield = index - keyframe_index

        for _ in range(frames_to_yield):
            next(self._decoder)
            self._next_idx += 1

    def _flush_writer(self):
        """Flush the filter and encoder

        This will reset the filter to `None` and send EoF to the encoder,
        i.e., after calling, no more frames may be written.

        """

        stream = self._video_stream

        if self._video_filter is not None:
            # flush encoder
            for av_frame in self._video_filter:
                if stream.frames == 0:
                    stream.width = av_frame.width
                    stream.height = av_frame.height
                for packet in stream.encode(av_frame):
                    self._container.mux(packet)
            self._video_filter = None

        # flush stream
        for packet in stream.encode():
            self._container.mux(packet)
        self._video_stream = None
