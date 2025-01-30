import plotly.graph_objs as go
from _plotly_utils.basevalidators import ColorscaleValidator
from ._core import apply_default_cascade, init_figure, configure_animation_controls
from .imshow_utils import rescale_intensity, _integer_ranges, _integer_types
import pandas as pd
import numpy as np
import itertools
from plotly.utils import image_array_to_data_uri

try:
    import xarray

    xarray_imported = True
except ImportError:
    xarray_imported = False

_float_types = []


def _vectorize_zvalue(z, mode="max"):
    alpha = 255 if mode == "max" else 0
    if z is None:
        return z
    elif np.isscalar(z):
        return [z] * 3 + [alpha]
    elif len(z) == 1:
        return list(z) * 3 + [alpha]
    elif len(z) == 3:
        return list(z) + [alpha]
    elif len(z) == 4:
        return z
    else:
        raise ValueError(
            "zmax can be a scalar, or an iterable of length 1, 3 or 4. "
            "A value of %s was passed for zmax." % str(z)
        )


def _infer_zmax_from_type(img):
    dt = img.dtype.type
    rtol = 1.05
    if dt in _integer_types:
        return _integer_ranges[dt][1]
    else:
        im_max = img[np.isfinite(img)].max()
        if im_max <= 1 * rtol:
            return 1
        elif im_max <= 255 * rtol:
            return 255
        elif im_max <= 65535 * rtol:
            return 65535
        else:
            return 2**32


def imshow(
    img,
    zmin=None,
    zmax=None,
    origin=None,
    labels={},
    x=None,
    y=None,
    animation_frame=None,
    facet_col=None,
    facet_col_wrap=None,
    facet_col_spacing=None,
    facet_row_spacing=None,
    color_continuous_scale=None,
    color_continuous_midpoint=None,
    range_color=None,
    title=None,
    template=None,
    width=None,
    height=None,
    aspect=None,
    contrast_rescaling=None,
    binary_string=None,
    binary_backend="auto",
    binary_compression_level=4,
    binary_format="png",
    text_auto=False,
) -> go.Figure:
    """
    Display an image, i.e. data on a 2D regular raster.

    Parameters
    ----------

    img: array-like image, or xarray
        The image data. Supported array shapes are

        - (M, N): an image with scalar data. The data is visualized
          using a colormap.
        - (M, N, 3): an image with RGB values.
        - (M, N, 4): an image with RGBA values, i.e. including transparency.

    zmin, zmax : scalar or iterable, optional
        zmin and zmax define the scalar range that the colormap covers. By default,
        zmin and zmax correspond to the min and max values of the datatype for integer
        datatypes (ie [0-255] for uint8 images, [0, 65535] for uint16 images, etc.). For
        a multichannel image of floats, the max of the image is computed and zmax is the
        smallest power of 256 (1, 255, 65535) greater than this max value,
        with a 5% tolerance. For a single-channel image, the max of the image is used.
        Overridden by range_color.

    origin : str, 'upper' or 'lower' (default 'upper')
        position of the [0, 0] pixel of the image array, in the upper left or lower left
        corner. The convention 'upper' is typically used for matrices and images.

    labels : dict with str keys and str values (default `{}`)
        Sets names used in the figure for axis titles (keys ``x`` and ``y``),
        colorbar title and hoverlabel (key ``color``). The values should correspond
        to the desired label to be displayed. If ``img`` is an xarray, dimension
        names are used for axis titles, and long name for the colorbar title
        (unless overridden in ``labels``). Possible keys are: x, y, and color.

    x, y: list-like, optional
        x and y are used to label the axes of single-channel heatmap visualizations and
        their lengths must match the lengths of the second and first dimensions of the
        img argument. They are auto-populated if the input is an xarray.

    animation_frame: int or str, optional (default None)
        axis number along which the image array is sliced to create an animation plot.
        If `img` is an xarray, `animation_frame` can be the name of one the dimensions.

    facet_col: int or str, optional (default None)
        axis number along which the image array is sliced to create a facetted plot.
        If `img` is an xarray, `facet_col` can be the name of one the dimensions.

    facet_col_wrap: int
        Maximum number of facet columns. Wraps the column variable at this width,
        so that the column facets span multiple rows.
        Ignored if `facet_col` is None.

    facet_col_spacing: float between 0 and 1
        Spacing between facet columns, in paper units. Default is 0.02.

    facet_row_spacing: float between 0 and 1
        Spacing between facet rows created when ``facet_col_wrap`` is used, in
        paper units. Default is 0.0.7.

    color_continuous_scale : str or list of str
        colormap used to map scalar data to colors (for a 2D image). This parameter is
        not used for RGB or RGBA images. If a string is provided, it should be the name
        of a known color scale, and if a list is provided, it should be a list of CSS-
        compatible colors.

    color_continuous_midpoint : number
        If set, computes the bounds of the continuous color scale to have the desired
        midpoint. Overridden by range_color or zmin and zmax.

    range_color : list of two numbers
        If provided, overrides auto-scaling on the continuous color scale, including
        overriding `color_continuous_midpoint`. Also overrides zmin and zmax. Used only
        for single-channel images.

    title : str
        The figure title.

    template : str or dict or plotly.graph_objects.layout.Template instance
        The figure template name or definition.

    width : number
        The figure width in pixels.

    height: number
        The figure height in pixels.

    aspect: 'equal', 'auto', or None
      - 'equal': Ensures an aspect ratio of 1 or pixels (square pixels)
      - 'auto': The axes is kept fixed and the aspect ratio of pixels is
        adjusted so that the data fit in the axes. In general, this will
        result in non-square pixels.
      - if None, 'equal' is used for numpy arrays and 'auto' for xarrays
        (which have typically heterogeneous coordinates)

    contrast_rescaling: 'minmax', 'infer', or None
        how to determine data values corresponding to the bounds of the color
        range, when zmin or zmax are not passed. If `minmax`, the min and max
        values of the image are used. If `infer`, a heuristic based on the image
        data type is used.

    binary_string: bool, default None
        if True, the image data are first rescaled and encoded as uint8 and
        then passed to plotly.js as a b64 PNG string. If False, data are passed
        unchanged as a numerical array. Setting to True may lead to performance
        gains, at the cost of a loss of precision depending on the original data
        type. If None, use_binary_string is set to True for multichannel (eg) RGB
        arrays, and to False for single-channel (2D) arrays. 2D arrays are
        represented as grayscale and with no colorbar if use_binary_string is
        True.

    binary_backend: str, 'auto' (default), 'pil' or 'pypng'
        Third-party package for the transformation of numpy arrays to
        png b64 strings. If 'auto', Pillow is used if installed,  otherwise
        pypng.

    binary_compression_level: int, between 0 and 9 (default 4)
        png compression level to be passed to the backend when transforming an
        array to a png b64 string. Increasing `binary_compression` decreases the
        size of the png string, but the compression step takes more time. For most
        images it is not worth using levels greater than 5, but it's possible to
        test `len(fig.data[0].source)` and to time the execution of `imshow` to
        tune the level of compression. 0 means no compression (not recommended).

    binary_format: str, 'png' (default) or 'jpg'
        compression format used to generate b64 string. 'png' is recommended
        since it uses lossless compression, but 'jpg' (lossy) compression can
        result if smaller binary strings for natural images.

    text_auto: bool or str (default `False`)
        If `True` or a string, single-channel `img` values will be displayed as text.
        A string like `'.2f'` will be interpreted as a `texttemplate` numeric formatting directive.

    Returns
    -------
    fig : graph_objects.Figure containing the displayed image

    See also
    --------

    plotly.graph_objects.Image : image trace
    plotly.graph_objects.Heatmap : heatmap trace

    Notes
    -----

    In order to update and customize the returned figure, use
    `go.Figure.update_traces` or `go.Figure.update_layout`.

    If an xarray is passed, dimensions names and coordinates are used for
    axes labels and ticks.
    """
    args = locals()
    apply_default_cascade(args)
    labels = labels.copy()
    nslices_facet = 1
    if facet_col is not None:
        if isinstance(facet_col, str):
            facet_col = img.dims.index(facet_col)
        nslices_facet = img.shape[facet_col]
        facet_slices = range(nslices_facet)
        ncols = int(facet_col_wrap) if facet_col_wrap is not None else nslices_facet
        nrows = (
            nslices_facet // ncols + 1
            if nslices_facet % ncols
            else nslices_facet // ncols
        )
    else:
        nrows = 1
        ncols = 1
    if animation_frame is not None:
        if isinstance(animation_frame, str):
            animation_frame = img.dims.index(animation_frame)
        nslices_animation = img.shape[animation_frame]
        animation_slices = range(nslices_animation)
    slice_dimensions = (facet_col is not None) + (
        animation_frame is not None
    )  # 0, 1, or 2
    facet_label = None
    animation_label = None
    img_is_xarray = False
    # ----- Define x and y, set labels if img is an xarray -------------------
    if xarray_imported and isinstance(img, xarray.DataArray):
        dims = list(img.dims)
        img_is_xarray = True
        pop_indexes = []
        if facet_col is not None:
            facet_slices = img.coords[img.dims[facet_col]].values
            pop_indexes.append(facet_col)
            facet_label = img.dims[facet_col]
        if animation_frame is not None:
            animation_slices = img.coords[img.dims[animation_frame]].values
            pop_indexes.append(animation_frame)
            animation_label = img.dims[animation_frame]
        # Remove indices in sorted order.
        for index in sorted(pop_indexes, reverse=True):
            _ = dims.pop(index)
        y_label, x_label = dims[0], dims[1]
        # np.datetime64 is not handled correctly by go.Heatmap
        for ax in [x_label, y_label]:
            if np.issubdtype(img.coords[ax].dtype, np.datetime64):
                img.coords[ax] = img.coords[ax].astype(str)
        if x is None:
            x = img.coords[x_label].values
        if y is None:
            y = img.coords[y_label].values
        if aspect is None:
            aspect = "auto"
        if labels.get("x", None) is None:
            labels["x"] = x_label
        if labels.get("y", None) is None:
            labels["y"] = y_label
        if labels.get("animation_frame", None) is None:
            labels["animation_frame"] = animation_label
        if labels.get("facet_col", None) is None:
            labels["facet_col"] = facet_label
        if labels.get("color", None) is None:
            labels["color"] = xarray.plot.utils.label_from_attrs(img)
            labels["color"] = labels["color"].replace("\n", "<br>")
    else:
        if hasattr(img, "columns") and hasattr(img.columns, "__len__"):
            if x is None:
                x = img.columns
            if labels.get("x", None) is None and hasattr(img.columns, "name"):
                labels["x"] = img.columns.name or ""
        if hasattr(img, "index") and hasattr(img.index, "__len__"):
            if y is None:
                y = img.index
            if labels.get("y", None) is None and hasattr(img.index, "name"):
                labels["y"] = img.index.name or ""

        if labels.get("x", None) is None:
            labels["x"] = ""
        if labels.get("y", None) is None:
            labels["y"] = ""
        if labels.get("color", None) is None:
            labels["color"] = ""
        if aspect is None:
            aspect = "equal"

    # --- Set the value of binary_string (forbidden for pandas)
    if isinstance(img, pd.DataFrame):
        if binary_string:
            raise ValueError("Binary strings cannot be used with pandas arrays")
        is_dataframe = True
    else:
        is_dataframe = False

    # --------------- Starting from here img is always a numpy array --------
    img = np.asanyarray(img)
    # Reshape array so that animation dimension comes first, then facets, then images
    if facet_col is not None:
        img = np.moveaxis(img, facet_col, 0)
        if animation_frame is not None and animation_frame < facet_col:
            animation_frame += 1
        facet_col = True
    if animation_frame is not None:
        img = np.moveaxis(img, animation_frame, 0)
        animation_frame = True
        args["animation_frame"] = (
            "animation_frame"
            if labels.get("animation_frame") is None
            else labels["animation_frame"]
        )
    iterables = ()
    if animation_frame is not None:
        iterables += (range(nslices_animation),)
    if facet_col is not None:
        iterables += (range(nslices_facet),)

    # Default behaviour of binary_string: True for RGB images, False for 2D
    if binary_string is None:
        binary_string = img.ndim >= (3 + slice_dimensions) and not is_dataframe

    # Cast bools to uint8 (also one byte)
    if img.dtype == bool:
        img = 255 * img.astype(np.uint8)

    if range_color is not None:
        zmin = range_color[0]
        zmax = range_color[1]

    # -------- Contrast rescaling: either minmax or infer ------------------
    if contrast_rescaling is None:
        contrast_rescaling = "minmax" if img.ndim == (2 + slice_dimensions) else "infer"

    # We try to set zmin and zmax only if necessary, because traces have good defaults
    if contrast_rescaling == "minmax":
        # When using binary_string and minmax we need to set zmin and zmax to rescale the image
        if (zmin is not None or binary_string) and zmax is None:
            zmax = img.max()
        if (zmax is not None or binary_string) and zmin is None:
            zmin = img.min()
    else:
        # For uint8 data and infer we let zmin and zmax to be None if passed as None
        if zmax is None and img.dtype != np.uint8:
            zmax = _infer_zmax_from_type(img)
        if zmin is None and zmax is not None:
            zmin = 0

    # For 2d data, use Heatmap trace, unless binary_string is True
    if img.ndim == 2 + slice_dimensions and not binary_string:
        y_index = slice_dimensions
        if y is not None and img.shape[y_index] != len(y):
            raise ValueError(
                "The length of the y vector must match the length of the first "
                + "dimension of the img matrix."
            )
        x_index = slice_dimensions + 1
        if x is not None and img.shape[x_index] != len(x):
            raise ValueError(
                "The length of the x vector must match the length of the second "
                + "dimension of the img matrix."
            )

        texttemplate = None
        if text_auto is True:
            texttemplate = "%{z}"
        elif text_auto is not False:
            texttemplate = "%{z:" + text_auto + "}"

        traces = [
            go.Heatmap(
                x=x,
                y=y,
                z=img[index_tup],
                coloraxis="coloraxis1",
                name=str(i),
                texttemplate=texttemplate,
            )
            for i, index_tup in enumerate(itertools.product(*iterables))
        ]
        autorange = True if origin == "lower" else "reversed"
        layout = dict(yaxis=dict(autorange=autorange))
        if aspect == "equal":
            layout["xaxis"] = dict(scaleanchor="y", constrain="domain")
            layout["yaxis"]["constrain"] = "domain"
        colorscale_validator = ColorscaleValidator("colorscale", "imshow")
        layout["coloraxis1"] = dict(
            colorscale=colorscale_validator.validate_coerce(
                args["color_continuous_scale"]
            ),
            cmid=color_continuous_midpoint,
            cmin=zmin,
            cmax=zmax,
        )
        if labels["color"]:
            layout["coloraxis1"]["colorbar"] = dict(title_text=labels["color"])

    # For 2D+RGB data, use Image trace
    elif (
        img.ndim >= 3
        and (img.shape[-1] in [3, 4] or slice_dimensions and binary_string)
    ) or (img.ndim == 2 and binary_string):
        rescale_image = True  # to check whether image has been modified
        if zmin is not None and zmax is not None:
            zmin, zmax = (
                _vectorize_zvalue(zmin, mode="min"),
                _vectorize_zvalue(zmax, mode="max"),
            )
        x0, y0, dx, dy = (None,) * 4
        error_msg_xarray = (
            "Non-numerical coordinates were passed with xarray `img`, but "
            "the Image trace cannot handle it. Please use `binary_string=False` "
            "for 2D data or pass instead the numpy array `img.values` to `px.imshow`."
        )
        if x is not None:
            x = np.asanyarray(x)
            if np.issubdtype(x.dtype, np.number):
                x0 = x[0]
                dx = x[1] - x[0]
            else:
                error_msg = (
                    error_msg_xarray
                    if img_is_xarray
                    else (
                        "Only numerical values are accepted for the `x` parameter "
                        "when an Image trace is used."
                    )
                )
                raise ValueError(error_msg)
        if y is not None:
            y = np.asanyarray(y)
            if np.issubdtype(y.dtype, np.number):
                y0 = y[0]
                dy = y[1] - y[0]
            else:
                error_msg = (
                    error_msg_xarray
                    if img_is_xarray
                    else (
                        "Only numerical values are accepted for the `y` parameter "
                        "when an Image trace is used."
                    )
                )
                raise ValueError(error_msg)
        if binary_string:
            if zmin is None and zmax is None:  # no rescaling, faster
                img_rescaled = img
                rescale_image = False
            elif img.ndim == 2 + slice_dimensions:  # single-channel image
                img_rescaled = rescale_intensity(
                    img, in_range=(zmin[0], zmax[0]), out_range=np.uint8
                )
            else:
                img_rescaled = np.stack(
                    [
                        rescale_intensity(
                            img[..., ch],
                            in_range=(zmin[ch], zmax[ch]),
                            out_range=np.uint8,
                        )
                        for ch in range(img.shape[-1])
                    ],
                    axis=-1,
                )
            img_str = [
                image_array_to_data_uri(
                    img_rescaled[index_tup],
                    backend=binary_backend,
                    compression=binary_compression_level,
                    ext=binary_format,
                )
                for index_tup in itertools.product(*iterables)
            ]

            traces = [
                go.Image(source=img_str_slice, name=str(i), x0=x0, y0=y0, dx=dx, dy=dy)
                for i, img_str_slice in enumerate(img_str)
            ]
        else:
            colormodel = "rgb" if img.shape[-1] == 3 else "rgba256"
            traces = [
                go.Image(
                    z=img[index_tup],
                    zmin=zmin,
                    zmax=zmax,
                    colormodel=colormodel,
                    x0=x0,
                    y0=y0,
                    dx=dx,
                    dy=dy,
                )
                for index_tup in itertools.product(*iterables)
            ]
        layout = {}
        if origin == "lower" or (dy is not None and dy < 0):
            layout["yaxis"] = dict(autorange=True)
        if dx is not None and dx < 0:
            layout["xaxis"] = dict(autorange="reversed")
    else:
        raise ValueError(
            "px.imshow only accepts 2D single-channel, RGB or RGBA images. "
            "An image of shape %s was provided. "
            "Alternatively, 3- or 4-D single or multichannel datasets can be "
            "visualized using the `facet_col` or/and `animation_frame` arguments."
            % str(img.shape)
        )

    # Now build figure
    col_labels = []
    if facet_col is not None:
        slice_label = (
            "facet_col" if labels.get("facet_col") is None else labels["facet_col"]
        )
        col_labels = [f"{slice_label}={i}" for i in facet_slices]
    fig = init_figure(args, "xy", [], nrows, ncols, col_labels, [])
    for attr_name in ["height", "width"]:
        if args[attr_name]:
            layout[attr_name] = args[attr_name]
    if args["title"]:
        layout["title_text"] = args["title"]
    elif args["template"].layout.margin.t is None:
        layout["margin"] = {"t": 60}

    frame_list = []
    for index, trace in enumerate(traces):
        if (facet_col and index < nrows * ncols) or index == 0:
            fig.add_trace(trace, row=nrows - index // ncols, col=index % ncols + 1)
    if animation_frame is not None:
        for i, index in zip(range(nslices_animation), animation_slices):
            frame_list.append(
                dict(
                    data=traces[nslices_facet * i : nslices_facet * (i + 1)],
                    layout=layout,
                    name=str(index),
                )
            )
    if animation_frame:
        fig.frames = frame_list
    fig.update_layout(layout)
    # Hover name, z or color
    if binary_string and rescale_image and not np.all(img == img_rescaled):
        # we rescaled the image, hence z is not displayed in hover since it does
        # not correspond to img values
        hovertemplate = "%s: %%{x}<br>%s: %%{y}<extra></extra>" % (
            labels["x"] or "x",
            labels["y"] or "y",
        )
    else:
        if trace["type"] == "heatmap":
            hover_name = "%{z}"
        elif img.ndim == 2:
            hover_name = "%{z[0]}"
        elif img.ndim == 3 and img.shape[-1] == 3:
            hover_name = "[%{z[0]}, %{z[1]}, %{z[2]}]"
        else:
            hover_name = "%{z}"
        hovertemplate = "%s: %%{x}<br>%s: %%{y}<br>%s: %s<extra></extra>" % (
            labels["x"] or "x",
            labels["y"] or "y",
            labels["color"] or "color",
            hover_name,
        )
    fig.update_traces(hovertemplate=hovertemplate)
    if labels["x"]:
        fig.update_xaxes(title_text=labels["x"], row=1)
    if labels["y"]:
        fig.update_yaxes(title_text=labels["y"], col=1)
    configure_animation_controls(args, go.Image, fig)
    fig.update_layout(template=args["template"], overwrite=True)
    return fig
