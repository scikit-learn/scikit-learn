import plotly.colors as clrs
from plotly.graph_objs import graph_objs as go
from plotly import exceptions
from plotly import optional_imports

np = optional_imports.get_module("numpy")
scipy_interp = optional_imports.get_module("scipy.interpolate")

from skimage import measure

# -------------------------- Layout ------------------------------


def _ternary_layout(
    title="Ternary contour plot", width=550, height=525, pole_labels=["a", "b", "c"]
):
    """
    Layout of ternary contour plot, to be passed to ``go.FigureWidget``
    object.

    Parameters
    ==========
    title : str or None
        Title of ternary plot
    width : int
        Figure width.
    height : int
        Figure height.
    pole_labels : str, default ['a', 'b', 'c']
        Names of the three poles of the triangle.
    """
    return dict(
        title=title,
        width=width,
        height=height,
        ternary=dict(
            sum=1,
            aaxis=dict(
                title=dict(text=pole_labels[0]), min=0.01, linewidth=2, ticks="outside"
            ),
            baxis=dict(
                title=dict(text=pole_labels[1]), min=0.01, linewidth=2, ticks="outside"
            ),
            caxis=dict(
                title=dict(text=pole_labels[2]), min=0.01, linewidth=2, ticks="outside"
            ),
        ),
        showlegend=False,
    )


# ------------- Transformations of coordinates -------------------


def _replace_zero_coords(ternary_data, delta=0.0005):
    """
    Replaces zero ternary coordinates with delta and normalize the new
    triplets (a, b, c).

    Parameters
    ----------

    ternary_data : ndarray of shape (N, 3)

    delta : float
        Small float to regularize logarithm.

    Notes
    -----
    Implements a method
    by J. A. Martin-Fernandez,  C. Barcelo-Vidal, V. Pawlowsky-Glahn,
    Dealing with zeros and missing values in compositional data sets
    using nonparametric imputation, Mathematical Geology 35 (2003),
    pp 253-278.
    """
    zero_mask = ternary_data == 0
    is_any_coord_zero = np.any(zero_mask, axis=0)

    unity_complement = 1 - delta * is_any_coord_zero
    if np.any(unity_complement) < 0:
        raise ValueError(
            "The provided value of delta led to negative"
            "ternary coords.Set a smaller delta"
        )
    ternary_data = np.where(zero_mask, delta, unity_complement * ternary_data)
    return ternary_data


def _ilr_transform(barycentric):
    """
    Perform Isometric Log-Ratio on barycentric (compositional) data.

    Parameters
    ----------
    barycentric: ndarray of shape (3, N)
        Barycentric coordinates.

    References
    ----------
    "An algebraic method to compute isometric logratio transformation and
    back transformation of compositional data", Jarauta-Bragulat, E.,
    Buenestado, P.; Hervada-Sala, C., in Proc. of the Annual Conf. of the
    Intl Assoc for Math Geology, 2003, pp 31-30.
    """
    barycentric = np.asarray(barycentric)
    x_0 = np.log(barycentric[0] / barycentric[1]) / np.sqrt(2)
    x_1 = (
        1.0 / np.sqrt(6) * np.log(barycentric[0] * barycentric[1] / barycentric[2] ** 2)
    )
    ilr_tdata = np.stack((x_0, x_1))
    return ilr_tdata


def _ilr_inverse(x):
    """
    Perform inverse Isometric Log-Ratio (ILR) transform to retrieve
    barycentric (compositional) data.

    Parameters
    ----------
    x : array of shape (2, N)
        Coordinates in ILR space.

    References
    ----------
    "An algebraic method to compute isometric logratio transformation and
    back transformation of compositional data", Jarauta-Bragulat, E.,
    Buenestado, P.; Hervada-Sala, C., in Proc. of the Annual Conf. of the
    Intl Assoc for Math Geology, 2003, pp 31-30.
    """
    x = np.array(x)
    matrix = np.array([[0.5, 1, 1.0], [-0.5, 1, 1.0], [0.0, 0.0, 1.0]])
    s = np.sqrt(2) / 2
    t = np.sqrt(3 / 2)
    Sk = np.einsum("ik, kj -> ij", np.array([[s, t], [-s, t]]), x)
    Z = -np.log(1 + np.exp(Sk).sum(axis=0))
    log_barycentric = np.einsum(
        "ik, kj -> ij", matrix, np.stack((2 * s * x[0], t * x[1], Z))
    )
    iilr_tdata = np.exp(log_barycentric)
    return iilr_tdata


def _transform_barycentric_cartesian():
    """
    Returns the transformation matrix from barycentric to Cartesian
    coordinates and conversely.
    """
    # reference triangle
    tri_verts = np.array([[0.5, np.sqrt(3) / 2], [0, 0], [1, 0]])
    M = np.array([tri_verts[:, 0], tri_verts[:, 1], np.ones(3)])
    return M, np.linalg.inv(M)


def _prepare_barycentric_coord(b_coords):
    """
    Check ternary coordinates and return the right barycentric coordinates.
    """
    if not isinstance(b_coords, (list, np.ndarray)):
        raise ValueError(
            "Data  should be either an array of shape (n,m),"
            "or a list of n m-lists, m=2 or 3"
        )
    b_coords = np.asarray(b_coords)
    if b_coords.shape[0] not in (2, 3):
        raise ValueError(
            "A point should have  2 (a, b) or 3 (a, b, c)" "barycentric coordinates"
        )
    if (
        (len(b_coords) == 3)
        and not np.allclose(b_coords.sum(axis=0), 1, rtol=0.01)
        and not np.allclose(b_coords.sum(axis=0), 100, rtol=0.01)
    ):
        msg = "The sum of coordinates should be 1 or 100 for all data points"
        raise ValueError(msg)

    if len(b_coords) == 2:
        A, B = b_coords
        C = 1 - (A + B)
    else:
        A, B, C = b_coords / b_coords.sum(axis=0)
    if np.any(np.stack((A, B, C)) < 0):
        raise ValueError("Barycentric coordinates should be positive.")
    return np.stack((A, B, C))


def _compute_grid(coordinates, values, interp_mode="ilr"):
    """
    Transform data points with Cartesian or ILR mapping, then Compute
    interpolation on a regular grid.

    Parameters
    ==========

    coordinates : array-like
        Barycentric coordinates of data points.
    values : 1-d array-like
        Data points, field to be represented as contours.
    interp_mode : 'ilr' (default) or 'cartesian'
        Defines how data are interpolated to compute contours.
    """
    if interp_mode == "cartesian":
        M, invM = _transform_barycentric_cartesian()
        coord_points = np.einsum("ik, kj -> ij", M, coordinates)
    elif interp_mode == "ilr":
        coordinates = _replace_zero_coords(coordinates)
        coord_points = _ilr_transform(coordinates)
    else:
        raise ValueError("interp_mode should be cartesian or ilr")
    xx, yy = coord_points[:2]
    x_min, x_max = xx.min(), xx.max()
    y_min, y_max = yy.min(), yy.max()
    n_interp = max(200, int(np.sqrt(len(values))))
    gr_x = np.linspace(x_min, x_max, n_interp)
    gr_y = np.linspace(y_min, y_max, n_interp)
    grid_x, grid_y = np.meshgrid(gr_x, gr_y)
    # We use cubic interpolation, except outside of the convex hull
    # of data points where we use nearest neighbor values.
    grid_z = scipy_interp.griddata(
        coord_points[:2].T, values, (grid_x, grid_y), method="cubic"
    )
    grid_z_other = scipy_interp.griddata(
        coord_points[:2].T, values, (grid_x, grid_y), method="nearest"
    )
    # mask_nan = np.isnan(grid_z)
    # grid_z[mask_nan] = grid_z_other[mask_nan]
    return grid_z, gr_x, gr_y


# ----------------------- Contour traces ----------------------


def _polygon_area(x, y):
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def _colors(ncontours, colormap=None):
    """
    Return a list of ``ncontours`` colors from the ``colormap`` colorscale.
    """
    if colormap in clrs.PLOTLY_SCALES.keys():
        cmap = clrs.PLOTLY_SCALES[colormap]
    else:
        raise exceptions.PlotlyError(
            "Colorscale must be a valid Plotly Colorscale."
            "The available colorscale names are {}".format(clrs.PLOTLY_SCALES.keys())
        )
    values = np.linspace(0, 1, ncontours)
    vals_cmap = np.array([pair[0] for pair in cmap])
    cols = np.array([pair[1] for pair in cmap])
    inds = np.searchsorted(vals_cmap, values)
    if "#" in cols[0]:  # for Viridis
        cols = [clrs.label_rgb(clrs.hex_to_rgb(col)) for col in cols]

    colors = [cols[0]]
    for ind, val in zip(inds[1:], values[1:]):
        val1, val2 = vals_cmap[ind - 1], vals_cmap[ind]
        interm = (val - val1) / (val2 - val1)
        col = clrs.find_intermediate_color(
            cols[ind - 1], cols[ind], interm, colortype="rgb"
        )
        colors.append(col)
    return colors


def _is_invalid_contour(x, y):
    """
    Utility function for _contour_trace

    Contours with an area of the order as 1 pixel are considered spurious.
    """
    too_small = np.all(np.abs(x - x[0]) < 2) and np.all(np.abs(y - y[0]) < 2)
    return too_small


def _extract_contours(im, values, colors):
    """
    Utility function for _contour_trace.

    In ``im`` only one part of the domain has valid values (corresponding
    to a subdomain where barycentric coordinates are well defined). When
    computing contours, we need to assign values outside of this domain.
    We can choose a value either smaller than all the values inside the
    valid domain, or larger. This value must be chose with caution so that
    no spurious contours are added. For example, if the boundary of the valid
    domain has large values and the outer value is set to a small one, all
    intermediate contours will be added at the boundary.

    Therefore, we compute the two sets of contours (with an outer value
    smaller of larger than all values in the valid domain), and choose
    the value resulting in a smaller total number of contours. There might
    be a faster way to do this, but it works...
    """
    mask_nan = np.isnan(im)
    im_min, im_max = (
        im[np.logical_not(mask_nan)].min(),
        im[np.logical_not(mask_nan)].max(),
    )
    zz_min = np.copy(im)
    zz_min[mask_nan] = 2 * im_min
    zz_max = np.copy(im)
    zz_max[mask_nan] = 2 * im_max
    all_contours1, all_values1, all_areas1, all_colors1 = [], [], [], []
    all_contours2, all_values2, all_areas2, all_colors2 = [], [], [], []
    for i, val in enumerate(values):
        contour_level1 = measure.find_contours(zz_min, val)
        contour_level2 = measure.find_contours(zz_max, val)
        all_contours1.extend(contour_level1)
        all_contours2.extend(contour_level2)
        all_values1.extend([val] * len(contour_level1))
        all_values2.extend([val] * len(contour_level2))
        all_areas1.extend(
            [_polygon_area(contour.T[1], contour.T[0]) for contour in contour_level1]
        )
        all_areas2.extend(
            [_polygon_area(contour.T[1], contour.T[0]) for contour in contour_level2]
        )
        all_colors1.extend([colors[i]] * len(contour_level1))
        all_colors2.extend([colors[i]] * len(contour_level2))
    if len(all_contours1) <= len(all_contours2):
        return all_contours1, all_values1, all_areas1, all_colors1
    else:
        return all_contours2, all_values2, all_areas2, all_colors2


def _add_outer_contour(
    all_contours,
    all_values,
    all_areas,
    all_colors,
    values,
    val_outer,
    v_min,
    v_max,
    colors,
    color_min,
    color_max,
):
    """
    Utility function for _contour_trace

    Adds the background color to fill gaps outside of computed contours.

    To compute the background color, the color of the contour with largest
    area (``val_outer``) is used. As background color, we choose the next
    color value in the direction of the extrema of the colormap.

    Then we add information for the outer contour for the different lists
    provided as arguments.

    A discrete colormap with all used colors is also returned (to be used
    by colorscale trace).
    """
    #  The exact value of outer contour is not used when defining the trace
    outer_contour = 20 * np.array([[0, 0, 1], [0, 1, 0.5]]).T
    all_contours = [outer_contour] + all_contours
    delta_values = np.diff(values)[0]
    values = np.concatenate(
        ([values[0] - delta_values], values, [values[-1] + delta_values])
    )
    colors = np.concatenate(([color_min], colors, [color_max]))
    index = np.nonzero(values == val_outer)[0][0]
    if index < len(values) / 2:
        index -= 1
    else:
        index += 1
    all_colors = [colors[index]] + all_colors
    all_values = [values[index]] + all_values
    all_areas = [0] + all_areas
    used_colors = [color for color in colors if color in all_colors]
    # Define discrete colorscale
    color_number = len(used_colors)
    scale = np.linspace(0, 1, color_number + 1)
    discrete_cm = []
    for i, color in enumerate(used_colors):
        discrete_cm.append([scale[i], used_colors[i]])
        discrete_cm.append([scale[i + 1], used_colors[i]])
    discrete_cm.append([scale[color_number], used_colors[color_number - 1]])

    return all_contours, all_values, all_areas, all_colors, discrete_cm


def _contour_trace(
    x,
    y,
    z,
    ncontours=None,
    colorscale="Electric",
    linecolor="rgb(150,150,150)",
    interp_mode="llr",
    coloring=None,
    v_min=0,
    v_max=1,
):
    """
    Contour trace in Cartesian coordinates.

    Parameters
    ==========

    x, y : array-like
        Cartesian coordinates
    z : array-like
        Field to be represented as contours.
    ncontours : int or None
        Number of contours to display (determined automatically if None).
    colorscale : None or str (Plotly colormap)
        colorscale of the contours.
    linecolor : rgb color
        Color used for lines. If ``colorscale`` is not None, line colors are
        determined from ``colorscale`` instead.
    interp_mode : 'ilr' (default) or 'cartesian'
        Defines how data are interpolated to compute contours. If 'irl',
        ILR (Isometric Log-Ratio) of compositional data is performed. If
        'cartesian', contours are determined in Cartesian space.
    coloring : None or 'lines'
        How to display contour. Filled contours if None, lines if ``lines``.
    vmin, vmax : float
        Bounds of interval of values used for the colorspace

    Notes
    =====
    """
    # Prepare colors
    # We do not take extrema, for example for one single contour
    # the color will be the middle point of the colormap
    colors = _colors(ncontours + 2, colorscale)
    # Values used for contours, extrema are not used
    # For example for a binary array [0, 1], the value of
    # the contour for ncontours=1 is 0.5.
    values = np.linspace(v_min, v_max, ncontours + 2)
    color_min, color_max = colors[0], colors[-1]
    colors = colors[1:-1]
    values = values[1:-1]

    # Color of line contours
    if linecolor is None:
        linecolor = "rgb(150, 150, 150)"
    else:
        colors = [linecolor] * ncontours

        # Retrieve all contours
    all_contours, all_values, all_areas, all_colors = _extract_contours(
        z, values, colors
    )

    # Now sort contours by decreasing area
    order = np.argsort(all_areas)[::-1]

    # Add outer contour
    all_contours, all_values, all_areas, all_colors, discrete_cm = _add_outer_contour(
        all_contours,
        all_values,
        all_areas,
        all_colors,
        values,
        all_values[order[0]],
        v_min,
        v_max,
        colors,
        color_min,
        color_max,
    )
    order = np.concatenate(([0], order + 1))

    # Compute traces, in the order of decreasing area
    traces = []
    M, invM = _transform_barycentric_cartesian()
    dx = (x.max() - x.min()) / x.size
    dy = (y.max() - y.min()) / y.size
    for index in order:
        y_contour, x_contour = all_contours[index].T
        val = all_values[index]
        if interp_mode == "cartesian":
            bar_coords = np.dot(
                invM,
                np.stack((dx * x_contour, dy * y_contour, np.ones(x_contour.shape))),
            )
        elif interp_mode == "ilr":
            bar_coords = _ilr_inverse(
                np.stack((dx * x_contour + x.min(), dy * y_contour + y.min()))
            )
        if index == 0:  # outer triangle
            a = np.array([1, 0, 0])
            b = np.array([0, 1, 0])
            c = np.array([0, 0, 1])
        else:
            a, b, c = bar_coords
        if _is_invalid_contour(x_contour, y_contour):
            continue

        _col = all_colors[index] if coloring == "lines" else linecolor
        trace = dict(
            type="scatterternary",
            a=a,
            b=b,
            c=c,
            mode="lines",
            line=dict(color=_col, shape="spline", width=1),
            fill="toself",
            fillcolor=all_colors[index],
            showlegend=True,
            hoverinfo="skip",
            name="%.3f" % val,
        )
        if coloring == "lines":
            trace["fill"] = None
        traces.append(trace)

    return traces, discrete_cm


# -------------------- Figure Factory for ternary contour -------------


def create_ternary_contour(
    coordinates,
    values,
    pole_labels=["a", "b", "c"],
    width=500,
    height=500,
    ncontours=None,
    showscale=False,
    coloring=None,
    colorscale="Bluered",
    linecolor=None,
    title=None,
    interp_mode="ilr",
    showmarkers=False,
):
    """
    Ternary contour plot.

    Parameters
    ----------

    coordinates : list or ndarray
        Barycentric coordinates of shape (2, N) or (3, N) where N is the
        number of data points. The sum of the 3 coordinates is expected
        to be 1 for all data points.
    values : array-like
        Data points of field to be represented as contours.
    pole_labels : str, default ['a', 'b', 'c']
        Names of the three poles of the triangle.
    width : int
        Figure width.
    height : int
        Figure height.
    ncontours : int or None
        Number of contours to display (determined automatically if None).
    showscale : bool, default False
        If True, a colorbar showing the color scale is displayed.
    coloring : None or 'lines'
        How to display contour. Filled contours if None, lines if ``lines``.
    colorscale : None or str (Plotly colormap)
        colorscale of the contours.
    linecolor : None or rgb color
        Color used for lines. ``colorscale`` has to be set to None, otherwise
        line colors are determined from ``colorscale``.
    title : str or None
        Title of ternary plot
    interp_mode : 'ilr' (default) or 'cartesian'
        Defines how data are interpolated to compute contours. If 'irl',
        ILR (Isometric Log-Ratio) of compositional data is performed. If
        'cartesian', contours are determined in Cartesian space.
    showmarkers : bool, default False
        If True, markers corresponding to input compositional points are
        superimposed on contours, using the same colorscale.

    Examples
    ========

    Example 1: ternary contour plot with filled contours

    >>> import plotly.figure_factory as ff
    >>> import numpy as np
    >>> # Define coordinates
    >>> a, b = np.mgrid[0:1:20j, 0:1:20j]
    >>> mask = a + b <= 1
    >>> a = a[mask].ravel()
    >>> b = b[mask].ravel()
    >>> c = 1 - a - b
    >>> # Values to be displayed as contours
    >>> z = a * b * c
    >>> fig = ff.create_ternary_contour(np.stack((a, b, c)), z)
    >>> fig.show()

    It is also possible to give only two barycentric coordinates for each
    point, since the sum of the three coordinates is one:

    >>> fig = ff.create_ternary_contour(np.stack((a, b)), z)


    Example 2: ternary contour plot with line contours

    >>> fig = ff.create_ternary_contour(np.stack((a, b, c)), z, coloring='lines')

    Example 3: customize number of contours

    >>> fig = ff.create_ternary_contour(np.stack((a, b, c)), z, ncontours=8)

    Example 4: superimpose contour plot and original data as markers

    >>> fig = ff.create_ternary_contour(np.stack((a, b, c)), z, coloring='lines',
    ...                                 showmarkers=True)

    Example 5: customize title and pole labels

    >>> fig = ff.create_ternary_contour(np.stack((a, b, c)), z,
    ...                                 title='Ternary plot',
    ...                                 pole_labels=['clay', 'quartz', 'fledspar'])
    """
    if scipy_interp is None:
        raise ImportError(
            """\
    The create_ternary_contour figure factory requires the scipy package"""
        )
    sk_measure = optional_imports.get_module("skimage")
    if sk_measure is None:
        raise ImportError(
            """\
    The create_ternary_contour figure factory requires the scikit-image
    package"""
        )
    if colorscale is None:
        showscale = False
    if ncontours is None:
        ncontours = 5
    coordinates = _prepare_barycentric_coord(coordinates)
    v_min, v_max = values.min(), values.max()
    grid_z, gr_x, gr_y = _compute_grid(coordinates, values, interp_mode=interp_mode)

    layout = _ternary_layout(
        pole_labels=pole_labels, width=width, height=height, title=title
    )

    contour_trace, discrete_cm = _contour_trace(
        gr_x,
        gr_y,
        grid_z,
        ncontours=ncontours,
        colorscale=colorscale,
        linecolor=linecolor,
        interp_mode=interp_mode,
        coloring=coloring,
        v_min=v_min,
        v_max=v_max,
    )

    fig = go.Figure(data=contour_trace, layout=layout)

    opacity = 1 if showmarkers else 0
    a, b, c = coordinates
    hovertemplate = (
        pole_labels[0]
        + ": %{a:.3f}<br>"
        + pole_labels[1]
        + ": %{b:.3f}<br>"
        + pole_labels[2]
        + ": %{c:.3f}<br>"
        "z: %{marker.color:.3f}<extra></extra>"
    )

    fig.add_scatterternary(
        a=a,
        b=b,
        c=c,
        mode="markers",
        marker={
            "color": values,
            "colorscale": colorscale,
            "line": {"color": "rgb(120, 120, 120)", "width": int(coloring != "lines")},
        },
        opacity=opacity,
        hovertemplate=hovertemplate,
    )
    if showscale:
        if not showmarkers:
            colorscale = discrete_cm
        colorbar = dict(
            {
                "type": "scatterternary",
                "a": [None],
                "b": [None],
                "c": [None],
                "marker": {
                    "cmin": values.min(),
                    "cmax": values.max(),
                    "colorscale": colorscale,
                    "showscale": True,
                },
                "mode": "markers",
            }
        )
        fig.add_trace(colorbar)

    return fig
