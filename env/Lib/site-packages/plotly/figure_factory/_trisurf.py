from plotly import exceptions, optional_imports
import plotly.colors as clrs
from plotly.graph_objs import graph_objs

np = optional_imports.get_module("numpy")


def map_face2color(face, colormap, scale, vmin, vmax):
    """
    Normalize facecolor values by vmin/vmax and return rgb-color strings

    This function takes a tuple color along with a colormap and a minimum
    (vmin) and maximum (vmax) range of possible mean distances for the
    given parametrized surface. It returns an rgb color based on the mean
    distance between vmin and vmax

    """
    if vmin >= vmax:
        raise exceptions.PlotlyError(
            "Incorrect relation between vmin "
            "and vmax. The vmin value cannot be "
            "bigger than or equal to the value "
            "of vmax."
        )
    if len(colormap) == 1:
        # color each triangle face with the same color in colormap
        face_color = colormap[0]
        face_color = clrs.convert_to_RGB_255(face_color)
        face_color = clrs.label_rgb(face_color)
        return face_color
    if face == vmax:
        # pick last color in colormap
        face_color = colormap[-1]
        face_color = clrs.convert_to_RGB_255(face_color)
        face_color = clrs.label_rgb(face_color)
        return face_color
    else:
        if scale is None:
            # find the normalized distance t of a triangle face between
            # vmin and vmax where the distance is between 0 and 1
            t = (face - vmin) / float((vmax - vmin))
            low_color_index = int(t / (1.0 / (len(colormap) - 1)))

            face_color = clrs.find_intermediate_color(
                colormap[low_color_index],
                colormap[low_color_index + 1],
                t * (len(colormap) - 1) - low_color_index,
            )

            face_color = clrs.convert_to_RGB_255(face_color)
            face_color = clrs.label_rgb(face_color)
        else:
            # find the face color for a non-linearly interpolated scale
            t = (face - vmin) / float((vmax - vmin))

            low_color_index = 0
            for k in range(len(scale) - 1):
                if scale[k] <= t < scale[k + 1]:
                    break
                low_color_index += 1

            low_scale_val = scale[low_color_index]
            high_scale_val = scale[low_color_index + 1]

            face_color = clrs.find_intermediate_color(
                colormap[low_color_index],
                colormap[low_color_index + 1],
                (t - low_scale_val) / (high_scale_val - low_scale_val),
            )

            face_color = clrs.convert_to_RGB_255(face_color)
            face_color = clrs.label_rgb(face_color)
        return face_color


def trisurf(
    x,
    y,
    z,
    simplices,
    show_colorbar,
    edges_color,
    scale,
    colormap=None,
    color_func=None,
    plot_edges=False,
    x_edge=None,
    y_edge=None,
    z_edge=None,
    facecolor=None,
):
    """
    Refer to FigureFactory.create_trisurf() for docstring
    """
    # numpy import check
    if not np:
        raise ImportError("FigureFactory._trisurf() requires " "numpy imported.")
    points3D = np.vstack((x, y, z)).T
    simplices = np.atleast_2d(simplices)

    # vertices of the surface triangles
    tri_vertices = points3D[simplices]

    # Define colors for the triangle faces
    if color_func is None:
        # mean values of z-coordinates of triangle vertices
        mean_dists = tri_vertices[:, :, 2].mean(-1)
    elif isinstance(color_func, (list, np.ndarray)):
        # Pre-computed list / array of values to map onto color
        if len(color_func) != len(simplices):
            raise ValueError(
                "If color_func is a list/array, it must "
                "be the same length as simplices."
            )

        # convert all colors in color_func to rgb
        for index in range(len(color_func)):
            if isinstance(color_func[index], str):
                if "#" in color_func[index]:
                    foo = clrs.hex_to_rgb(color_func[index])
                    color_func[index] = clrs.label_rgb(foo)

            if isinstance(color_func[index], tuple):
                foo = clrs.convert_to_RGB_255(color_func[index])
                color_func[index] = clrs.label_rgb(foo)

        mean_dists = np.asarray(color_func)
    else:
        # apply user inputted function to calculate
        # custom coloring for triangle vertices
        mean_dists = []
        for triangle in tri_vertices:
            dists = []
            for vertex in triangle:
                dist = color_func(vertex[0], vertex[1], vertex[2])
                dists.append(dist)
            mean_dists.append(np.mean(dists))
        mean_dists = np.asarray(mean_dists)

    # Check if facecolors are already strings and can be skipped
    if isinstance(mean_dists[0], str):
        facecolor = mean_dists
    else:
        min_mean_dists = np.min(mean_dists)
        max_mean_dists = np.max(mean_dists)

        if facecolor is None:
            facecolor = []
        for index in range(len(mean_dists)):
            color = map_face2color(
                mean_dists[index], colormap, scale, min_mean_dists, max_mean_dists
            )
            facecolor.append(color)

    # Make sure facecolor is a list so output is consistent across Pythons
    facecolor = np.asarray(facecolor)
    ii, jj, kk = simplices.T

    triangles = graph_objs.Mesh3d(
        x=x, y=y, z=z, facecolor=facecolor, i=ii, j=jj, k=kk, name=""
    )

    mean_dists_are_numbers = not isinstance(mean_dists[0], str)

    if mean_dists_are_numbers and show_colorbar is True:
        # make a colorscale from the colors
        colorscale = clrs.make_colorscale(colormap, scale)
        colorscale = clrs.convert_colorscale_to_rgb(colorscale)

        colorbar = graph_objs.Scatter3d(
            x=x[:1],
            y=y[:1],
            z=z[:1],
            mode="markers",
            marker=dict(
                size=0.1,
                color=[min_mean_dists, max_mean_dists],
                colorscale=colorscale,
                showscale=True,
            ),
            hoverinfo="none",
            showlegend=False,
        )

    # the triangle sides are not plotted
    if plot_edges is False:
        if mean_dists_are_numbers and show_colorbar is True:
            return [triangles, colorbar]
        else:
            return [triangles]

    # define the lists x_edge, y_edge and z_edge, of x, y, resp z
    # coordinates of edge end points for each triangle
    # None separates data corresponding to two consecutive triangles
    is_none = [ii is None for ii in [x_edge, y_edge, z_edge]]
    if any(is_none):
        if not all(is_none):
            raise ValueError(
                "If any (x_edge, y_edge, z_edge) is None, " "all must be None"
            )
        else:
            x_edge = []
            y_edge = []
            z_edge = []

    # Pull indices we care about, then add a None column to separate tris
    ixs_triangles = [0, 1, 2, 0]
    pull_edges = tri_vertices[:, ixs_triangles, :]
    x_edge_pull = np.hstack(
        [pull_edges[:, :, 0], np.tile(None, [pull_edges.shape[0], 1])]
    )
    y_edge_pull = np.hstack(
        [pull_edges[:, :, 1], np.tile(None, [pull_edges.shape[0], 1])]
    )
    z_edge_pull = np.hstack(
        [pull_edges[:, :, 2], np.tile(None, [pull_edges.shape[0], 1])]
    )

    # Now unravel the edges into a 1-d vector for plotting
    x_edge = np.hstack([x_edge, x_edge_pull.reshape([1, -1])[0]])
    y_edge = np.hstack([y_edge, y_edge_pull.reshape([1, -1])[0]])
    z_edge = np.hstack([z_edge, z_edge_pull.reshape([1, -1])[0]])

    if not (len(x_edge) == len(y_edge) == len(z_edge)):
        raise exceptions.PlotlyError(
            "The lengths of x_edge, y_edge and " "z_edge are not the same."
        )

    # define the lines for plotting
    lines = graph_objs.Scatter3d(
        x=x_edge,
        y=y_edge,
        z=z_edge,
        mode="lines",
        line=graph_objs.scatter3d.Line(color=edges_color, width=1.5),
        showlegend=False,
    )

    if mean_dists_are_numbers and show_colorbar is True:
        return [triangles, lines, colorbar]
    else:
        return [triangles, lines]


def create_trisurf(
    x,
    y,
    z,
    simplices,
    colormap=None,
    show_colorbar=True,
    scale=None,
    color_func=None,
    title="Trisurf Plot",
    plot_edges=True,
    showbackground=True,
    backgroundcolor="rgb(230, 230, 230)",
    gridcolor="rgb(255, 255, 255)",
    zerolinecolor="rgb(255, 255, 255)",
    edges_color="rgb(50, 50, 50)",
    height=800,
    width=800,
    aspectratio=None,
):
    """
    Returns figure for a triangulated surface plot

    :param (array) x: data values of x in a 1D array
    :param (array) y: data values of y in a 1D array
    :param (array) z: data values of z in a 1D array
    :param (array) simplices: an array of shape (ntri, 3) where ntri is
        the number of triangles in the triangularization. Each row of the
        array contains the indicies of the verticies of each triangle
    :param (str|tuple|list) colormap: either a plotly scale name, an rgb
        or hex color, a color tuple or a list of colors. An rgb color is
        of the form 'rgb(x, y, z)' where x, y, z belong to the interval
        [0, 255] and a color tuple is a tuple of the form (a, b, c) where
        a, b and c belong to [0, 1]. If colormap is a list, it must
        contain the valid color types aforementioned as its members
    :param (bool) show_colorbar: determines if colorbar is visible
    :param (list|array) scale: sets the scale values to be used if a non-
        linearly interpolated colormap is desired. If left as None, a
        linear interpolation between the colors will be excecuted
    :param (function|list) color_func: The parameter that determines the
        coloring of the surface. Takes either a function with 3 arguments
        x, y, z or a list/array of color values the same length as
        simplices. If None, coloring will only depend on the z axis
    :param (str) title: title of the plot
    :param (bool) plot_edges: determines if the triangles on the trisurf
        are visible
    :param (bool) showbackground: makes background in plot visible
    :param (str) backgroundcolor: color of background. Takes a string of
        the form 'rgb(x,y,z)' x,y,z are between 0 and 255 inclusive
    :param (str) gridcolor: color of the gridlines besides the axes. Takes
        a string of the form 'rgb(x,y,z)' x,y,z are between 0 and 255
        inclusive
    :param (str) zerolinecolor: color of the axes. Takes a string of the
        form 'rgb(x,y,z)' x,y,z are between 0 and 255 inclusive
    :param (str) edges_color: color of the edges, if plot_edges is True
    :param (int|float) height: the height of the plot (in pixels)
    :param (int|float) width: the width of the plot (in pixels)
    :param (dict) aspectratio: a dictionary of the aspect ratio values for
        the x, y and z axes. 'x', 'y' and 'z' take (int|float) values

    Example 1: Sphere

    >>> # Necessary Imports for Trisurf
    >>> import numpy as np
    >>> from scipy.spatial import Delaunay

    >>> from plotly.figure_factory import create_trisurf
    >>> from plotly.graph_objs import graph_objs

    >>> # Make data for plot
    >>> u = np.linspace(0, 2*np.pi, 20)
    >>> v = np.linspace(0, np.pi, 20)
    >>> u,v = np.meshgrid(u,v)
    >>> u = u.flatten()
    >>> v = v.flatten()

    >>> x = np.sin(v)*np.cos(u)
    >>> y = np.sin(v)*np.sin(u)
    >>> z = np.cos(v)

    >>> points2D = np.vstack([u,v]).T
    >>> tri = Delaunay(points2D)
    >>> simplices = tri.simplices

    >>> # Create a figure
    >>> fig1 = create_trisurf(x=x, y=y, z=z, colormap="Rainbow",
    ...                       simplices=simplices)

    Example 2: Torus

    >>> # Necessary Imports for Trisurf
    >>> import numpy as np
    >>> from scipy.spatial import Delaunay

    >>> from plotly.figure_factory import create_trisurf
    >>> from plotly.graph_objs import graph_objs

    >>> # Make data for plot
    >>> u = np.linspace(0, 2*np.pi, 20)
    >>> v = np.linspace(0, 2*np.pi, 20)
    >>> u,v = np.meshgrid(u,v)
    >>> u = u.flatten()
    >>> v = v.flatten()

    >>> x = (3 + (np.cos(v)))*np.cos(u)
    >>> y = (3 + (np.cos(v)))*np.sin(u)
    >>> z = np.sin(v)

    >>> points2D = np.vstack([u,v]).T
    >>> tri = Delaunay(points2D)
    >>> simplices = tri.simplices

    >>> # Create a figure
    >>> fig1 = create_trisurf(x=x, y=y, z=z, colormap="Viridis",
    ...                       simplices=simplices)

    Example 3: Mobius Band

    >>> # Necessary Imports for Trisurf
    >>> import numpy as np
    >>> from scipy.spatial import Delaunay

    >>> from plotly.figure_factory import create_trisurf
    >>> from plotly.graph_objs import graph_objs

    >>> # Make data for plot
    >>> u = np.linspace(0, 2*np.pi, 24)
    >>> v = np.linspace(-1, 1, 8)
    >>> u,v = np.meshgrid(u,v)
    >>> u = u.flatten()
    >>> v = v.flatten()

    >>> tp = 1 + 0.5*v*np.cos(u/2.)
    >>> x = tp*np.cos(u)
    >>> y = tp*np.sin(u)
    >>> z = 0.5*v*np.sin(u/2.)

    >>> points2D = np.vstack([u,v]).T
    >>> tri = Delaunay(points2D)
    >>> simplices = tri.simplices

    >>> # Create a figure
    >>> fig1 = create_trisurf(x=x, y=y, z=z, colormap=[(0.2, 0.4, 0.6), (1, 1, 1)],
    ...                       simplices=simplices)

    Example 4: Using a Custom Colormap Function with Light Cone

    >>> # Necessary Imports for Trisurf
    >>> import numpy as np
    >>> from scipy.spatial import Delaunay

    >>> from plotly.figure_factory import create_trisurf
    >>> from plotly.graph_objs import graph_objs

    >>> # Make data for plot
    >>> u=np.linspace(-np.pi, np.pi, 30)
    >>> v=np.linspace(-np.pi, np.pi, 30)
    >>> u,v=np.meshgrid(u,v)
    >>> u=u.flatten()
    >>> v=v.flatten()

    >>> x = u
    >>> y = u*np.cos(v)
    >>> z = u*np.sin(v)

    >>> points2D = np.vstack([u,v]).T
    >>> tri = Delaunay(points2D)
    >>> simplices = tri.simplices

    >>> # Define distance function
    >>> def dist_origin(x, y, z):
    ...     return np.sqrt((1.0 * x)**2 + (1.0 * y)**2 + (1.0 * z)**2)

    >>> # Create a figure
    >>> fig1 = create_trisurf(x=x, y=y, z=z,
    ...                       colormap=['#FFFFFF', '#E4FFFE',
    ...                                 '#A4F6F9', '#FF99FE',
    ...                                 '#BA52ED'],
    ...                       scale=[0, 0.6, 0.71, 0.89, 1],
    ...                       simplices=simplices,
    ...                       color_func=dist_origin)

    Example 5: Enter color_func as a list of colors

    >>> # Necessary Imports for Trisurf
    >>> import numpy as np
    >>> from scipy.spatial import Delaunay
    >>> import random

    >>> from plotly.figure_factory import create_trisurf
    >>> from plotly.graph_objs import graph_objs

    >>> # Make data for plot
    >>> u=np.linspace(-np.pi, np.pi, 30)
    >>> v=np.linspace(-np.pi, np.pi, 30)
    >>> u,v=np.meshgrid(u,v)
    >>> u=u.flatten()
    >>> v=v.flatten()

    >>> x = u
    >>> y = u*np.cos(v)
    >>> z = u*np.sin(v)

    >>> points2D = np.vstack([u,v]).T
    >>> tri = Delaunay(points2D)
    >>> simplices = tri.simplices


    >>> colors = []
    >>> color_choices = ['rgb(0, 0, 0)', '#6c4774', '#d6c7dd']

    >>> for index in range(len(simplices)):
    ...     colors.append(random.choice(color_choices))

    >>> fig = create_trisurf(
    ...     x, y, z, simplices,
    ...     color_func=colors,
    ...     show_colorbar=True,
    ...     edges_color='rgb(2, 85, 180)',
    ...     title=' Modern Art'
    ... )
    """
    if aspectratio is None:
        aspectratio = {"x": 1, "y": 1, "z": 1}

    # Validate colormap
    clrs.validate_colors(colormap)
    colormap, scale = clrs.convert_colors_to_same_type(
        colormap, colortype="tuple", return_default_colors=True, scale=scale
    )

    data1 = trisurf(
        x,
        y,
        z,
        simplices,
        show_colorbar=show_colorbar,
        color_func=color_func,
        colormap=colormap,
        scale=scale,
        edges_color=edges_color,
        plot_edges=plot_edges,
    )

    axis = dict(
        showbackground=showbackground,
        backgroundcolor=backgroundcolor,
        gridcolor=gridcolor,
        zerolinecolor=zerolinecolor,
    )
    layout = graph_objs.Layout(
        title=title,
        width=width,
        height=height,
        scene=graph_objs.layout.Scene(
            xaxis=graph_objs.layout.scene.XAxis(**axis),
            yaxis=graph_objs.layout.scene.YAxis(**axis),
            zaxis=graph_objs.layout.scene.ZAxis(**axis),
            aspectratio=dict(
                x=aspectratio["x"], y=aspectratio["y"], z=aspectratio["z"]
            ),
        ),
    )

    return graph_objs.Figure(data=data1, layout=layout)
