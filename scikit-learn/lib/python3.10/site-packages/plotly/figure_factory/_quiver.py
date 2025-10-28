import math

from plotly import exceptions
from plotly.graph_objs import graph_objs
from plotly.figure_factory import utils


def create_quiver(
    x, y, u, v, scale=0.1, arrow_scale=0.3, angle=math.pi / 9, scaleratio=None, **kwargs
):
    """
    Returns data for a quiver plot.

    :param (list|ndarray) x: x coordinates of the arrow locations
    :param (list|ndarray) y: y coordinates of the arrow locations
    :param (list|ndarray) u: x components of the arrow vectors
    :param (list|ndarray) v: y components of the arrow vectors
    :param (float in [0,1]) scale: scales size of the arrows(ideally to
        avoid overlap). Default = .1
    :param (float in [0,1]) arrow_scale: value multiplied to length of barb
        to get length of arrowhead. Default = .3
    :param (angle in radians) angle: angle of arrowhead. Default = pi/9
    :param (positive float) scaleratio: the ratio between the scale of the y-axis
        and the scale of the x-axis (scale_y / scale_x). Default = None, the
        scale ratio is not fixed.
    :param kwargs: kwargs passed through plotly.graph_objs.Scatter
        for more information on valid kwargs call
        help(plotly.graph_objs.Scatter)

    :rtype (dict): returns a representation of quiver figure.

    Example 1: Trivial Quiver

    >>> from plotly.figure_factory import create_quiver
    >>> import math

    >>> # 1 Arrow from (0,0) to (1,1)
    >>> fig = create_quiver(x=[0], y=[0], u=[1], v=[1], scale=1)
    >>> fig.show()


    Example 2: Quiver plot using meshgrid

    >>> from plotly.figure_factory import create_quiver

    >>> import numpy as np
    >>> import math

    >>> # Add data
    >>> x,y = np.meshgrid(np.arange(0, 2, .2), np.arange(0, 2, .2))
    >>> u = np.cos(x)*y
    >>> v = np.sin(x)*y

    >>> #Create quiver
    >>> fig = create_quiver(x, y, u, v)
    >>> fig.show()


    Example 3: Styling the quiver plot

    >>> from plotly.figure_factory import create_quiver
    >>> import numpy as np
    >>> import math

    >>> # Add data
    >>> x, y = np.meshgrid(np.arange(-np.pi, math.pi, .5),
    ...                    np.arange(-math.pi, math.pi, .5))
    >>> u = np.cos(x)*y
    >>> v = np.sin(x)*y

    >>> # Create quiver
    >>> fig = create_quiver(x, y, u, v, scale=.2, arrow_scale=.3, angle=math.pi/6,
    ...                     name='Wind Velocity', line=dict(width=1))

    >>> # Add title to layout
    >>> fig.update_layout(title='Quiver Plot') # doctest: +SKIP
    >>> fig.show()


    Example 4: Forcing a fix scale ratio to maintain the arrow length

    >>> from plotly.figure_factory import create_quiver
    >>> import numpy as np

    >>> # Add data
    >>> x,y = np.meshgrid(np.arange(0.5, 3.5, .5), np.arange(0.5, 4.5, .5))
    >>> u = x
    >>> v = y
    >>> angle = np.arctan(v / u)
    >>> norm = 0.25
    >>> u = norm * np.cos(angle)
    >>> v = norm * np.sin(angle)

    >>> # Create quiver with a fix scale ratio
    >>> fig = create_quiver(x, y, u, v, scale = 1, scaleratio = 0.5)
    >>> fig.show()
    """
    utils.validate_equal_length(x, y, u, v)
    utils.validate_positive_scalars(arrow_scale=arrow_scale, scale=scale)

    if scaleratio is None:
        quiver_obj = _Quiver(x, y, u, v, scale, arrow_scale, angle)
    else:
        quiver_obj = _Quiver(x, y, u, v, scale, arrow_scale, angle, scaleratio)

    barb_x, barb_y = quiver_obj.get_barbs()
    arrow_x, arrow_y = quiver_obj.get_quiver_arrows()

    quiver_plot = graph_objs.Scatter(
        x=barb_x + arrow_x, y=barb_y + arrow_y, mode="lines", **kwargs
    )

    data = [quiver_plot]

    if scaleratio is None:
        layout = graph_objs.Layout(hovermode="closest")
    else:
        layout = graph_objs.Layout(
            hovermode="closest", yaxis=dict(scaleratio=scaleratio, scaleanchor="x")
        )

    return graph_objs.Figure(data=data, layout=layout)


class _Quiver(object):
    """
    Refer to FigureFactory.create_quiver() for docstring
    """

    def __init__(self, x, y, u, v, scale, arrow_scale, angle, scaleratio=1, **kwargs):
        try:
            x = utils.flatten(x)
        except exceptions.PlotlyError:
            pass

        try:
            y = utils.flatten(y)
        except exceptions.PlotlyError:
            pass

        try:
            u = utils.flatten(u)
        except exceptions.PlotlyError:
            pass

        try:
            v = utils.flatten(v)
        except exceptions.PlotlyError:
            pass

        self.x = x
        self.y = y
        self.u = u
        self.v = v
        self.scale = scale
        self.scaleratio = scaleratio
        self.arrow_scale = arrow_scale
        self.angle = angle
        self.end_x = []
        self.end_y = []
        self.scale_uv()
        barb_x, barb_y = self.get_barbs()
        arrow_x, arrow_y = self.get_quiver_arrows()

    def scale_uv(self):
        """
        Scales u and v to avoid overlap of the arrows.

        u and v are added to x and y to get the
        endpoints of the arrows so a smaller scale value will
        result in less overlap of arrows.
        """
        self.u = [i * self.scale * self.scaleratio for i in self.u]
        self.v = [i * self.scale for i in self.v]

    def get_barbs(self):
        """
        Creates x and y startpoint and endpoint pairs

        After finding the endpoint of each barb this zips startpoint and
        endpoint pairs to create 2 lists: x_values for barbs and y values
        for barbs

        :rtype: (list, list) barb_x, barb_y: list of startpoint and endpoint
            x_value pairs separated by a None to create the barb of the arrow,
            and list of startpoint and endpoint y_value pairs separated by a
            None to create the barb of the arrow.
        """
        self.end_x = [i + j for i, j in zip(self.x, self.u)]
        self.end_y = [i + j for i, j in zip(self.y, self.v)]
        empty = [None] * len(self.x)
        barb_x = utils.flatten(zip(self.x, self.end_x, empty))
        barb_y = utils.flatten(zip(self.y, self.end_y, empty))
        return barb_x, barb_y

    def get_quiver_arrows(self):
        """
        Creates lists of x and y values to plot the arrows

        Gets length of each barb then calculates the length of each side of
        the arrow. Gets angle of barb and applies angle to each side of the
        arrowhead. Next uses arrow_scale to scale the length of arrowhead and
        creates x and y values for arrowhead point1 and point2. Finally x and y
        values for point1, endpoint and point2s for each arrowhead are
        separated by a None and zipped to create lists of x and y values for
        the arrows.

        :rtype: (list, list) arrow_x, arrow_y: list of point1, endpoint, point2
            x_values separated by a None to create the arrowhead and list of
            point1, endpoint, point2 y_values separated by a None to create
            the barb of the arrow.
        """
        dif_x = [i - j for i, j in zip(self.end_x, self.x)]
        dif_y = [i - j for i, j in zip(self.end_y, self.y)]

        # Get barb lengths(default arrow length = 30% barb length)
        barb_len = [None] * len(self.x)
        for index in range(len(barb_len)):
            barb_len[index] = math.hypot(dif_x[index] / self.scaleratio, dif_y[index])

        # Make arrow lengths
        arrow_len = [None] * len(self.x)
        arrow_len = [i * self.arrow_scale for i in barb_len]

        # Get barb angles
        barb_ang = [None] * len(self.x)
        for index in range(len(barb_ang)):
            barb_ang[index] = math.atan2(dif_y[index], dif_x[index] / self.scaleratio)

        # Set angles to create arrow
        ang1 = [i + self.angle for i in barb_ang]
        ang2 = [i - self.angle for i in barb_ang]

        cos_ang1 = [None] * len(ang1)
        for index in range(len(ang1)):
            cos_ang1[index] = math.cos(ang1[index])
        seg1_x = [i * j for i, j in zip(arrow_len, cos_ang1)]

        sin_ang1 = [None] * len(ang1)
        for index in range(len(ang1)):
            sin_ang1[index] = math.sin(ang1[index])
        seg1_y = [i * j for i, j in zip(arrow_len, sin_ang1)]

        cos_ang2 = [None] * len(ang2)
        for index in range(len(ang2)):
            cos_ang2[index] = math.cos(ang2[index])
        seg2_x = [i * j for i, j in zip(arrow_len, cos_ang2)]

        sin_ang2 = [None] * len(ang2)
        for index in range(len(ang2)):
            sin_ang2[index] = math.sin(ang2[index])
        seg2_y = [i * j for i, j in zip(arrow_len, sin_ang2)]

        # Set coordinates to create arrow
        for index in range(len(self.end_x)):
            point1_x = [i - j * self.scaleratio for i, j in zip(self.end_x, seg1_x)]
            point1_y = [i - j for i, j in zip(self.end_y, seg1_y)]
            point2_x = [i - j * self.scaleratio for i, j in zip(self.end_x, seg2_x)]
            point2_y = [i - j for i, j in zip(self.end_y, seg2_y)]

        # Combine lists to create arrow
        empty = [None] * len(self.end_x)
        arrow_x = utils.flatten(zip(point1_x, self.end_x, point2_x, empty))
        arrow_y = utils.flatten(zip(point1_y, self.end_y, point2_y, empty))
        return arrow_x, arrow_y
