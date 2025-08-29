from plotly.utils import _list_repr_elided


class InputDeviceState:
    def __init__(
        self, ctrl=None, alt=None, shift=None, meta=None, button=None, buttons=None, **_
    ):
        self._ctrl = ctrl
        self._alt = alt
        self._meta = meta
        self._shift = shift
        self._button = button
        self._buttons = buttons

    def __repr__(self):
        return """\
InputDeviceState(
    ctrl={ctrl},
    alt={alt},
    shift={shift},
    meta={meta},
    button={button},
    buttons={buttons})""".format(
            ctrl=repr(self.ctrl),
            alt=repr(self.alt),
            meta=repr(self.meta),
            shift=repr(self.shift),
            button=repr(self.button),
            buttons=repr(self.buttons),
        )

    @property
    def alt(self):
        """
        Whether alt key pressed

        Returns
        -------
        bool
        """
        return self._alt

    @property
    def ctrl(self):
        """
        Whether ctrl key pressed

        Returns
        -------
        bool
        """
        return self._ctrl

    @property
    def shift(self):
        """
        Whether shift key pressed

        Returns
        -------
        bool
        """
        return self._shift

    @property
    def meta(self):
        """
        Whether meta key pressed

        Returns
        -------
        bool
        """
        return self._meta

    @property
    def button(self):
        """
        Integer code for the button that was pressed on the mouse to trigger
        the event

        - 0: Main button pressed, usually the left button or the
             un-initialized state
        - 1: Auxiliary button pressed, usually the wheel button or the middle
             button (if present)
        - 2: Secondary button pressed, usually the right button
        - 3: Fourth button, typically the Browser Back button
        - 4: Fifth button, typically the Browser Forward button

        Returns
        -------
        int
        """
        return self._button

    @property
    def buttons(self):
        """
        Integer code for which combination of buttons are pressed on the
        mouse when the event is triggered.

        -  0: No button or un-initialized
        -  1: Primary button (usually left)
        -  2: Secondary button (usually right)
        -  4: Auxilary button (usually middle or mouse wheel button)
        -  8: 4th button (typically the "Browser Back" button)
        - 16: 5th button (typically the "Browser Forward" button)

        Combinations of buttons are represented as the decimal form of the
        bitmask of the values above.

        For example, pressing both the primary (1) and auxilary (4) buttons
        will result in a code of 5

        Returns
        -------
        int
        """
        return self._buttons


class Points:
    def __init__(self, point_inds=[], xs=[], ys=[], trace_name=None, trace_index=None):
        self._point_inds = point_inds
        self._xs = xs
        self._ys = ys
        self._trace_name = trace_name
        self._trace_index = trace_index

    def __repr__(self):
        return """\
Points(point_inds={point_inds},
       xs={xs},
       ys={ys},
       trace_name={trace_name},
       trace_index={trace_index})""".format(
            point_inds=_list_repr_elided(
                self.point_inds, indent=len("Points(point_inds=")
            ),
            xs=_list_repr_elided(self.xs, indent=len("       xs=")),
            ys=_list_repr_elided(self.ys, indent=len("       ys=")),
            trace_name=repr(self.trace_name),
            trace_index=repr(self.trace_index),
        )

    @property
    def point_inds(self):
        """
        List of selected indexes into the trace's points

        Returns
        -------
        list[int]
        """
        return self._point_inds

    @property
    def xs(self):
        """
        List of x-coordinates of selected points

        Returns
        -------
        list[float]
        """
        return self._xs

    @property
    def ys(self):
        """
        List of y-coordinates of selected points

        Returns
        -------
        list[float]
        """
        return self._ys

    @property
    def trace_name(self):
        """
        Name of the trace

        Returns
        -------
        str
        """
        return self._trace_name

    @property
    def trace_index(self):
        """
        Index of the trace in the figure

        Returns
        -------
        int
        """
        return self._trace_index


class BoxSelector:
    def __init__(self, xrange=None, yrange=None, **_):
        self._type = "box"
        self._xrange = xrange
        self._yrange = yrange

    def __repr__(self):
        return """\
BoxSelector(xrange={xrange},
            yrange={yrange})""".format(xrange=self.xrange, yrange=self.yrange)

    @property
    def type(self):
        """
        The selector's type

        Returns
        -------
        str
        """
        return self._type

    @property
    def xrange(self):
        """
        x-axis range extents of the box selection

        Returns
        -------
        (float, float)
        """
        return self._xrange

    @property
    def yrange(self):
        """
        y-axis range extents of the box selection

        Returns
        -------
        (float, float)
        """
        return self._yrange


class LassoSelector:
    def __init__(self, xs=None, ys=None, **_):
        self._type = "lasso"
        self._xs = xs
        self._ys = ys

    def __repr__(self):
        return """\
LassoSelector(xs={xs},
              ys={ys})""".format(
            xs=_list_repr_elided(self.xs, indent=len("LassoSelector(xs=")),
            ys=_list_repr_elided(self.ys, indent=len("              ys=")),
        )

    @property
    def type(self):
        """
        The selector's type

        Returns
        -------
        str
        """
        return self._type

    @property
    def xs(self):
        """
        list of x-axis coordinates of each point in the lasso selection
        boundary

        Returns
        -------
        list[float]
        """
        return self._xs

    @property
    def ys(self):
        """
        list of y-axis coordinates of each point in the lasso selection
        boundary

        Returns
        -------
        list[float]
        """
        return self._ys
