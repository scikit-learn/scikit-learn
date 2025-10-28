#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Stream(_BaseTraceHierarchyType):
    _parent_path_str = "barpolar"
    _path_str = "barpolar.stream"
    _valid_props = {"maxpoints", "token"}

    @property
    def maxpoints(self):
        """
        Sets the maximum number of points to keep on the plots from an
        incoming stream. If `maxpoints` is set to 50, only the newest
        50 points will be displayed on the plot.

        The 'maxpoints' property is a number and may be specified as:
          - An int or float in the interval [0, 10000]

        Returns
        -------
        int|float
        """
        return self["maxpoints"]

    @maxpoints.setter
    def maxpoints(self, val):
        self["maxpoints"] = val

    @property
    def token(self):
        """
        The stream id number links a data trace on a plot with a
        stream. See https://chart-studio.plotly.com/settings for more
        details.

        The 'token' property is a string and must be specified as:
          - A non-empty string

        Returns
        -------
        str
        """
        return self["token"]

    @token.setter
    def token(self, val):
        self["token"] = val

    @property
    def _prop_descriptions(self):
        return """\
        maxpoints
            Sets the maximum number of points to keep on the plots
            from an incoming stream. If `maxpoints` is set to 50,
            only the newest 50 points will be displayed on the
            plot.
        token
            The stream id number links a data trace on a plot with
            a stream. See https://chart-studio.plotly.com/settings
            for more details.
        """

    def __init__(self, arg=None, maxpoints=None, token=None, **kwargs):
        """
        Construct a new Stream object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.barpolar.Stream`
        maxpoints
            Sets the maximum number of points to keep on the plots
            from an incoming stream. If `maxpoints` is set to 50,
            only the newest 50 points will be displayed on the
            plot.
        token
            The stream id number links a data trace on a plot with
            a stream. See https://chart-studio.plotly.com/settings
            for more details.

        Returns
        -------
        Stream
        """
        super().__init__("stream")
        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError("""\
The first argument to the plotly.graph_objs.barpolar.Stream
constructor must be a dict or
an instance of :class:`plotly.graph_objs.barpolar.Stream`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("maxpoints", arg, maxpoints)
        self._set_property("token", arg, token)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
