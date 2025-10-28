#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Selection(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.selection"
    _valid_props = {
        "line",
        "name",
        "opacity",
        "path",
        "templateitemname",
        "type",
        "x0",
        "x1",
        "xref",
        "y0",
        "y1",
        "yref",
    }

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.selection.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.layout.selection.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def name(self):
        """
        When used in a template, named items are created in the output
        figure in addition to any items the figure already has in this
        array. You can modify these items in the output figure by
        making your own item with `templateitemname` matching this
        `name` alongside your modifications (including `visible: false`
        or `enabled: false` to hide it). Has no effect outside of a
        template.

        The 'name' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["name"]

    @name.setter
    def name(self, val):
        self["name"] = val

    @property
    def opacity(self):
        """
        Sets the opacity of the selection.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def path(self):
        """
        For `type` "path" - a valid SVG path similar to `shapes.path`
        in data coordinates. Allowed segments are: M, L and Z.

        The 'path' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["path"]

    @path.setter
    def path(self, val):
        self["path"] = val

    @property
    def templateitemname(self):
        """
        Used to refer to a named item in this array in the template.
        Named items from the template will be created even without a
        matching item in the input figure, but you can modify one by
        making an item with `templateitemname` matching its `name`,
        alongside your modifications (including `visible: false` or
        `enabled: false` to hide it). If there is no template or no
        matching item, this item will be hidden unless you explicitly
        show it with `visible: true`.

        The 'templateitemname' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["templateitemname"]

    @templateitemname.setter
    def templateitemname(self, val):
        self["templateitemname"] = val

    @property
    def type(self):
        """
        Specifies the selection type to be drawn. If "rect", a
        rectangle is drawn linking (`x0`,`y0`), (`x1`,`y0`),
        (`x1`,`y1`) and (`x0`,`y1`). If "path", draw a custom SVG path
        using `path`.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['rect', 'path']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    @property
    def x0(self):
        """
        Sets the selection's starting x position.

        The 'x0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["x0"]

    @x0.setter
    def x0(self, val):
        self["x0"] = val

    @property
    def x1(self):
        """
        Sets the selection's end x position.

        The 'x1' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["x1"]

    @x1.setter
    def x1(self, val):
        self["x1"] = val

    @property
    def xref(self):
        """
        Sets the selection's x coordinate axis. If set to a x axis id
        (e.g. "x" or "x2"), the `x` position refers to a x coordinate.
        If set to "paper", the `x` position refers to the distance from
        the left of the plotting area in normalized coordinates where 0
        (1) corresponds to the left (right). If set to a x axis ID
        followed by "domain" (separated by a space), the position
        behaves like for "paper", but refers to the distance in
        fractions of the domain length from the left of the domain of
        that axis: e.g., *x2 domain* refers to the domain of the second
        x  axis and a x position of 0.5 refers to the point between the
        left and the right of the domain of the second x axis.

        The 'xref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['paper']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["xref"]

    @xref.setter
    def xref(self, val):
        self["xref"] = val

    @property
    def y0(self):
        """
        Sets the selection's starting y position.

        The 'y0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["y0"]

    @y0.setter
    def y0(self, val):
        self["y0"] = val

    @property
    def y1(self):
        """
        Sets the selection's end y position.

        The 'y1' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["y1"]

    @y1.setter
    def y1(self, val):
        self["y1"] = val

    @property
    def yref(self):
        """
        Sets the selection's x coordinate axis. If set to a y axis id
        (e.g. "y" or "y2"), the `y` position refers to a y coordinate.
        If set to "paper", the `y` position refers to the distance from
        the bottom of the plotting area in normalized coordinates where
        0 (1) corresponds to the bottom (top). If set to a y axis ID
        followed by "domain" (separated by a space), the position
        behaves like for "paper", but refers to the distance in
        fractions of the domain length from the bottom of the domain of
        that axis: e.g., *y2 domain* refers to the domain of the second
        y  axis and a y position of 0.5 refers to the point between the
        bottom and the top of the domain of the second y axis.

        The 'yref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['paper']
          - A string that matches one of the following regular expressions:
                ['^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["yref"]

    @yref.setter
    def yref(self, val):
        self["yref"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.layout.selection.Line`
            instance or dict with compatible properties
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the selection.
        path
            For `type` "path" - a valid SVG path similar to
            `shapes.path` in data coordinates. Allowed segments
            are: M, L and Z.
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        type
            Specifies the selection type to be drawn. If "rect", a
            rectangle is drawn linking (`x0`,`y0`), (`x1`,`y0`),
            (`x1`,`y1`) and (`x0`,`y1`). If "path", draw a custom
            SVG path using `path`.
        x0
            Sets the selection's starting x position.
        x1
            Sets the selection's end x position.
        xref
            Sets the selection's x coordinate axis. If set to a x
            axis id (e.g. "x" or "x2"), the `x` position refers to
            a x coordinate. If set to "paper", the `x` position
            refers to the distance from the left of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        y0
            Sets the selection's starting y position.
        y1
            Sets the selection's end y position.
        yref
            Sets the selection's x coordinate axis. If set to a y
            axis id (e.g. "y" or "y2"), the `y` position refers to
            a y coordinate. If set to "paper", the `y` position
            refers to the distance from the bottom of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis.
        """

    def __init__(
        self,
        arg=None,
        line=None,
        name=None,
        opacity=None,
        path=None,
        templateitemname=None,
        type=None,
        x0=None,
        x1=None,
        xref=None,
        y0=None,
        y1=None,
        yref=None,
        **kwargs,
    ):
        """
        Construct a new Selection object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Selection`
        line
            :class:`plotly.graph_objects.layout.selection.Line`
            instance or dict with compatible properties
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the selection.
        path
            For `type` "path" - a valid SVG path similar to
            `shapes.path` in data coordinates. Allowed segments
            are: M, L and Z.
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        type
            Specifies the selection type to be drawn. If "rect", a
            rectangle is drawn linking (`x0`,`y0`), (`x1`,`y0`),
            (`x1`,`y1`) and (`x0`,`y1`). If "path", draw a custom
            SVG path using `path`.
        x0
            Sets the selection's starting x position.
        x1
            Sets the selection's end x position.
        xref
            Sets the selection's x coordinate axis. If set to a x
            axis id (e.g. "x" or "x2"), the `x` position refers to
            a x coordinate. If set to "paper", the `x` position
            refers to the distance from the left of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        y0
            Sets the selection's starting y position.
        y1
            Sets the selection's end y position.
        yref
            Sets the selection's x coordinate axis. If set to a y
            axis id (e.g. "y" or "y2"), the `y` position refers to
            a y coordinate. If set to "paper", the `y` position
            refers to the distance from the bottom of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis.

        Returns
        -------
        Selection
        """
        super().__init__("selections")
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
The first argument to the plotly.graph_objs.layout.Selection
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Selection`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._set_property("name", arg, name)
        self._set_property("opacity", arg, opacity)
        self._set_property("path", arg, path)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("type", arg, type)
        self._set_property("x0", arg, x0)
        self._set_property("x1", arg, x1)
        self._set_property("xref", arg, xref)
        self._set_property("y0", arg, y0)
        self._set_property("y1", arg, y1)
        self._set_property("yref", arg, yref)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
