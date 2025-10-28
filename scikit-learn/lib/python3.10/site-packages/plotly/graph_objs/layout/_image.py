#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Image(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.image"
    _valid_props = {
        "layer",
        "name",
        "opacity",
        "sizex",
        "sizey",
        "sizing",
        "source",
        "templateitemname",
        "visible",
        "x",
        "xanchor",
        "xref",
        "y",
        "yanchor",
        "yref",
    }

    @property
    def layer(self):
        """
        Specifies whether images are drawn below or above traces. When
        `xref` and `yref` are both set to `paper`, image is drawn below
        the entire plot area.

        The 'layer' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['below', 'above']

        Returns
        -------
        Any
        """
        return self["layer"]

    @layer.setter
    def layer(self, val):
        self["layer"] = val

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
        Sets the opacity of the image.

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
    def sizex(self):
        """
        Sets the image container size horizontally. The image will be
        sized based on the `position` value. When `xref` is set to
        `paper`, units are sized relative to the plot width. When
        `xref` ends with ` domain`, units are sized relative to the
        axis width.

        The 'sizex' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["sizex"]

    @sizex.setter
    def sizex(self, val):
        self["sizex"] = val

    @property
    def sizey(self):
        """
        Sets the image container size vertically. The image will be
        sized based on the `position` value. When `yref` is set to
        `paper`, units are sized relative to the plot height. When
        `yref` ends with ` domain`, units are sized relative to the
        axis height.

        The 'sizey' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["sizey"]

    @sizey.setter
    def sizey(self, val):
        self["sizey"] = val

    @property
    def sizing(self):
        """
        Specifies which dimension of the image to constrain.

        The 'sizing' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['fill', 'contain', 'stretch']

        Returns
        -------
        Any
        """
        return self["sizing"]

    @sizing.setter
    def sizing(self, val):
        self["sizing"] = val

    @property
    def source(self):
        """
        Specifies the URL of the image to be used. The URL must be
        accessible from the domain where the plot code is run, and can
        be either relative or absolute.

        The 'source' property is an image URI that may be specified as:
          - A remote image URI string
            (e.g. 'http://www.somewhere.com/image.png')
          - A data URI image string
            (e.g. 'data:image/png;base64,iVBORw0KGgoAAAANSU')
          - A PIL.Image.Image object which will be immediately converted
            to a data URI image string
            See http://pillow.readthedocs.io/en/latest/reference/Image.html

        Returns
        -------
        str
        """
        return self["source"]

    @source.setter
    def source(self, val):
        self["source"] = val

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
    def visible(self):
        """
        Determines whether or not this image is visible.

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    @property
    def x(self):
        """
        Sets the image's x position. When `xref` is set to `paper`,
        units are sized relative to the plot height. See `xref` for
        more info

        The 'x' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def xanchor(self):
        """
        Sets the anchor for the x position

        The 'xanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'center', 'right']

        Returns
        -------
        Any
        """
        return self["xanchor"]

    @xanchor.setter
    def xanchor(self, val):
        self["xanchor"] = val

    @property
    def xref(self):
        """
        Sets the images's x coordinate axis. If set to a x axis id
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
    def y(self):
        """
        Sets the image's y position. When `yref` is set to `paper`,
        units are sized relative to the plot height. See `yref` for
        more info

        The 'y' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def yanchor(self):
        """
        Sets the anchor for the y position.

        The 'yanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'middle', 'bottom']

        Returns
        -------
        Any
        """
        return self["yanchor"]

    @yanchor.setter
    def yanchor(self, val):
        self["yanchor"] = val

    @property
    def yref(self):
        """
        Sets the images's y coordinate axis. If set to a y axis id
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
        layer
            Specifies whether images are drawn below or above
            traces. When `xref` and `yref` are both set to `paper`,
            image is drawn below the entire plot area.
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
            Sets the opacity of the image.
        sizex
            Sets the image container size horizontally. The image
            will be sized based on the `position` value. When
            `xref` is set to `paper`, units are sized relative to
            the plot width. When `xref` ends with ` domain`, units
            are sized relative to the axis width.
        sizey
            Sets the image container size vertically. The image
            will be sized based on the `position` value. When
            `yref` is set to `paper`, units are sized relative to
            the plot height. When `yref` ends with ` domain`, units
            are sized relative to the axis height.
        sizing
            Specifies which dimension of the image to constrain.
        source
            Specifies the URL of the image to be used. The URL must
            be accessible from the domain where the plot code is
            run, and can be either relative or absolute.
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
        visible
            Determines whether or not this image is visible.
        x
            Sets the image's x position. When `xref` is set to
            `paper`, units are sized relative to the plot height.
            See `xref` for more info
        xanchor
            Sets the anchor for the x position
        xref
            Sets the images's x coordinate axis. If set to a x axis
            id (e.g. "x" or "x2"), the `x` position refers to a x
            coordinate. If set to "paper", the `x` position refers
            to the distance from the left of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        y
            Sets the image's y position. When `yref` is set to
            `paper`, units are sized relative to the plot height.
            See `yref` for more info
        yanchor
            Sets the anchor for the y position.
        yref
            Sets the images's y coordinate axis. If set to a y axis
            id (e.g. "y" or "y2"), the `y` position refers to a y
            coordinate. If set to "paper", the `y` position refers
            to the distance from the bottom of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            bottom (top). If set to a y axis ID followed by
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
        layer=None,
        name=None,
        opacity=None,
        sizex=None,
        sizey=None,
        sizing=None,
        source=None,
        templateitemname=None,
        visible=None,
        x=None,
        xanchor=None,
        xref=None,
        y=None,
        yanchor=None,
        yref=None,
        **kwargs,
    ):
        """
        Construct a new Image object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Image`
        layer
            Specifies whether images are drawn below or above
            traces. When `xref` and `yref` are both set to `paper`,
            image is drawn below the entire plot area.
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
            Sets the opacity of the image.
        sizex
            Sets the image container size horizontally. The image
            will be sized based on the `position` value. When
            `xref` is set to `paper`, units are sized relative to
            the plot width. When `xref` ends with ` domain`, units
            are sized relative to the axis width.
        sizey
            Sets the image container size vertically. The image
            will be sized based on the `position` value. When
            `yref` is set to `paper`, units are sized relative to
            the plot height. When `yref` ends with ` domain`, units
            are sized relative to the axis height.
        sizing
            Specifies which dimension of the image to constrain.
        source
            Specifies the URL of the image to be used. The URL must
            be accessible from the domain where the plot code is
            run, and can be either relative or absolute.
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
        visible
            Determines whether or not this image is visible.
        x
            Sets the image's x position. When `xref` is set to
            `paper`, units are sized relative to the plot height.
            See `xref` for more info
        xanchor
            Sets the anchor for the x position
        xref
            Sets the images's x coordinate axis. If set to a x axis
            id (e.g. "x" or "x2"), the `x` position refers to a x
            coordinate. If set to "paper", the `x` position refers
            to the distance from the left of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        y
            Sets the image's y position. When `yref` is set to
            `paper`, units are sized relative to the plot height.
            See `yref` for more info
        yanchor
            Sets the anchor for the y position.
        yref
            Sets the images's y coordinate axis. If set to a y axis
            id (e.g. "y" or "y2"), the `y` position refers to a y
            coordinate. If set to "paper", the `y` position refers
            to the distance from the bottom of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis.

        Returns
        -------
        Image
        """
        super().__init__("images")
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
The first argument to the plotly.graph_objs.layout.Image
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Image`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("layer", arg, layer)
        self._set_property("name", arg, name)
        self._set_property("opacity", arg, opacity)
        self._set_property("sizex", arg, sizex)
        self._set_property("sizey", arg, sizey)
        self._set_property("sizing", arg, sizing)
        self._set_property("source", arg, source)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("visible", arg, visible)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("xref", arg, xref)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._set_property("yref", arg, yref)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
