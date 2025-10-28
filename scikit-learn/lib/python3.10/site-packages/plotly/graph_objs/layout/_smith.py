#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Smith(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.smith"
    _valid_props = {"bgcolor", "domain", "imaginaryaxis", "realaxis"}

    @property
    def bgcolor(self):
        """
        Set the background color of the subplot

        The 'bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.smith.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

        Returns
        -------
        plotly.graph_objs.layout.smith.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    @property
    def imaginaryaxis(self):
        """
        The 'imaginaryaxis' property is an instance of Imaginaryaxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.smith.Imaginaryaxis`
          - A dict of string/value properties that will be passed
            to the Imaginaryaxis constructor

        Returns
        -------
        plotly.graph_objs.layout.smith.Imaginaryaxis
        """
        return self["imaginaryaxis"]

    @imaginaryaxis.setter
    def imaginaryaxis(self, val):
        self["imaginaryaxis"] = val

    @property
    def realaxis(self):
        """
        The 'realaxis' property is an instance of Realaxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.smith.Realaxis`
          - A dict of string/value properties that will be passed
            to the Realaxis constructor

        Returns
        -------
        plotly.graph_objs.layout.smith.Realaxis
        """
        return self["realaxis"]

    @realaxis.setter
    def realaxis(self, val):
        self["realaxis"] = val

    @property
    def _prop_descriptions(self):
        return """\
        bgcolor
            Set the background color of the subplot
        domain
            :class:`plotly.graph_objects.layout.smith.Domain`
            instance or dict with compatible properties
        imaginaryaxis
            :class:`plotly.graph_objects.layout.smith.Imaginaryaxis
            ` instance or dict with compatible properties
        realaxis
            :class:`plotly.graph_objects.layout.smith.Realaxis`
            instance or dict with compatible properties
        """

    def __init__(
        self,
        arg=None,
        bgcolor=None,
        domain=None,
        imaginaryaxis=None,
        realaxis=None,
        **kwargs,
    ):
        """
        Construct a new Smith object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Smith`
        bgcolor
            Set the background color of the subplot
        domain
            :class:`plotly.graph_objects.layout.smith.Domain`
            instance or dict with compatible properties
        imaginaryaxis
            :class:`plotly.graph_objects.layout.smith.Imaginaryaxis
            ` instance or dict with compatible properties
        realaxis
            :class:`plotly.graph_objects.layout.smith.Realaxis`
            instance or dict with compatible properties

        Returns
        -------
        Smith
        """
        super().__init__("smith")
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
The first argument to the plotly.graph_objs.layout.Smith
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Smith`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("domain", arg, domain)
        self._set_property("imaginaryaxis", arg, imaginaryaxis)
        self._set_property("realaxis", arg, realaxis)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
