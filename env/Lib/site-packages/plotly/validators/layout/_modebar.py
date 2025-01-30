import _plotly_utils.basevalidators


class ModebarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="modebar", parent_name="layout", **kwargs):
        super(ModebarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Modebar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            activecolor
                Sets the color of the active or hovered on
                icons in the modebar.
            add
                Determines which predefined modebar buttons to
                add. Please note that these buttons will only
                be shown if they are compatible with all trace
                types used in a graph. Similar to
                `config.modeBarButtonsToAdd` option. This may
                include "v1hovermode", "hoverclosest",
                "hovercompare", "togglehover",
                "togglespikelines", "drawline", "drawopenpath",
                "drawclosedpath", "drawcircle", "drawrect",
                "eraseshape".
            addsrc
                Sets the source reference on Chart Studio Cloud
                for `add`.
            bgcolor
                Sets the background color of the modebar.
            color
                Sets the color of the icons in the modebar.
            orientation
                Sets the orientation of the modebar.
            remove
                Determines which predefined modebar buttons to
                remove. Similar to
                `config.modeBarButtonsToRemove` option. This
                may include "autoScale2d", "autoscale",
                "editInChartStudio", "editinchartstudio",
                "hoverCompareCartesian", "hovercompare",
                "lasso", "lasso2d", "orbitRotation",
                "orbitrotation", "pan", "pan2d", "pan3d",
                "reset", "resetCameraDefault3d",
                "resetCameraLastSave3d", "resetGeo",
                "resetSankeyGroup", "resetScale2d",
                "resetViewMap", "resetViewMapbox",
                "resetViews", "resetcameradefault",
                "resetcameralastsave", "resetsankeygroup",
                "resetscale", "resetview", "resetviews",
                "select", "select2d", "sendDataToCloud",
                "senddatatocloud", "tableRotation",
                "tablerotation", "toImage", "toggleHover",
                "toggleSpikelines", "togglehover",
                "togglespikelines", "toimage", "zoom",
                "zoom2d", "zoom3d", "zoomIn2d", "zoomInGeo",
                "zoomInMap", "zoomInMapbox", "zoomOut2d",
                "zoomOutGeo", "zoomOutMap", "zoomOutMapbox",
                "zoomin", "zoomout".
            removesrc
                Sets the source reference on Chart Studio Cloud
                for `remove`.
            uirevision
                Controls persistence of user-driven changes
                related to the modebar, including `hovermode`,
                `dragmode`, and `showspikes` at both the root
                level and inside subplots. Defaults to
                `layout.uirevision`.
""",
            ),
            **kwargs,
        )
