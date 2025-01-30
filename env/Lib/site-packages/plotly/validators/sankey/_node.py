import _plotly_utils.basevalidators


class NodeValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="node", parent_name="sankey", **kwargs):
        super(NodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Node"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            align
                Sets the alignment method used to position the
                nodes along the horizontal axis.
            color
                Sets the `node` color. It can be a single
                value, or an array for specifying color for
                each `node`. If `node.color` is omitted, then
                the default `Plotly` color palette will be
                cycled through to have a variety of colors.
                These defaults are not fully opaque, to allow
                some visibility of what is beneath the node.
            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
            customdata
                Assigns extra data to each node.
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            groups
                Groups of nodes. Each group is defined by an
                array with the indices of the nodes it
                contains. Multiple groups can be specified.
            hoverinfo
                Determines which trace information appear when
                hovering nodes. If `none` or `skip` are set, no
                information is displayed upon hovering. But, if
                `none` is set, click and hover events are still
                fired.
            hoverlabel
                :class:`plotly.graph_objects.sankey.node.Hoverl
                abel` instance or dict with compatible
                properties
            hovertemplate
                Template string used for rendering the
                information that appear on hover box. Note that
                this will override `hoverinfo`. Variables are
                inserted using %{variable}, for example "y:
                %{y}" as well as %{xother}, {%_xother},
                {%_xother_}, {%xother_}. When showing info for
                several points, "xother" will be added to those
                with different x positions from the first
                point. An underscore before or after
                "(x|y)other" will add a space on that side,
                only when this field is shown. Numbers are
                formatted using d3-format's syntax
                %{variable:d3-format}, for example "Price:
                %{y:$.2f}". https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format for details on the
                formatting syntax. Dates are formatted using
                d3-time-format's syntax %{variable|d3-time-
                format}, for example "Day: %{2019-01-01|%A}".
                https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format for details on
                the date formatting syntax. The variables
                available in `hovertemplate` are the ones
                emitted as event data described at this link
                https://plotly.com/javascript/plotlyjs-
                events/#event-data. Additionally, every
                attributes that can be specified per-point (the
                ones that are `arrayOk: true`) are available.
                Variables `sourceLinks` and `targetLinks` are
                arrays of link objects.Finally, the template
                string has access to variables `value` and
                `label`. Anything contained in tag `<extra>` is
                displayed in the secondary box, for example
                "<extra>{fullData.name}</extra>". To hide the
                secondary box completely, use an empty tag
                `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
            label
                The shown name of the node.
            labelsrc
                Sets the source reference on Chart Studio Cloud
                for `label`.
            line
                :class:`plotly.graph_objects.sankey.node.Line`
                instance or dict with compatible properties
            pad
                Sets the padding (in px) between the `nodes`.
            thickness
                Sets the thickness (in px) of the `nodes`.
            x
                The normalized horizontal position of the node.
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            y
                The normalized vertical position of the node.
            ysrc
                Sets the source reference on Chart Studio Cloud
                for `y`.
""",
            ),
            **kwargs,
        )
