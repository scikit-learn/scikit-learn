import _plotly_utils.basevalidators


class ParcatsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="parcats", parent_name="", **kwargs):
        super(ParcatsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Parcats"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            arrangement
                Sets the drag interaction mode for categories
                and dimensions. If `perpendicular`, the
                categories can only move along a line
                perpendicular to the paths. If `freeform`, the
                categories can freely move on the plane. If
                `fixed`, the categories and dimensions are
                stationary.
            bundlecolors
                Sort paths so that like colors are bundled
                together within each category.
            counts
                The number of observations represented by each
                state. Defaults to 1 so that each state
                represents one observation
            countssrc
                Sets the source reference on Chart Studio Cloud
                for `counts`.
            dimensions
                The dimensions (variables) of the parallel
                categories diagram.
            dimensiondefaults
                When used in a template (as layout.template.dat
                a.parcats.dimensiondefaults), sets the default
                property values to use for elements of
                parcats.dimensions
            domain
                :class:`plotly.graph_objects.parcats.Domain`
                instance or dict with compatible properties
            hoverinfo
                Determines which trace information appear on
                hover. If `none` or `skip` are set, no
                information is displayed upon hovering. But, if
                `none` is set, click and hover events are still
                fired.
            hoveron
                Sets the hover interaction mode for the parcats
                diagram. If `category`, hover interaction take
                place per category. If `color`, hover
                interactions take place per color per category.
                If `dimension`, hover interactions take place
                across all categories per dimension.
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
                This value here applies when hovering over
                dimensions. Note that `*categorycount`,
                "colorcount" and "bandcolorcount" are only
                available when `hoveron` contains the "color"
                flagFinally, the template string has access to
                variables `count`, `probability`, `category`,
                `categorycount`, `colorcount` and
                `bandcolorcount`. Anything contained in tag
                `<extra>` is displayed in the secondary box,
                for example "<extra>{fullData.name}</extra>".
                To hide the secondary box completely, use an
                empty tag `<extra></extra>`.
            labelfont
                Sets the font for the `dimension` labels.
            legendgrouptitle
                :class:`plotly.graph_objects.parcats.Legendgrou
                ptitle` instance or dict with compatible
                properties
            legendwidth
                Sets the width (in px or fraction) of the
                legend for this trace.
            line
                :class:`plotly.graph_objects.parcats.Line`
                instance or dict with compatible properties
            meta
                Assigns extra meta information associated with
                this trace that can be used in various text
                attributes. Attributes such as trace `name`,
                graph, axis and colorbar `title.text`,
                annotation `text` `rangeselector`,
                `updatemenues` and `sliders` `label` text all
                support `meta`. To access the trace `meta`
                values in an attribute in the same trace,
                simply use `%{meta[i]}` where `i` is the index
                or key of the `meta` item in question. To
                access trace `meta` in layout attributes, use
                `%{data[n[.meta[i]}` where `i` is the index or
                key of the `meta` and `n` is the trace index.
            metasrc
                Sets the source reference on Chart Studio Cloud
                for `meta`.
            name
                Sets the trace name. The trace name appears as
                the legend item and on hover.
            sortpaths
                Sets the path sorting algorithm. If `forward`,
                sort paths based on dimension categories from
                left to right. If `backward`, sort paths based
                on dimensions categories from right to left.
            stream
                :class:`plotly.graph_objects.parcats.Stream`
                instance or dict with compatible properties
            tickfont
                Sets the font for the `category` labels.
            uid
                Assign an id to this trace, Use this to provide
                object constancy between traces during
                animations and transitions.
            uirevision
                Controls persistence of some user-driven
                changes to the trace: `constraintrange` in
                `parcoords` traces, as well as some `editable:
                true` modifications such as `name` and
                `colorbar.title`. Defaults to
                `layout.uirevision`. Note that other user-
                driven trace attribute changes are controlled
                by `layout` attributes: `trace.visible` is
                controlled by `layout.legend.uirevision`,
                `selectedpoints` is controlled by
                `layout.selectionrevision`, and
                `colorbar.(x|y)` (accessible with `config:
                {editable: true}`) is controlled by
                `layout.editrevision`. Trace changes are
                tracked by `uid`, which only falls back on
                trace index if no `uid` is provided. So if your
                app can add/remove traces before the end of the
                `data` array, such that the same trace has a
                different index, you can still preserve user-
                driven changes if you give each trace a `uid`
                that stays with it as it moves.
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
""",
            ),
            **kwargs,
        )
