import _plotly_utils.basevalidators


class Mesh3DValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="mesh3d", parent_name="", **kwargs):
        super(Mesh3DValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Mesh3d"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            alphahull
                Determines how the mesh surface triangles are
                derived from the set of vertices (points)
                represented by the `x`, `y` and `z` arrays, if
                the `i`, `j`, `k` arrays are not supplied. For
                general use of `mesh3d` it is preferred that
                `i`, `j`, `k` are supplied. If "-1", Delaunay
                triangulation is used, which is mainly suitable
                if the mesh is a single, more or less layer
                surface that is perpendicular to
                `delaunayaxis`. In case the `delaunayaxis`
                intersects the mesh surface at more than one
                point it will result triangles that are very
                long in the dimension of `delaunayaxis`. If
                ">0", the alpha-shape algorithm is used. In
                this case, the positive `alphahull` value
                signals the use of the alpha-shape algorithm,
                _and_ its value acts as the parameter for the
                mesh fitting. If 0,  the convex-hull algorithm
                is used. It is suitable for convex bodies or if
                the intention is to enclose the `x`, `y` and
                `z` point set into a convex hull.
            autocolorscale
                Determines whether the colorscale is a default
                palette (`autocolorscale: true`) or the palette
                determined by `colorscale`. In case
                `colorscale` is unspecified or `autocolorscale`
                is true, the default palette will be chosen
                according to whether numbers in the `color`
                array are all positive, all negative or mixed.
            cauto
                Determines whether or not the color domain is
                computed with respect to the input data (here
                `intensity`) or the bounds set in `cmin` and
                `cmax` Defaults to `false` when `cmin` and
                `cmax` are set by the user.
            cmax
                Sets the upper bound of the color domain. Value
                should have the same units as `intensity` and
                if set, `cmin` must be set as well.
            cmid
                Sets the mid-point of the color domain by
                scaling `cmin` and/or `cmax` to be equidistant
                to this point. Value should have the same units
                as `intensity`. Has no effect when `cauto` is
                `false`.
            cmin
                Sets the lower bound of the color domain. Value
                should have the same units as `intensity` and
                if set, `cmax` must be set as well.
            color
                Sets the color of the whole mesh
            coloraxis
                Sets a reference to a shared color axis.
                References to these shared color axes are
                "coloraxis", "coloraxis2", "coloraxis3", etc.
                Settings for these shared color axes are set in
                the layout, under `layout.coloraxis`,
                `layout.coloraxis2`, etc. Note that multiple
                color scales can be linked to the same color
                axis.
            colorbar
                :class:`plotly.graph_objects.mesh3d.ColorBar`
                instance or dict with compatible properties
            colorscale
                Sets the colorscale. The colorscale must be an
                array containing arrays mapping a normalized
                value to an rgb, rgba, hex, hsl, hsv, or named
                color string. At minimum, a mapping for the
                lowest (0) and highest (1) values are required.
                For example, `[[0, 'rgb(0,0,255)'], [1,
                'rgb(255,0,0)']]`. To control the bounds of the
                colorscale in color space, use `cmin` and
                `cmax`. Alternatively, `colorscale` may be a
                palette name string of the following list: Blac
                kbody,Bluered,Blues,Cividis,Earth,Electric,Gree
                ns,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,R
                eds,Viridis,YlGnBu,YlOrRd.
            contour
                :class:`plotly.graph_objects.mesh3d.Contour`
                instance or dict with compatible properties
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            delaunayaxis
                Sets the Delaunay axis, which is the axis that
                is perpendicular to the surface of the Delaunay
                triangulation. It has an effect if `i`, `j`,
                `k` are not provided and `alphahull` is set to
                indicate Delaunay triangulation.
            facecolor
                Sets the color of each face Overrides "color"
                and "vertexcolor".
            facecolorsrc
                Sets the source reference on Chart Studio Cloud
                for `facecolor`.
            flatshading
                Determines whether or not normal smoothing is
                applied to the meshes, creating meshes with an
                angular, low-poly look via flat reflections.
            hoverinfo
                Determines which trace information appear on
                hover. If `none` or `skip` are set, no
                information is displayed upon hovering. But, if
                `none` is set, click and hover events are still
                fired.
            hoverinfosrc
                Sets the source reference on Chart Studio Cloud
                for `hoverinfo`.
            hoverlabel
                :class:`plotly.graph_objects.mesh3d.Hoverlabel`
                instance or dict with compatible properties
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
                Anything contained in tag `<extra>` is
                displayed in the secondary box, for example
                "<extra>{fullData.name}</extra>". To hide the
                secondary box completely, use an empty tag
                `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
            hovertext
                Same as `text`.
            hovertextsrc
                Sets the source reference on Chart Studio Cloud
                for `hovertext`.
            i
                A vector of vertex indices, i.e. integer values
                between 0 and the length of the vertex vectors,
                representing the "first" vertex of a triangle.
                For example, `{i[m], j[m], k[m]}` together
                represent face m (triangle m) in the mesh,
                where `i[m] = n` points to the triplet `{x[n],
                y[n], z[n]}` in the vertex arrays. Therefore,
                each element in `i` represents a point in
                space, which is the first vertex of a triangle.
            ids
                Assigns id labels to each datum. These ids for
                object constancy of data points during
                animation. Should be an array of strings, not
                numbers or any other type.
            idssrc
                Sets the source reference on Chart Studio Cloud
                for `ids`.
            intensity
                Sets the intensity values for vertices or cells
                as defined by `intensitymode`. It can be used
                for plotting fields on meshes.
            intensitymode
                Determines the source of `intensity` values.
            intensitysrc
                Sets the source reference on Chart Studio Cloud
                for `intensity`.
            isrc
                Sets the source reference on Chart Studio Cloud
                for `i`.
            j
                A vector of vertex indices, i.e. integer values
                between 0 and the length of the vertex vectors,
                representing the "second" vertex of a triangle.
                For example, `{i[m], j[m], k[m]}`  together
                represent face m (triangle m) in the mesh,
                where `j[m] = n` points to the triplet `{x[n],
                y[n], z[n]}` in the vertex arrays. Therefore,
                each element in `j` represents a point in
                space, which is the second vertex of a
                triangle.
            jsrc
                Sets the source reference on Chart Studio Cloud
                for `j`.
            k
                A vector of vertex indices, i.e. integer values
                between 0 and the length of the vertex vectors,
                representing the "third" vertex of a triangle.
                For example, `{i[m], j[m], k[m]}` together
                represent face m (triangle m) in the mesh,
                where `k[m] = n` points to the triplet  `{x[n],
                y[n], z[n]}` in the vertex arrays. Therefore,
                each element in `k` represents a point in
                space, which is the third vertex of a triangle.
            ksrc
                Sets the source reference on Chart Studio Cloud
                for `k`.
            legend
                Sets the reference to a legend to show this
                trace in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgroup
                Sets the legend group for this trace. Traces
                and shapes part of the same legend group
                hide/show at the same time when toggling legend
                items.
            legendgrouptitle
                :class:`plotly.graph_objects.mesh3d.Legendgroup
                title` instance or dict with compatible
                properties
            legendrank
                Sets the legend rank for this trace. Items and
                groups with smaller ranks are presented on
                top/left side while with "reversed"
                `legend.traceorder` they are on bottom/right
                side. The default legendrank is 1000, so that
                you can use ranks less than 1000 to place
                certain items before all unranked items, and
                ranks greater than 1000 to go after all
                unranked items. When having unranked or equal
                rank items shapes would be displayed after
                traces i.e. according to their order in data
                and layout.
            legendwidth
                Sets the width (in px or fraction) of the
                legend for this trace.
            lighting
                :class:`plotly.graph_objects.mesh3d.Lighting`
                instance or dict with compatible properties
            lightposition
                :class:`plotly.graph_objects.mesh3d.Lightpositi
                on` instance or dict with compatible properties
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
            opacity
                Sets the opacity of the surface. Please note
                that in the case of using high `opacity` values
                for example a value greater than or equal to
                0.5 on two surfaces (and 0.25 with four
                surfaces), an overlay of multiple transparent
                surfaces may not perfectly be sorted in depth
                by the webgl API. This behavior may be improved
                in the near future and is subject to change.
            reversescale
                Reverses the color mapping if true. If true,
                `cmin` will correspond to the last color in the
                array and `cmax` will correspond to the first
                color.
            scene
                Sets a reference between this trace's 3D
                coordinate system and a 3D scene. If "scene"
                (the default value), the (x,y,z) coordinates
                refer to `layout.scene`. If "scene2", the
                (x,y,z) coordinates refer to `layout.scene2`,
                and so on.
            showlegend
                Determines whether or not an item corresponding
                to this trace is shown in the legend.
            showscale
                Determines whether or not a colorbar is
                displayed for this trace.
            stream
                :class:`plotly.graph_objects.mesh3d.Stream`
                instance or dict with compatible properties
            text
                Sets the text elements associated with the
                vertices. If trace `hoverinfo` contains a
                "text" flag and "hovertext" is not set, these
                elements will be seen in the hover labels.
            textsrc
                Sets the source reference on Chart Studio Cloud
                for `text`.
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
            vertexcolor
                Sets the color of each vertex Overrides
                "color". While Red, green and blue colors are
                in the range of 0 and 255; in the case of
                having vertex color data in RGBA format, the
                alpha color should be normalized to be between
                0 and 1.
            vertexcolorsrc
                Sets the source reference on Chart Studio Cloud
                for `vertexcolor`.
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
            x
                Sets the X coordinates of the vertices. The nth
                element of vectors `x`, `y` and `z` jointly
                represent the X, Y and Z coordinates of the nth
                vertex.
            xcalendar
                Sets the calendar system to use with `x` date
                data.
            xhoverformat
                Sets the hover text formatting rulefor `x`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format. And for dates
                see: https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                *09~15~23.46*By default the values are
                formatted using `xaxis.hoverformat`.
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            y
                Sets the Y coordinates of the vertices. The nth
                element of vectors `x`, `y` and `z` jointly
                represent the X, Y and Z coordinates of the nth
                vertex.
            ycalendar
                Sets the calendar system to use with `y` date
                data.
            yhoverformat
                Sets the hover text formatting rulefor `y`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format. And for dates
                see: https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                *09~15~23.46*By default the values are
                formatted using `yaxis.hoverformat`.
            ysrc
                Sets the source reference on Chart Studio Cloud
                for `y`.
            z
                Sets the Z coordinates of the vertices. The nth
                element of vectors `x`, `y` and `z` jointly
                represent the X, Y and Z coordinates of the nth
                vertex.
            zcalendar
                Sets the calendar system to use with `z` date
                data.
            zhoverformat
                Sets the hover text formatting rulefor `z`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format. And for dates
                see: https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                *09~15~23.46*By default the values are
                formatted using `zaxis.hoverformat`.
            zsrc
                Sets the source reference on Chart Studio Cloud
                for `z`.
""",
            ),
            **kwargs,
        )
