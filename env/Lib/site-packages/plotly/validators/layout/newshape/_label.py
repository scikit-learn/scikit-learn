import _plotly_utils.basevalidators


class LabelValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="label", parent_name="layout.newshape", **kwargs):
        super(LabelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Label"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets the new shape label text font.
            padding
                Sets padding (in px) between edge of label and
                edge of new shape.
            text
                Sets the text to display with the new shape. It
                is also used for legend item if `name` is not
                provided.
            textangle
                Sets the angle at which the label text is drawn
                with respect to the horizontal. For lines,
                angle "auto" is the same angle as the line. For
                all other shapes, angle "auto" is horizontal.
            textposition
                Sets the position of the label text relative to
                the new shape. Supported values for rectangles,
                circles and paths are *top left*, *top center*,
                *top right*, *middle left*, *middle center*,
                *middle right*, *bottom left*, *bottom center*,
                and *bottom right*. Supported values for lines
                are "start", "middle", and "end". Default:
                *middle center* for rectangles, circles, and
                paths; "middle" for lines.
            texttemplate
                Template string used for rendering the new
                shape's label. Note that this will override
                `text`. Variables are inserted using
                %{variable}, for example "x0: %{x0}". Numbers
                are formatted using d3-format's syntax
                %{variable:d3-format}, for example "Price:
                %{x0:$.2f}". See https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format for details on the
                formatting syntax. Dates are formatted using
                d3-time-format's syntax %{variable|d3-time-
                format}, for example "Day: %{x0|%m %b %Y}". See
                https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format for details on
                the date formatting syntax. A single
                multiplication or division operation may be
                applied to numeric variables, and combined with
                d3 number formatting, for example "Length in
                cm: %{x0*2.54}", "%{slope*60:.1f} meters per
                second." For log axes, variable values are
                given in log units. For date axes, x/y
                coordinate variables and center variables use
                datetimes, while all other variable values use
                values in ms. Finally, the template string has
                access to variables `x0`, `x1`, `y0`, `y1`,
                `slope`, `dx`, `dy`, `width`, `height`,
                `length`, `xcenter` and `ycenter`.
            xanchor
                Sets the label's horizontal position anchor
                This anchor binds the specified `textposition`
                to the "left", "center" or "right" of the label
                text. For example, if `textposition` is set to
                *top right* and `xanchor` to "right" then the
                right-most portion of the label text lines up
                with the right-most edge of the new shape.
            yanchor
                Sets the label's vertical position anchor This
                anchor binds the specified `textposition` to
                the "top", "middle" or "bottom" of the label
                text. For example, if `textposition` is set to
                *top right* and `yanchor` to "top" then the
                top-most portion of the label text lines up
                with the top-most edge of the new shape.
""",
            ),
            **kwargs,
        )
