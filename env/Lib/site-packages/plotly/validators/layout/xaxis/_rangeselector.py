import _plotly_utils.basevalidators


class RangeselectorValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="rangeselector", parent_name="layout.xaxis", **kwargs
    ):
        super(RangeselectorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Rangeselector"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            activecolor
                Sets the background color of the active range
                selector button.
            bgcolor
                Sets the background color of the range selector
                buttons.
            bordercolor
                Sets the color of the border enclosing the
                range selector.
            borderwidth
                Sets the width (in px) of the border enclosing
                the range selector.
            buttons
                Sets the specifications for each buttons. By
                default, a range selector comes with no
                buttons.
            buttondefaults
                When used in a template (as layout.template.lay
                out.xaxis.rangeselector.buttondefaults), sets
                the default property values to use for elements
                of layout.xaxis.rangeselector.buttons
            font
                Sets the font of the range selector button
                text.
            visible
                Determines whether or not this range selector
                is visible. Note that range selectors are only
                available for x axes of `type` set to or auto-
                typed to "date".
            x
                Sets the x position (in normalized coordinates)
                of the range selector.
            xanchor
                Sets the range selector's horizontal position
                anchor. This anchor binds the `x` position to
                the "left", "center" or "right" of the range
                selector.
            y
                Sets the y position (in normalized coordinates)
                of the range selector.
            yanchor
                Sets the range selector's vertical position
                anchor This anchor binds the `y` position to
                the "top", "middle" or "bottom" of the range
                selector.
""",
            ),
            **kwargs,
        )
