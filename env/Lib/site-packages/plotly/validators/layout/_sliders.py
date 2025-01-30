import _plotly_utils.basevalidators


class SlidersValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="sliders", parent_name="layout", **kwargs):
        super(SlidersValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Slider"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            active
                Determines which button (by index starting from
                0) is considered active.
            activebgcolor
                Sets the background color of the slider grip
                while dragging.
            bgcolor
                Sets the background color of the slider.
            bordercolor
                Sets the color of the border enclosing the
                slider.
            borderwidth
                Sets the width (in px) of the border enclosing
                the slider.
            currentvalue
                :class:`plotly.graph_objects.layout.slider.Curr
                entvalue` instance or dict with compatible
                properties
            font
                Sets the font of the slider step labels.
            len
                Sets the length of the slider This measure
                excludes the padding of both ends. That is, the
                slider's length is this length minus the
                padding on both ends.
            lenmode
                Determines whether this slider length is set in
                units of plot "fraction" or in *pixels. Use
                `len` to set the value.
            minorticklen
                Sets the length in pixels of minor step tick
                marks
            name
                When used in a template, named items are
                created in the output figure in addition to any
                items the figure already has in this array. You
                can modify these items in the output figure by
                making your own item with `templateitemname`
                matching this `name` alongside your
                modifications (including `visible: false` or
                `enabled: false` to hide it). Has no effect
                outside of a template.
            pad
                Set the padding of the slider component along
                each side.
            steps
                A tuple of :class:`plotly.graph_objects.layout.
                slider.Step` instances or dicts with compatible
                properties
            stepdefaults
                When used in a template (as
                layout.template.layout.slider.stepdefaults),
                sets the default property values to use for
                elements of layout.slider.steps
            templateitemname
                Used to refer to a named item in this array in
                the template. Named items from the template
                will be created even without a matching item in
                the input figure, but you can modify one by
                making an item with `templateitemname` matching
                its `name`, alongside your modifications
                (including `visible: false` or `enabled: false`
                to hide it). If there is no template or no
                matching item, this item will be hidden unless
                you explicitly show it with `visible: true`.
            tickcolor
                Sets the color of the border enclosing the
                slider.
            ticklen
                Sets the length in pixels of step tick marks
            tickwidth
                Sets the tick width (in px).
            transition
                :class:`plotly.graph_objects.layout.slider.Tran
                sition` instance or dict with compatible
                properties
            visible
                Determines whether or not the slider is
                visible.
            x
                Sets the x position (in normalized coordinates)
                of the slider.
            xanchor
                Sets the slider's horizontal position anchor.
                This anchor binds the `x` position to the
                "left", "center" or "right" of the range
                selector.
            y
                Sets the y position (in normalized coordinates)
                of the slider.
            yanchor
                Sets the slider's vertical position anchor This
                anchor binds the `y` position to the "top",
                "middle" or "bottom" of the range selector.
""",
            ),
            **kwargs,
        )
