import _plotly_utils.basevalidators


class StepsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="steps", parent_name="layout.slider", **kwargs):
        super(StepsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Step"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            args
                Sets the arguments values to be passed to the
                Plotly method set in `method` on slide.
            execute
                When true, the API method is executed. When
                false, all other behaviors are the same and
                command execution is skipped. This may be
                useful when hooking into, for example, the
                `plotly_sliderchange` method and executing the
                API command manually without losing the benefit
                of the slider automatically binding to the
                state of the plot through the specification of
                `method` and `args`.
            label
                Sets the text label to appear on the slider
            method
                Sets the Plotly method to be called when the
                slider value is changed. If the `skip` method
                is used, the API slider will function as normal
                but will perform no API calls and will not bind
                automatically to state updates. This may be
                used to create a component interface and attach
                to slider events manually via JavaScript.
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
            value
                Sets the value of the slider step, used to
                refer to the step programatically. Defaults to
                the slider label if not provided.
            visible
                Determines whether or not this step is included
                in the slider.
""",
            ),
            **kwargs,
        )
