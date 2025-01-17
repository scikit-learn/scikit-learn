import _plotly_utils.basevalidators


class ButtonsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(
        self, plotly_name="buttons", parent_name="layout.xaxis.rangeselector", **kwargs
    ):
        super(ButtonsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Button"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            count
                Sets the number of steps to take to update the
                range. Use with `step` to specify the update
                interval.
            label
                Sets the text label to appear on the button.
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
            step
                The unit of measurement that the `count` value
                will set the range by.
            stepmode
                Sets the range update mode. If "backward", the
                range update shifts the start of range back
                "count" times "step" milliseconds. If "todate",
                the range update shifts the start of range back
                to the first timestamp from "count" times
                "step" milliseconds back. For example, with
                `step` set to "year" and `count` set to 1 the
                range update shifts the start of the range back
                to January 01 of the current year. Month and
                year "todate" are currently available only for
                the built-in (Gregorian) calendar.
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
            visible
                Determines whether or not this button is
                visible.
""",
            ),
            **kwargs,
        )
