import _plotly_utils.basevalidators


class TransitionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="transition", parent_name="layout", **kwargs):
        super(TransitionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Transition"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            duration
                The duration of the transition, in
                milliseconds. If equal to zero, updates are
                synchronous.
            easing
                The easing function used for the transition
            ordering
                Determines whether the figure's layout or
                traces smoothly transitions during updates that
                make both traces and layout change.
""",
            ),
            **kwargs,
        )
